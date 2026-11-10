# P7 计划：Trio binning 前端 + Hi-C 前端

## Context

rshifiasm 的默认 HiFi 流水线（P0–P6）已逐字节一致。P7 的目标是把 trio 分区和 Hi-C 分区两个前端补齐，使 `-1..-4`（trio）与 `--h1/--h2`（Hi-C）达到与 hifiasm 0.25.0-r726 相同的逐字节输出。深挖结论：

- **Trio 的下游机制已全部移植**（`src/assembly/trio_clean.rs` 29 个函数、`bp_graph.rs` 的 `output_trio_graph_joint`/`adjust_utg_by_trio`、`purge_dups.rs` 的 `mc_solve_edgeless`）。缺口只有：binning 前端（名单模式 + yak dump 模式）、`drop_edges_by_trio`、读图层 trio 切割器（`*_trio` 变体 + `is_trio` 分支，现为 `gfa_clean.rs:810/1121` 的 assert）、checkpoint。估算 ~3k 行 Rust，**确定性风险极低**（kt_for 线程无关、逐读 f64）。
- **Hi-C 是 6 阶段工程**（live ~10.9k 行 C++ → 估 ~13–14k 行 Rust，新 `src/hic/` 模块）。inter.cpp 链化、kalloc、线程池都**不需要**；horder.cpp 只需 misjoin 簇（~1200 行），scaffolding 是死代码不移植。唯一硬风险：Phase 5 的 x87 `long double` + `logl` 权重管线。
- **测试数据**：合成。门禁本质是 Rust vs C++ 同输入逐字节比对，数据只需覆盖执行路径。trio 名单用读名哈希拆分 chr11-2M；yak dump 用拆分后的伪单倍型读段生成（克隆编译 lh3/yak，失败则按已摸清的 dump 格式写生成器：magic `"YAK\2"` + `t[3]={k,pre,10}` + 每子表 `u32 cap,size` + `size×u64` 键，键低 10 位为计数）；Hi-C 用确定性脚本从 chr11-2M 模拟配对读段（顺式偏置的邻位对，~2×150bp，几十万对）。生成脚本入 `scripts/`（提交），数据落 `reference/`（git-ignored）。

实施顺序：**先 Trio（A→B→C），后 Hi-C（Phase 0→6）**。每阶段一个 commit，附专属字节门 + 既有 9-profile 默认门禁回归。

## Part 1 — Trio binning（三阶段）

### Stage A：名单模式 + 分发接线（`-3/-4`，无需外部数据）
- 前端：`ha_triobin`/`ha_triobin_list`/`ha_triobin_set_list`（Trio.cpp:450–457/328–350/294–326）→ 写 `ReadStore::trio_labels`（字段已存在，read_store.rs:25，全管线已穿线）。调用点对应 Assembly.cpp:2113：最终 overlap 轮之后、`build_string_graph_without_clean` 之前（`src/assembly/pipeline.rs`）。名单名→id 哈希保持"首现优先"。日志三行逐字复刻。
- `drop_edges_by_trio`（Overlaps.cpp:1590–1621）：插到 `normalize_ma_hit_t_single_side_advance` 之后（pipeline.rs ~160，注释已预留 else 半支）。
- 读图切割器：`get_tip_trio_infor`（Overlaps.cpp:6890，0.9/0.1 阈值常量已在 trio_clean.rs:42/48）+ `asg_arc_del_short_diploid_by_exact_trio`（:8241）+ `_by_length_trio`（:6937）；`gfa_clean.rs` 的 `asg_arc_cut_inexact`/`asg_arc_cut_length` 增补 `is_trio` 分支（gfa_ut.cpp:1007–1041/1326–1366，量小）并拆除 assert；`ul_clean_gfa`（pipeline.rs:355）加 `is_trio` 参数。
- 分发修复（pipeline.rs:286–292 的过时 Err）：`trio_binning_requested()` 且非 Hi-C 且非 `--trio-dual` → `output_trio_graph_joint`（bp_graph.rs:1438），参数照抄 Overlaps.cpp:39636–39638。
- **门 A**：C++ vs Rust `-3 l1 -4 l2 -t1`：三条 binning 日志行、`trio_labels` 探针 dump（扩展 `probe.rs` 加 `d00_after_triobin` 标签）、受 trio 影响的 d01/d02 阶段 dump、`*.hap1/hap2.p_ctg.gfa`（+noseq/lowQ.bed）与 `bp.p_ctg.gfa` 逐字节一致。

### Stage B：yak dump 模式（`-1/-2`）
- dump 读取器（`yak_ch_restore_core`，Trio.cpp:66–137：`cnt=key&0x3FF`，label=`cnt>=mid?2:(cnt>=min?1:drop)`，父本 shift 0 / 母本 shift 2 合入同一 khashl 表）+ `yak_hash64`/`yak_hash_long`（htab.h:138–166，平凡）+ `tb_worker`/`tb_classify`（Trio.cpp:173–259，`ratio_thres=0.33`）+ `ha_triobin_yak`（:261–292）。`-c/-d` 是装载期阈值；`--t-occ` 属下游（已移植）。
- 数据：`scripts/gen_trio_corpus.sh`（读名哈希拆两个伪单倍型 → 名单 + yak dump），先验证 C++ 能读我们生成的 dump。
- **门 B**：`-1 f.yak -2 m.yak -t1`：yak 日志行 + trio_labels 探针（应与名单真值一致）+ 下游产物同门 A；错误路径（`ERROR: wrong file magic.`）。

### Stage C：checkpoint + `--bin-only`
- `test_yak_binning`/`write_yak_binning`/`load_yak_binning`（Overlaps.cpp:23501–23565）：纯文本，首行 `name\t<yak路径>` 指纹；命中后走名单路径重放。接线到现有 `read_cache.rs` 的写/读流程；补 `binning_only` 的两个退出点。
- **门 C**：`.hap1/.hap2.phase.bin` 逐字节（相同 argv）+ 断点续跑产物一致。
- **保持显式 Err**：`--trio-dual`、`-5` 多倍体、`benchmark_hic_graph`（随 Hi-C trio 后置）。

## Part 2 — Hi-C（六阶段，新 `src/hic/` 模块）

- **Phase 0 胶水**（~150 行 C++）：接管 `output_hic_graph` 分发 Err（pipeline.rs）；`probe.rs` 加 Hi-C 阶段 dump；确认默认门禁不回归。
- **Phase 1 单元格索引**（hic.cpp:306–860 → `src/hic/ug_index.rs`，~800 行）：`ha_ug_index`——unitig 序列全精确 k-mer（k=31，`hc_hash_long` 128-bit 滚动规范化），2^12 桶、键低 12 位计数（`HIC_MAX_COUNT=4095`）、位置打包 `(rev<<63)|(uID<<...)|pos`、`hap_occ` 过滤；复用 `container/khash.rs` 的 `U64HashTable::with_ignored_low_bits` 与 `container/ksort.rs` radix（勿改 `PositionIndex`，语义不同）。**门 A：`.hic.tlb.bin` 逐字节**（自带图指纹）。
- **Phase 2 Hi-C 读比对**（hic.cpp:867–2205 live 子集 → `src/hic/map.rs`，~1800 行）：`get_alignment` k-mer 投票 + `set_pe_pos` R1/R2 合票（kt_pipeline 顺序化）+ `dedup_hits`/`sort_hits` + IO。**门 B：`.hic.lk.bin` 逐字节**。
- **Phase 3 气泡/转录/距离**（~2900 行 C++）：`identify_bubbles`（hic.cpp:2354–2607）、`resolve_tangles_hic`（:16259–16316）、`idx_hc_links`（:15321）、距离簇（:3387–5018：`measure_distance`/`all_pair_shortest_path`/`fill_utg_distance_multi`/`update_containment_distance`）、`refine_hic_trans`（Overlaps.cpp:15581–16231）；在 Hi-C 语料上重验 `hic_clean_adv` 跳过是否成立（不成立则移植）。**门 C：refine 后 `k_trans` 探针 dump**。
- **Phase 4 mc_solve**（rcut.cpp ~3500 行 → `src/hic/maxcut.rs`）：splitmix64/`kr_drand_r` RNG 逐字移植、radix 键打包（含 tie-break 索引）、贪心初始化的 bug 兼容分支（rcut.cpp:1565–1587 注释 "may have bugs"）、扰动循环与 `n_iter%flush` 复位。权重为整数值 double，f64 精确。**门 D：每轮 `label_unitigs_sm` 后 spin/trio_flag 探针**。
- **Phase 5 重加权 + 驱动**（~2500 行 C++）：`get_trans_rate_function_advance`、`LeastSquare_advance`、`get_trans_weight_advance`（`logl`）——**风险决策点**：先在 Hi-C 语料上探针对比 f64 与 C++ 分区结果；若翻转则仅对这两个函数实现 x87-80 软浮点（add/sub/mul/div + fyl2x 语义 log）。再加 `weight_kv_u_trans`/`renew_kv_u_trans`、`hic_short_align`（hic.cpp:17016–17164）、`hic_analysis`、`label_unitigs_sm`、misjoin 簇（horder.cpp:435–2084 `update_switch_unitig`，因 `--l-msjoin` 默认开）。**门 E：每轮 spin + 终态 per-read trio_flag**。
- **Phase 6 端到端**（~200 行胶水）：d_edges 删除、`print_utg`/`output_trio_graph_joint` 接线、`{o}.hic.*` 命名。**门 F：`{o}.hic.p_ctg.gfa`、`{o}.hic.hap1/hap2.p_ctg.gfa`（+noseq/bed）与两个 checkpoint 全量逐字节**。
- **后置**（保持 Err）：多倍体 `output_hic_graph_mmhap`/`mc_solve_general`、trio+Hi-C `benchmark_hic_graph`、`--dual-scaf`、`bp_solve`、死代码 H_partition/G_partition。

## 关键文件

- 新：`src/hic/`（ug_index.rs / map.rs / distance.rs / maxcut.rs / driver.rs）、`scripts/gen_trio_corpus.sh`、`scripts/gen_hic_corpus.py`、`src/command_line/trio_bin.rs`（或并入 trio_clean.rs）
- 改：`src/assembly/pipeline.rs`（triobin 调用点、drop_edges_by_trio、分发）、`src/assembly/gfa_clean.rs`（is_trio 分支）、`src/assembly/read_cache.rs`（phase.bin + bin_only）、`src/assembly/probe.rs`（新探针标签）、`src/reads/read_store.rs`（如需）、`src/command_line/parser.rs`（拆除 trio/Hi-C 的 parse 期拒绝——`check_supported` 现不拦 trio，Hi-C 拦截在分发处）、`docs/consistency.md`
- 复用：`trio_clean.rs` 全部、`bp_graph.rs:1438` `output_trio_graph_joint`、`unitig.rs` hic_clean/ma_ug_gen/ma_ug_print_bed、`container/khash.rs`、`container/ksort.rs`、`counting/minimizer.rs` 哈希

## 验证

每阶段 commit：(1) `cargo test` 全绿；(2) 该阶段专属字节门（上列各门，C++ 基线先跑入 `reference/`，`scripts/compare_outputs.sh` 比对，`ec.bin` 豁免惯例不变）；(3) 9-profile 默认门禁回归（default/write_paf_ec/primary/bcov/telo/dbg/e/i/load 全部 SHA-256 不变）；(4) 门 A/C 之前的阶段同时验证 C++ 接受合成输入（yak dump 魔数/格式、Hi-C fastq）。
