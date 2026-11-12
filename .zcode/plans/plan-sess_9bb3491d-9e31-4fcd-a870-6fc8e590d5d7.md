# hifiasm 迁移收尾：U7 e2e 门禁 + Hi-C 完整迁移

## Context

rshifiasm 是 hifiasm 0.25.0-r726 的逐字节一致 Rust 移植（C++ 只读参考在 `/home/keith/projects/hifiasm`）。主线（P0–P8 默认 HiFi、trio 全分支、ONT、`--sec-in`）已完成并门禁通过；`--ul` 的 U1–U6 机制层全部移植并接线（checkpoint 续十四，HEAD `2f0b613`），解析解锁已提交（`1085527`）。

**剩余工作**（用户要求一次性完成、不中断）：
1. **U7 e2e 门禁**：--ul 全链首次端到端执行，与 `/tmp/ulref` 参考逐字节对拍。
2. **Hi-C（`--h1`/`--h2`）H0–H7**：完全未开始；涉及 hic.cpp 活代码 ~6k 行 + Overlaps.cpp 驱动 ~800 行 + horder.cpp misjoin/scaffolding ~3.4k 行。

语料现状（已核实存在）：`/tmp/ul.fq`（89MB，739 条合成 UL 读）、`/tmp/ulref.*` 全套参考产物（11-10 生成）、`target/release/rshifiasm`（11-12 构建）、`/home/keith/Downloads/chr11-2M.fa.gz`。

## 方法（沿用既有流程，每阶段循环）

字面转录 C++ → `/tmp` 探针（`g++ -O2 -I /home/keith/projects/hifiasm xxx_probe.cpp $(ls /home/keith/projects/hifiasm/*.o | grep -v '/main.o$') -lz -lpthread`）生成黄金值 → Rust 回放单测（`#[ignored]`，关键黄金内嵌常驻）→ 接线 → 门禁（`scripts/compare_outputs.sh` + `cmp` + 日志锚点）→ 更新 `docs/consistency.md` → git commit（沿用 per-feature "Port ..." / "docs: checkpoint ..." 交替风格）。

已知豁免沿用：`.ec.bin` padding、`*.ul.ovlp.bin`/`*.re.ul.ovlp.bin` 未初始化字段、日志时间/CPU/内存字段。

## Part A：U7 e2e 门禁（完成 --ul）

1. `cargo build --release` → `target/release/rshifiasm -t1 --ul /tmp/ul.fq -o /tmp/ulrsm /home/keith/Downloads/chr11-2M.fa.gz`
2. 对拍：8 个 GFA/BED（`/tmp/ulref.bp{,.hap1,.hap2}.{p_ctg,p_utg,r_utg}.{gfa,noseq.gfa,lowQ.bed}`）逐字节；`re.ul.msk.bin` 逐字节；`re.uidx.bin`/`re.ul.ovlp.bin` 按豁免规则字段级对拍；日志锚点 vs `/tmp/ulref_run.log`（`rescall_ul_pipeline`/`print_integert_ovlp_stat`/`realloc_rdb_adv`/`update_paf`/`renew_ng`）。
3. 预期需调试（checkpoint 预判：错误多半在 `renew_ng`/`gen_ul2ul` 断言区）。用探针/阶段 dump 二分定位，修复后重跑至全对。
4. 回归：`cargo test --release` 全绿（345+13 ignored）；默认 9-profile、trio、`-5`、`--dual-scaf`、ONT 门禁无回归。
5. 更新 consistency.md（续十五）+ commit。

## Part B：Hi-C H0–H7

**H0 语料**：新 `scripts/gen_hic_reads.py`——从 chr11-2M 默认组装 contig 确定性合成 Hi-C 双端读对（2×150bp，接触衰减距离分布 + 随机碰撞背景，固定种子）；跑 C++ `hifiasm -t1 --h1 r1 --h2 r2 -o reference/hic chr11-2M.fa.gz` 固定参考（含 `.hic.tlb.bin`/`.hic.lk.bin`）。

**H1 hic_clean_adv**（Overlaps.cpp:14403 + 依赖的 hic.cpp 泡泡机制 `cal_bub_best_by_len` 等）：移植并填入 `pipeline.rs:373` 现跳过点与 `ul_renew_driver.rs` 的 no-op 桩（gradually_renew_g 与 HiFi-only 尾部都调用它）；`deep_graph_clean`（gfa_ut.cpp:18184）随之补齐。默认语料 no-op 验证 + 逐字段对拍。

**H2 hic_pt 单位图 k-mer 索引**（hic.cpp 349–917）：`hc_pt_t`（12 位计数器打包键+位置表）、`hc_pt_t_gen_single`/`worker_for`/`parallel_count_hc_pt1`、`build_unitig_index`(815)、`write/load_hc_pt_index`(`.hic.tlb.bin`)、`interpret_pos`/`check_exact_match`。黄金：索引内容对拍 C++ dump。新模块 `src/assembly/hic_pt.rs`。

**H3 短读投票映射器**（hic.cpp 936–2204）：`collect_votes`(958)、`get_longest_hit`(1058)、`compress_mapped_pos_advance`(1135)、`get_alignment`(1257)、`set_pe_pos_hap/set_pe_pos`(1695/1789)、`get_5_3_list`、`is_unreliable_hits`、`worker_for_alignment`(1918)/`worker_pipeline`(1976，单线程顺序等价)、`dedup_hits`/`sort_hits`、`write/load_hc_hits`(5236/5399)。

**H4 距离模型**（hic.cpp 3710–5170 活跃部）：`push_hc_edge`(3256)、`get_enzyme_occ`(3334)、`collect_hc_links`(4576)、`LCA_distance`/`get_LCA_bubble`/`get_LCA`(3710/3730/3806)、`worker_for_dis`(3884)、`fill_utg_distance_multi`(3968)、`get_trans_ug_arch`/`update_ug_by_trans`(4007/4024)、`push_LCA_edges(_rev)`(4104/4162)、`up_contain`(4237)、`update_containment_distance`(4487)、`measure_distance`(4952)、`append_tig_link`(4711)、`get_hc_edge`(6744)、`get_trans_rate_function_advance`(9986)。跳过死代码（min-cut 5434–6711、H_partition 10216–12805，~3.9k 行）。

**H5 tangle/权重/驱动**（hic.cpp 13510–17162 活跃部）：`clean_sub_tangle`(13705)、`resolve_bubble_chain_by_hic`(13990)、`resolve_tangles_hic`(16259)、`dfs_weight_hic`/`get_chain_weight_hic`(13315/13350)、`idx_hc_links`(15321)、`trans_checking_pass`(15360)、`weight_kv_u_trans`(15434)、`interpr_hit`(15516)、`get_interval_weight`/`get_hits_weight`/`get_trans_interval_weight`、`append_trans_hits`(15797)、`merge_u_trans_list`(15869)、`adjust_weight_kv_u_trans_advance`(15893)、`renew_kv_u_trans`(16003，含 `--unskew` 位)、`verbose_het_stat`(16146)、`renew_idx_para`(16754)、`hic_short_align`(17016)。

**H6 输出接线（全量）**：`hic_analysis`(17706) → `output_hic_graph`(Overlaps.cpp:16244，替换 `pipeline.rs:506` Err 桩；复用已移植的 `mc_solve`/`identify_bubbles`/`adjust_utg_by_primary`/`output_trio_graph_joint`)；多倍体 `output_hic_graph_mmhap`(16504) + mmhap 簇（hic.cpp 17433–17704，`bp_solve`/`purge_phase`）；trio+Hi-C `benchmark_hic_graph`(17793，替换 `pipeline.rs:435` 桩)；Hi-C+`--dual-scaf` scaffolding（horder.cpp `update_switch_unitig` misjoin 簇 ~1.2k 行 + `ha_aware_order`/`layout_scg` ~2.2k 行）。门禁：`.hic.tlb.bin`/`.hic.lk.bin`（未初始化字段按既有豁免）+ 最终 GFA/BED vs `reference/hic.*` + 日志锚点。

**H7 解锁**：`check_supported` 移除 `--h1`/`--h2` 拒绝（`parser.rs:953`）+ 全量回归 + 文档收尾。

## 关键文件

- Rust 新增：`src/assembly/hic_pt.rs`、`src/assembly/hic_hits.rs`（或合并为 `hic.rs`，随 H2–H5 逐步成文）、`scripts/gen_hic_reads.py`
- Rust 修改：`src/assembly/pipeline.rs`（373 跳过点、435/506 Err 桩、H6 接线）、`src/assembly/ul_renew_driver.rs`（hic_clean_adv/deep_graph_clean no-op → 实装）、`src/command_line/parser.rs`（953 解锁）
- 参考对照（只读）：`/home/keith/projects/hifiasm/{hic.cpp,horder.cpp,Overlaps.cpp,gfa_ut.cpp}`
- 文档：`docs/consistency.md`（每阶段更新，终态「仍未迁移的分支」仅剩 `--dbg-gfa`）

## 验证

- 每簇：`cargo test --release` 全绿（含新增 `--ignored` 黄金回放）。
- 每阶段：`scripts/compare_outputs.sh` + `cmp` 逐字节；UL 用 `/tmp/ulref`；Hi-C 用 H0 新参考 `reference/hic.*`。
- 最终门禁矩阵：`--ul`、Hi-C 端到端逐字节；默认 9-profile、trio、`-5`、`--trio-dual`、`--dual-scaf`、yak、ONT、`--fast`、`--sec-in`（已知 UB 豁免）无回归。
- 终态：consistency.md 残留仅 `--dbg-gfa`（维持显式拒绝）。