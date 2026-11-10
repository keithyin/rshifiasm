# hifiasm → Rust 迁移：输出层收尾（break_ug_contig + 端粒 + dbg-ovec）

## Context

默认 HiFi 链路 P0–P6 与选项空间基础包（--primary/-e/重启缓存）已全部逐字节一致。本轮完成三个"现有数据即可门禁"的小项，使输出层不再有任何显式 Err 路径、默认 HiFi 选项空间收尾：
1. `break_ug_contig`（`--b-cov/--h-cov/--m-rate`，Overlaps.cpp 11503–12565，~560 行）——移除 bp_graph.rs 三处 "not ported" 守卫
2. 端粒 `--telo-m/-p/-d/-s`（~155 行读段扫描）——消费侧 `TeloEnd` 已全部移植于 gfa_clean.rs，只缺生产侧与参数穿线
3. `--dbg-ovec`（~150 行）——纠错链化/比对层的字节级验证检查点

探索确认：三者均单线程确定性；break_ug_contig 的全部下游 helper 已移植；三处守卫在 bp_graph.rs:417（print_utg）、917（output_trio_unitig_graph）、1734（output_contig_graph_primary）。

## 阶段一：break_ug_contig

新模块 `src/assembly/break_contig.rs`，按 C++ 逐函数移植：
- `InSubInterval`（kv_in_sub_t，Overlaps.cpp 11593）、`push_cov_interval_advance`(:11503)、`get_break_point_cov_advance`(:11537，经 sources+get_R_to_U 的同 unitig 覆盖深度扫描)
- `detect_break_point_advance`(:11866，两遍扫描 + radix sort + 区间合并 + 谷/峰跃迁 + 圆形回绕)
- `get_overlap_contig_dir`(:11620) + `get_break_point_idx`(:11659，m_rate 精确匹配比例判定)
- `push_sub_unitig`(:12313，复用已移植 fill_unitig unitig.rs:2849) + `break_all_contigs`(:12340，拆分/弧重映射/renew_utg) + `break_ug_contig`(:12521 入口)
- 移除 bp_graph.rs 三处 Err，替换为 `break_ug_contig(...)` 调用（参数：options.breakpoint_low_coverage/breakpoint_high_coverage/exact_match_rate）
- 注意：`push_cov_interval_advance` 的返回值门控第二遍的贡献（圆形回绕用例），`c_break` 拒绝非圆形端点切割——这些细节决定字节一致性
- 门禁：C++ `hifiasm --b-cov 20 -t1`（及 `--h-cov`、`--primary --b-cov 20` 变体视需要）参考运行 → Rust 同命令 → compare_outputs.sh 全部 gfa 逐字节一致

## 阶段二：端粒

- `gen_telo_ends`（mark_telo_ends + gen_telo_end_t，Overlaps.cpp 191–347，~155 行）——建议放 gfa_clean.rs（`TeloEnd` 定义处），单线程化（-t1 门禁）
- 穿线：`options.telomere_motif/penalty/drop/min_score` → `ul_clean_gfa`（pipeline.rs，新增参数）→ `asg_arc_cut_tips`/`asg_iterative_semi_circ`/`asg_arc_cut_chimeric`（已接受 `Option<&TeloEnd>`，全部现成）
- 门禁：C++ `hifiasm --telo-m CCCTAA -t1` vs Rust。已知风险：chr11-2M 是内部区域，可能无端粒读段（C++ 会打 5'/3' 计数摘要，若为 0 则该语料上为 no-op——如实记录，默认奇偶性已由现有门禁保证）

## 阶段三：--dbg-ovec

- 实现 `error_correction_debug`（src/assembly/mod.rs 桩位）：无 filter table 建索引 → 按 2000 读分块 → 每读复用已移植的 `lchain_qgen_mcopy_fast`（ec_driver.rs:354）与 `gen_hc_r_alin_ea`（hc_r_alin.rs:291，chem_drop=0，gap -1/-1）→ `print_ov_dbg_paf` 输出带 `cg:Z` 的 PAF（对照 ecovlp.cpp 6311–6452 移植 worker 与打印）
- 门禁：C++ `hifiasm --dbg-ovec -t1` 的 `<out>.ovlp.paf` vs Rust 逐字节 cmp
- 价值：把 C++/Rust 分歧隔离在共识之前的链化/比对状态

## 阶段四：收尾

1. `cargo test --release` 全量 + golden_pipeline 2 条 e2e（串行）
2. `docs/consistency.md`：新增本轮三节门禁实测 + 更新"选项空间剩余缺口"（移除已完成项，剩 --fast/--dual-scaf/trio/Hi-C/UL/pair）
3. 提交（按阶段 1–3 个提交）

## 关键文件

- 新 `src/assembly/break_contig.rs`（~700 行）
- `src/assembly/bp_graph.rs`（3 处守卫替换）
- `src/assembly/gfa_clean.rs`（端粒扫描）+ `src/assembly/pipeline.rs`（端粒穿线）
- `src/assembly/mod.rs` + 纠错驱动复用（dbg-ovec）
- `docs/consistency.md`

## 风险与注意

- break_ug_contig 的 b_low/b_high 共用一个区间缓冲但调用点从不同时启用；radix sort 的 64 位键含符号位约定（gl 案例教训）
- 端粒门禁可能在本语料上为 no-op（如实记录）
- dbg-ovec 的 .ovlp.paf 文件名与 --write-paf 冲突（模式不共存，无碍）
- C++ 参考需 3 次运行（--b-cov 20、--telo-m CCCTAA、--dbg-ovec），各 ~2 分钟

## 验证

1. 三个 profile 的 `scripts/compare_outputs.sh`（gfa 全绿；dbg 用 .ovlp.paf cmp）
2. `cargo test --release` + `cargo test --release --test golden_pipeline -- --ignored`
3. 每阶段推进同步 docs/consistency.md