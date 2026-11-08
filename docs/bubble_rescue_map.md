# 气泡 rescue 系移植地图（trio Stage A 门 A 的唯一阻塞 → Hi-C Phase 3 共用）

目标：移植 `rescue_bubble_by_chain`（Overlaps.cpp:36832）及其全部传递依赖，
使 trio 门 A（`-3/-4` hap1/hap2 逐字节）通过。全部机械随后被 Hi-C Phase 3
复用（hic.cpp 气泡半区）。当前已验证：清洗 rounds/semi_circ 之前与 C++
逐字节一致（见 consistency.md「Trio Stage A 进行中」节），分歧只在
post_rescue 的这个函数（C++ 2387→2396 活弧，Rust 停在 2390）。

## 已确认的 Rust 可复用件

- `Asg`（push_arc/cleanup/delete_arc/symmetrize/delete_sequence/sequence_visited）
- `ReadToUnitig::{set,set_raw,get}`（hit.rs:1130-1162，= set_R_to_U/get_R_to_U）
- `ma_ug_gen`（unitig.rs:80，返回 UnitigGraph；C++ ma_ug_gen_primary 对应
  `ma_ug_gen_primary`——unitig.rs 也有）
- `copy_read_graph`/`copy_untig_graph`（unitig.rs）
- `adjust_utg_by_primary`（bp_graph.rs:169；注意 rescue_bubble_by_chain 的
  C++ 调用第三参传 TRIO_THRES=0.9——对照 C++ 原型 Overlaps.h:1067 确认
  Rust 形参位置）
- `asg_bub_pop1_primary_trio`（bubble_identify.rs:440；返回 1/0=是否成泡；
  **缺 path_base_len 出参**——由独立的
  `asg_bub_backtrack_primary_length`（Overlaps.cpp:24135-24182，需新移植）
  在成功后补算，勿改既有函数）
- `asg_bub_pop1_primary_trio_switch_check`（hamming.rs:957）
- `BubbleBuffer`（bubble.rs:32：sources=b.S、visited=b.b、visited_edges=b.e）、
  `BubbleInfo.parent`（= b.a[v].p）
- `get_s_bub_pop_max_dist_advance`（bubble_identify.rs:269）
- `pre_clean`（rescue.rs:226）、`lable_all_bubbles`（rescue.rs:808）、
  `get_edge_from_source`（unitig.rs:603）
- `HapCov`/`TransChain`（hap_cov.rs:187/105，t_ch 含 k_trans/k_trans_idx/
  ir_het/bed——由 adjust_utg_by_primary 填充）

## 移植清单（自底向上，每簇编译一次）

### 簇 1：bubble_type 与基础件（✅ 已移植，src/assembly/bubble_chain.rs，
### 2026-11-08：BubbleType/Pdq(reset/push/pop/shortest_path)/get_bubbles/
### dfs_bubble(_broken)/get_unitig_het_arb/asg_bub_backtrack_primary_length/
### update_bub_b_s_idx/get_bub_id/connect_bub_occ/IF_HET/IF_HOM 全部就位，
### cargo test 294 全绿）
### 簇 1 原始清单：
| Rust 项 | C++ 位置 | 备注 |
|---|---|---|
| `BubbleType` | hic.h:37-53 | index:Vec<u32>(n_vtx)、round_id、n_round、ug、list/num/pathLen/b_s_idx:Vec<u64>、s_bub/f_bub/b_bub/b_end_bub/tangle_bub/cross_bub/mess_bub:u64、check_het:u8、b_g/b_ug、chain_weight（trio 路径恒空）、c_w（Hi-C 专用，先省略） |
| `M_het/P_het/IF_HET` | hic.h:82/83/88 | M=num.n、P=num.n+1、IF_HET(ID)=index[ID]==f_bub+1 |
| `Pdq` | hic.h:66-78 | vis/x/dis 三数组 + uID_mode/uID_shift/tmp_v/tmp_d；init/destory（hic.h:109-110，实现 hic.cpp 搜 init_pdq） |
| `get_bubbles` | hic.cpp:2217 | 返回 (beg,sink,&list[num[id]+2..],pathLen[id]) |
| `dfs_bubble` | hic.cpp:2286 | 双向栈 DFS，结果=途经节点（cur>>1） |
| `dfs_bubble_broken` | hic.cpp:2232 | 同上带 beg/sink 屏蔽 |
| `get_unitig_het_arb` | hic.cpp:2327 | P_HET=2/C_HET=1（Overlaps.h:962-963）；ref→u_trans 判空即可（trio 传 ref=NULL） |
| `asg_bub_backtrack_primary_length` | Overlaps.cpp:24135 | 经 BubbleInfo.parent 回溯；首跳加 seq.len、后续加 arc.ul 低 32 |
| `update_bub_b_s_idx` | hic.cpp:7663 | b_s_idx 高低 32 位各存一个 bubble id（注意 C 的 `<<=32` 移位语义） |

### 簇 2：identify_bubbles（✅ 已移植 bubble_chain.rs，含 P_HET/C_HET 掩码、
### 三遍扫描、dfs_bubble 单链校验、radix_sort_by_key_32、het 修正与
### b_s_idx 登记；pathLen 经 asg_bub_backtrack_primary_length 补算；
### round_id!=0 分支同步移植。尾部调用的 build_bub_graph 目前是
### unimplemented! 桩（无调用者，不可达））
### 簇 2 原始清单：
- 两遍 asg_bub_pop1_primary_trio 扫描（第一遍标 index 低 2 位：体内+1、
  beg+2、sink^1+3；第二遍收集 num/list/pathLen，b.b.n∈{2,3,5} 时用
  dfs_bubble 校验单链，失败置 num 元素 |1<<31）；
- s_bub 计数、f_bub=num.n-1；
- 之后 get_unitig_het_arb 填 index（m_het_occ=20）、按 f_bub 循环做
  beg/sink 的 index 修正与 b_s_idx 登记（2455-2500）、M_het→P_het 归并；
- 入参 r_het_flag = cov.t_ch.ir_het（trio 传 ref=NULL）；
- 尾部 `build_bub_graph`（→ 簇 3）。

### 簇 3：bubble graph 与链（✅ 已移植 2026-11-08，bubble_chain.rs ~2670 行；
### ma_2_bub_arc/set_bub_idx/determine_bub_idx/detect_bub_graph/get_bub_graph/
### is_simple_broken_bubble/double_check_broken_bubble/is_local_simple_circle/
### update_bubble_graph/get_related_bub_nodes/calculate_chain_weight/get_btg_occ/
### set_b_utg_weight_flag/check_bubble_tangle/find_bubble_tangle/get_utg_end_from_btg/
### drop_g_edges_by_utg/update_bsg/resolve_bubble_chain_tangle/update_bubble_chain/
### build_bub_graph（真实现，桩已除）。要点：BubbleType 增加 bubble_chain_graph
### 字段（b_ug->g）；链数据以参数传递规避借用；chain_weight 排序对 trio 输出
### 不可见（仅 Hi-C c_w 消费），用稳定 sort_by。cargo test 294 全绿）
### 簇 3 原始清单：
| Rust 项 | C++ 位置 | 行数 |
|---|---|---|
| `detect_bub_graph` | hic.cpp:7271 | 70（依赖 determine_bub_idx hic.cpp:7254、IF_HET、arc_first/arc_cnt——asg.h 内联，Rust 对应 arc_range） |
| `get_bub_graph` | hic.cpp:7341 | 79（依赖 connect_bub_occ:6962、ma_2_bub_arc:6972、get_shortest_path:3495、Pdq、asg_seq_set） |
| `is_simple_broken_bubble` | hic.cpp:7476 | 49 |
| `double_check_broken_bubble` | hic.cpp:7525 | 116 |
| `is_local_simple_circle` | hic.cpp:7641 | 22 |
| `update_bubble_graph` | hic.cpp:7705 | 102 |
| `get_related_bub_nodes` | hic.cpp:7807 | 85 |
| `check_bubble_tangle` | hic.cpp:8092 | 93 |
| `find_bubble_tangle` | hic.cpp:8185 | 196 |
| `drop_g_edges_by_utg` | hic.cpp:8381 | 63 |
| `update_bsg` | hic.cpp:8513 | 49 |
| `resolve_bubble_chain_tangle` | hic.cpp:8562 | 174 |
| `update_bubble_chain` | hic.cpp:8736 | 224（依赖 calculate_chain_weight/double_check_broken_bubble/get_related_bub_nodes/set_b_utg_weight_flag/update_bubble_graph/update_bub_b_s_idx/pdq；**不含 hc_links**） |
| `set_b_utg_weight_flag` | hic.cpp:8960 | 33 |

明确不移植（Hi-C 专用，trio 不可达）：print_bubble_chain(7420, 调试)、
resolve_bubble_chain_tangle_back(7913, hc_links)、debug_tangle_bubble(8444)、
set_b_utg_weight_flag_xor(8993)、dfs_weight(9025)、if_conflict_utg(9060)、
get_chain_weight(9105)、double_check_bub_branch(9270)、
clean_bubble_chain_by_HiC(9324)、append_boundary_chain(9516)、
init_chain_hic_warp(9687)。
（若 resolve_bubble_chain_tangle/update_bubble_chain 体内实际调用上述任一，
以 C++ 调用图为准补入——移植时逐一核对。）

### 簇 4：三个 rescue + 驱动（Overlaps.cpp）
| Rust 项 | C++ 位置 | 行数 |
|---|---|---|
| `get_edge2existing_node_advance_by_broken_bub` | Overlaps.cpp:33409(inline) | ~90 |
| `get_contained_reads_chain_by_broken_bub` | Overlaps.cpp:34276 | 102 |
| `set_rtg_flag_by_bubble` | Overlaps.cpp:35321 | 55 |
| `minor_transitive_reduction_r_g` | Overlaps.cpp:35509 | 32 |
| `if_recoverable` | Overlaps.cpp:35541 | 79 |
| `bub_complex` | Overlaps.cpp:36549 | 171（用 switch_check + R_INF.trio_flag↔trio_labels） |
| `rescue_bubbles_by_contained_reads` | Overlaps.cpp:35620 | 263 |
| `rescue_bubbles_by_missing_ovlp` | Overlaps.cpp:36092 | 154 |
| `rescue_bubbles_by_missing_ovlp_backward` | Overlaps.cpp:36247 | 287（用 update_unitig 36891/36246） |
| `reset_bub` | Overlaps.cpp:36534 | 15（destory+memset+identify_bubbles） |
| `rescue_bubble_by_chain` | Overlaps.cpp:36832 | 59（ma_ug_gen_primary → copy_read_graph/copy_untig_graph → adjust_utg_by_primary(第三参 TRIO_THRES) 出 cov → reset_bub ×3 + 三 rescue；trio_recover 块已注释；cmk=NULL 跳过第 4 个） |
| `update_unitig` | Overlaps.cpp:36891 | 60 |

rescue 系的 R_INF.trio_flag 全部改传 `trio_labels: &mut [u8]`（bub_complex
会写 DROP→AMBIGU 回置，需要 &mut）；`r_g->seq[].c`（coverage_class）与
`coverage_cut[].del/.c` 的读写对齐 Rust Asg/CoverageCut 字段名。

### 簇 4（✅ 已移植 2026-11-08，~1900 行 Rust：get_num_edges2existing/
### get_edge2existing 的 broken_bub 变体、get_contained_reads_chain_by_broken_bub、
### set_rtg_flag_by_bubble、minor_transitive_reduction_r_g、三 rescue、reset_bub、
### rescue_bubble_by_chain 驱动；已接入 pipeline post_rescue d15 后——暂限
### trio_binning_requested() 才进入，非 trio 语料 C++ 证明 no-op）。
###
### ⚠️ 新发现的下一个缺口（门 A 的现行阻塞）：
### rescue 内部的 adjust_utg_by_primary 在 trio 标签下触发 Purge_Dups 的
### hap_alignment_advance_worker 候选扫描（Purge_Dups.cpp:2964-3179，215 行，
### 默认语料不触发）——含 get_xy_pos/deduplicate_edge/get_candidate_hap_alignment
### (DP 链化)/calculate_pair_hap_similarity_advance/prefilter/
### filter_secondary_chain/cmp_hap_candidates 助手簇，估 ~2000 行 C++。
### 这是此前 purge_dups 深链 stub 清单中的"hap_alignment_advance_worker
### candidate scan"项。另：簇 3 有一个未解 panic（默认语料曾触发
### shortest_path len=92 vs src=217345——trio 门控后默认不进入；trio 语料
### 走到该处前会先撞上述 Err，故暂不可达；簇 4 完成后 trio 全链路会暴露它，
### 需用 RUST_BACKTRACE=1 重新定位）。
###
### 簇 4 原始清单（ historical）：Overlaps.cpp 33409-36891，~1400 行 C++：
### get_edge2existing_node_advance_by_broken_bub(33409,inline ~90)/
### get_contained_reads_chain_by_broken_bub(34276,102)/set_rtg_flag_by_bubble(35321,55)/
### minor_transitive_reduction_r_g(35509,32)/if_recoverable(35541,79)/bub_complex(36549,171)/
### rescue_bubbles_by_contained_reads(35620,263 已通读)/rescue_bubbles_by_cmk_reads(35883,cmk 门控跳过)/
### rescue_bubbles_by_missing_ovlp(36092,154)/rescue_bubbles_by_missing_ovlp_backward(36247,287)/
### update_unitig(36891,60)/reset_bub(36534,15)/rescue_bubble_by_chain(36832,59 已通读)。
### R_INF.trio_flag → trio_labels: &mut [u8]。

### 接入点（簇 4 完成后执行）
pipeline.rs post_rescue 等价段的 d15_after_set_hom_cov dump 之后
（pipeline.rs:833 的注释处）：
`rescue_bubble_by_chain(graph, coverage_cut, sources, reverse_sources,
read_to_unitig, options, trio_labels, cov?, bubble_labels)`——对照 C++
post_rescue(gfa_ut.cpp:3013) 末行实参 `(asm_opt.max_short_tip*2, 0.15, 3,
rI, 0.05, 0.9, max_hang, min_ovlp, 10, gap_fuzz, b_mask_t, no_trio_recover=0, cmk=NULL)`。
注意 cov 来源：C++ 在函数内部经 adjust_utg_by_primary 生成（出参 &cov），
Rust 同样在函数内生成，不从外部带入。

## 验证

1. 每簇编译 + 现有单测不回归。
2. 门 A 重跑：`-3 l1 -4 l2 -t1` vs C++ reference_trio，GFA/bed 全部逐字节；
   中途可用 gdb 断点（gfa_ut.cpp:3188 = post_rescue 后）对比 p14 状态。
3. 9-profile 默认回归（该函数在默认语料是 no-op，必须保持 d15==d16）。
4. 完成后：拆掉 pipeline.rs 的临时探针 d00a/d00b（或转正），Stage A commit。
