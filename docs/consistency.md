# 结果一致性记录

本文件记录与 C++ 版 hifiasm 的逐字节比对基线、各阶段门禁状态与已知差异。
每次阶段推进都需要更新本文件。

## 基准环境

| 项目 | 值 |
| --- | --- |
| C++ 源码 | `/home/keith/projects/hifiasm`（只读参考） |
| C++ 版本 | `0.25.0-r726` |
| 测试数据 | `/home/keith/Downloads/chr11-2M.fa.gz`（5626 条 HiFi reads，84042542 bp） |
| 线程数 | 1（单线程，保证确定性） |
| 基准命令 | `hifiasm -t1 -o reference/reference chr11-2M.fa.gz` |

## 基准产物校验和

由 `scripts/reference_run.sh 1` 生成，结果写入 `reference/`（已被 `.gitignore` 忽略）。

| 文件 | SHA-256 |
| --- | --- |
| `reference.bp.p_ctg.gfa` | `2c0441e881aebaab32373211f0b1030fe6c0d0b8280e7c72663f15c81001ac0e` |
| `reference.bp.p_utg.gfa` | `dfba9927e9e0a8161622bbbd345730bd6af35c919db0eea75b80a2e1ba4d0016` |
| `reference.bp.hap1.p_ctg.gfa` | `9fa6a9263a5afa48191e9eebf71b6c854c5d8b2ceb24786633fbea4111c8bb09` |
| `write_paf_ec.ec.fa` | `c719794bd60738766fc46b3ca2b2ee72a2800a9d7a01e38b422b826af28383db` |
| `write_paf_ec.ovlp.paf` | `dc2b4d1d4452dd54f6071ea323bf702f76084d9aabc954f311665491e4508eb2` |

## 已确认的事实

1. **单线程输出可复现**：同一命令运行两次，`bp.p_ctg.gfa`、`bp.p_utg.gfa`、
   `bp.hap1.p_ctg.gfa`、`bp.hap2.p_ctg.gfa` 完全一致（`cmp` 无差异）。
   因此逐字节比对标准是可达到的，不需要放宽。
2. **运行时间**：单线程约 115–120 秒，峰值 RSS 约 16.5 GB。
3. **读段统计**：Rust 版读入后得到 `5626 reads, 84042542 bp`，与 C++ 日志中的
   `# bases: 84042542` 完全一致，说明读段过滤（adapter 裁剪、长度阈值）语义等价。

## 阶段门禁状态

| 阶段 | 门禁 | 状态 |
| --- | --- | --- |
| P0 基准 | 生成 `-t1` 基准产物与比对脚本 | ✅ 完成 |
| P1 基础设施 | `-h` / `--version` 文本与 C++ 逐字节一致 | ✅ 完成 |
| P1 基础设施 | 读段数量与总碱基数与 C++ 一致 | ✅ 完成 |
| P1 基础设施 | khashl 桶布局与 C++ 一致（黄金布局对拍） | ✅ 完成 |
| P2 k-mer 计数 | `peak_hom` / `peak_het` / 过滤 k-mer 数与 C++ 一致 | ✅ 完成（60 行直方图日志逐字节一致） |
| P3 minimizer | `collected N minimizers` / 位置索引前的计数与直方图 | ✅ 完成（第二轮 60 行直方图逐字节一致） |
| P3 位置索引 | `indexed 2224445 positions` | ✅ 完成 |
| P3 编辑距离 | 64 位位并行 `ed_band_cal_*_w` 与 C 逐例一致 | ✅ 完成（非回溯 4 支 + global/semi trace，21 组黄金用例含 CIGAR） |
| P3 链化 | `chain_DP` 与真实 `Hash_Table.o` 逐字段一致 | ✅ 完成（4 组黄金用例） |
| P3 候选种子 | `ha_get_new_candidates` 的种子收集与权重映射 | ✅ 完成（单测覆盖镜像位置与计数权重） |
| P3 overlap 追加 | `append_inexact_overlap_region_alloc` 与真实 `.o` 逐字段一致 | ✅ 完成（4 组黄金用例，含 f_cigar 全部条目） |
| P3 每条读的驱动 | `ha_get_new_candidates` 的链数限制与 `ha_ov_type` 分类 | ✅ 完成（单测覆盖 4 类 containment） |
| P3 overlap 落盘 | `push_overlaps` 与 `Output_PAF` 的字段/格式 | ✅ 完成（含反向坐标镜像的单测） |
| P3 overlap 精化 | `append_window_list` / `ovlp_chain_gen` / `worker_ovec` | 旧路径遗留，非缺口（2026-11-07 核实：仅 `ha_overlap_cal`（`--sec-in` pair 模式）与 C++ 死代码 `rescue_hp_reads`/`ha_overlap_and_correct` 走此路；默认纠错 `cal_ec_r` 与 `--dbg-ovec` 的 `cal_ec_r_dbg` 均不经此，Rust 已由新路径逐字节覆盖。随 `--sec-in` 前端一并实现） |
| P4 纠错 cigar | `Cigar_record` 打包、合并、`lost_base` 累积 | ✅ 完成（与真实 `.o` 三组黄金用例一致） |
| P4 POA 图 | `Node`/`Edge` 布局与四类边权重更新 | ✅ 完成（整图逐节点逐边对拍 C 探针） |
| P4 窗口 cigar | `window_list` 的 `op<<14\|len` 打包与共享 cell 缓冲 | ✅ 完成 |
| P4 滑窗调度 | `Window_Pool` / `get_interval` 的分组与游标 | ✅ 完成（对拍真实 `Correct.cpp`） |
| P4 目标区域映射 | `determine_overlap_region` 的 padding 与拒绝条件 | ✅ 完成（9 组黄金用例） |
| P4 带状比对 | `Reserve_Banded_BPM` 的 Myers 位并行与 band 滑动 | ✅ 完成（50 组黄金用例，含 N） |
| P4 窗口验证 | `y_start_offset` + `verify_window` 组合 | ✅ 完成（对拍真实 `Hash_Table.h`，4 路 SSE 已证明等价于 4 次标量） |
| P4 阈值与重算 | `get_init_err_thres` / `double_error_threshold` / `append_unmatched_wins` / `simi_pass` | ✅ 完成（阈值 7 组、窗口散列 4 组黄金用例） |
| P4 比对回溯 | `Reserve_Banded_BPM_PATH` 的 DP 矩阵与路径回溯 | ✅ 完成（15 组黄金用例，含 indel 与失败用例） |
| P4 cigar 生成 | `move_gap_greedy` / `generate_cigar` 的裁剪、贪心挪 gap、游程编码 | ✅ 完成（5 组黄金用例，逐 cell 比对） |
| P4 边界修正 | `fix_boundary` 的左右移窗与"仅接受更优"判定 | ✅ 完成（6 组黄金用例） |
| P4 同聚物/可用区间 | `if_is_homopolymer_strict` / `get_available_interval` | ✅ 完成（11+3 组黄金用例） |
| P4 POA 图操作 | 双向边定位、visited 复用 `length`、节点删除与三类边遍历 | ✅ 完成 |
| P4 insertion DAG | `build_DAGCon` / `Merge_DAGCon` / 最优路径 | ✅ 完成（4 组黄金用例） |
| P4 一致性回溯 | `get_seq_from_Graph` / `add_path_to_correct_read_new` | ✅ 完成（4 组黄金用例） |
| P4 窗口共识 | `window_consensus` + `generate_consensus` 单轮主循环 | ✅ 完成（6 组长度黄金用例） |
| P4 二轮纠错 | `process_boundary` / `add_cigar_to_cigar` / `merge_cigars` | ✅ 完成（merge_cigars 3 组黄金用例） |
| P4 读段写回 | `worker_sl_ec` 等价的 cigar 展开、`Output_corrected_reads` / `Output_PAF` | ✅ 完成（cigar 展开 3 组黄金用例） |
| P4 新链路起步 | `set_lchain_dp_op` / `minimizers_qgen0` 种子权重与排序 | ✅ 完成（参数 10 组黄金用例） |
| P4 延伸比对 | `Reserve_Banded_BPM_Extension` / `alignment_extension` 的双向延伸 | ✅ 完成（16 组黄金用例，双向） |
| P4 空洞重算 | `verify_sub_window` / `get_init_paras` / `non_trim_error_rate` | ✅ 完成（5 组误差率黄金用例 + 双向窗口坐标） |
| P4 覆盖判定 | `get_available_fully_covered_interval` / `check_if_fully_covered` | ✅ 完成（7 组黄金用例，含空 overlap 的 `-2` 提前返回） |
| P4 边界 cigar | `scan_cigar` / `scan_cigar_interval` / `calculate_boundary_cigars` | ✅ 完成（3 组黄金用例，含逐 cell 比对与 `-2` 跳过） |
| P4 窗口重算主循环 | `recalcate_window_advance` 的补洞、补 cigar、`is_match` 分类 | ✅ 完成（4 组结构用例；`align_length` 累加语义与 C 一致） |
| P4 单读驱动 | `correct_overlap`：验证 → 重算 → 共识 → 覆盖判定 | ✅ 完成（4 组用例，含 6 条 overlap 的满覆盖） |
| P4 纠错流水线 | `ha_ec`：逐读建链 → 纠错 → 写回 `ReadStore` | ✅ 完成（`--write-ec` 全量跑通，162.9s / 7.5GB） |
| P4 `--write-ec` 门禁 | 与 `reference/write_paf_ec.ec.fa` 逐字节一致 | ✅ 完成（84,247,384 B，`cmp` 无差异） |
| P4 新链化打分 | `cal_bw` / `comput_sc_ch_ec` | ✅ 完成（8 组黄金用例，含带宽拒绝与惩罚分支） |
| P4 新链化包装 | `lchain_qgen_mcopy_fast` 的分组/限链/覆盖率窗/r485 | ✅ 完成（W1–W5 五组黄金用例，逐字段一致） |
| P4 新链路比对 | `align_hc_ed_post_extz` / `push_hc_wlst_exz` / `aln_wlst_adv_exz` / `recal_boundary_exz` | ✅ 完成（12 组黄金用例，含反向链、左右边界与 backtrace CIGAR） |
| P4 延伸误差 | `gen_extend_err_exz` / `gen_extend_err_0_exz` / `ed_cut`（含反向端点约定） | ✅ 完成（11 组黄金用例；含双端补洞、补洞 backtrace 与 `e_max` 提前中止） |
| P4 exz 比对层 | `init_waln` / `adjust_ext_offset` / `update_semi_coord` / `cal_estimate_err(_hc)` / `cal_exact_exz` / `cal_exz_infi_adv` | ✅ 完成（4 组场景、381 行黄金回放逐字节一致） |
| P4 比对驱动 | `hc_aln_exz_adv` / `hc_aln_exz_adv_hc` 的阈值阶梯与早期返回 | ✅ 完成（4 组场景 × 16 组参数、182 行黄金回放逐字节一致） |
| P4 整段比对 | `hc_ovlp_base_direct` 的精确窗口捷径、`pre_mode` 分段与未比对窗口 | ✅ 完成（10 组场景 × 25 次调用、70 行黄金回放逐字节一致） |
| P4 窗口重链 | `rechain_aln_hc` 的重试链化、整段比对与旧窗替换 | ✅ 完成（3 组读段 × 4 组重试、691 行黄金回放逐字节一致） |
| P4 快速比对驱动 | `gen_hc_fast_cigar0` 的链分组、首轮整段比对、重试与窗口整理 | ✅ 完成（6 组读段 × 4 组参数、216 行黄金回放逐字节一致） |
| P4 空洞重排 | `reassign_gaps` / `move_wins` / `adjust_gap` / `ajust_end_cigar` | ✅ 完成（9 组对齐场景 + 6 组手写 cigar、62 行黄金回放逐字节一致） |
| P4 整段精化 | `gen_hc_r_alin` / `ff_lunalign` / `ff_tend` 的逐 overlap 流程与就地压缩 | ✅ 完成（24 组场景 + 6 组终端过滤、111 行黄金回放逐字节一致） |
| P4 精确重叠回收 | `gen_hc_r_alin_ea` / `gen_hc_r_alin_nec` 的 `el` 捷径、首条命中查询与已匹配跳过 | ✅ 完成（18 组场景、106 行黄金回放逐字节一致） |
| P4 已存 overlap 辅助 | `dedup_chains` / `push_nec_re` / `extract_max_exact` / `push_ne_ovlp` / `check_well_cal` | ✅ 完成（22 组场景、149 行黄金回放逐字节一致） |
| P4 纠错 | `--write-ec` 纠错读段逐字节一致 | ✅ 完成（三轮全部一致，见下文实测） |
| P5 组装图 | `p_utg.gfa`、`p_ctg.gfa` 逐字节一致 | ✅ 完成（r_utg/p_utg 全链路一致；p_ctg 探针级 + e2e 一致，见下文实测） |
| P6 单倍型 | `hap1/hap2.p_ctg.gfa`、`a_ctg.gfa` 一致 | ✅ 完成（`hap1/hap2.p_ctg.*` 探针级 + e2e 逐字节一致，见下文实测；`--dual-scaf` 的完整路径见下文 `trans_base_infer` 一节） |
| P7 Hi-C/trio | trio 全部分支逐字节；Hi-C 需额外测试数据 | ✅ trio（名单模式/`--trio-dual`/`-5`/`--dual-scaf`）；⏳ Hi-C |
| P8 ONT/UL | 需额外测试数据 | ⏳ 待实现 |

## 已验证的翻译细节

- `ketopt` 语义：短选项聚合、长选项前缀缩写、`--key=value`、参数置换，
  以及 `[ERROR] unknown option in "..."` 等报错文本。
- `printf("%.2g")` / `("%.3g")` 语义：`0.8` 而不是 `0.80`，已由
  `utilities::printf_format::format_general` 覆盖并有单测。
- 4-bit 碱基压缩：非 A/C/G/T 存为 `A` 并记录位置，正反链子区间还原一致。
- khashl：Fibonacci 哈希 `hash * 2654435769 >> (32 - bits)`、线性探测、
  0.75 装载因子扩容、kick-out 重散列；桶布局与 C++ 逐桶一致。
- 计数表的相等性判定忽略低 12 位（计数打包在键里），修改键即计数自增。
- blocked bloom filter：每分区 `1<<(bf_shift-pre-3)` 字节、64 字节块、
  4 次探测，`hash>>shift` 与 `hash>>(shift-9)` 的分工与 C++ 相同。
- `yak_hash64_64` 的移位/加法序列与 C++ 完全一致（含 `~key + (key<<21)` 的
  求值顺序），这是直方图能逐字节重合的关键。

### P4 纠错 cigar（`correction::cigar`）

`CigarRecord` 逐条移植 `Cigar_record` 的四个 `add_*` 函数：低 2 位放操作类型、
高 30 位放长度；只有 match 直接累加，mismatch/insertion/deletion 还要把碱基
追加到 `lost_base`，且 deletion 不增加 `new_read_length`。黄金值与真实
`Correct.o` 对拍：

```
matches_merge            new_read_length=6 cells=1  [0] op=0 len=6
mismatches               new_read_length=6 cells=2  [0] op=0 len=4  [1] op=1 len=2   lost=GT
deletion_then_insertion  new_read_length=5 cells=3  [0] op=0 len=4  [1] op=3 len=2  [2] op=2 len=1  lost=GGC
```

`clear_Cigar_record` 把 `current_operation` 复位为哨兵值（C 用 `127`，Rust 用
`-1`，两者都不会与 0..3 的操作码相等，语义等价）。

### P4 POA 图（`correction::poa`）

三个边表沿用 C 的宏命名而非直观方向，这一点在多处容易读错，已在注释里写明：

- `deletion_edges` 存**出边**（`Output_Edges` 宏），
- `insertion_edges` 同时存**入边**（`Input_Edges` 宏）与该节点挂载的 insertion bubble，
- `mismatch_edges` 存该节点出发的替代碱基边。

删除边不移动位置，只把七个字段全部写成 `u64::MAX` 并递增 `delete_length`，
因此边下标（`self_edge_ID` / `reverse_edge_ID`）在删除后依然有效——这是
`clear`/`remove` 语义能被逐边复现的前提。

黄金值来自真实 `POA.h` 内联函数的探针，覆盖 `addUnmatchedSeqToGraph`、
重复折叠、mismatch、insertion bubble、deletion 五步，落地为整图渲染比对
（节点数、节点基/权重、每条边的 `in/out/weight/length/num_insertions`）。
两个易错点在测试中被锁定：新建 mismatch 节点后会补一条
`新节点 -> in_node + 1` 的返回边；单碱基 insertion 的返回边挂在**插入节点**
上而不是锚点节点上。

### P4 窗口与区域映射（`correction::window` / `correction::interval` / `correction::region`）

`window_list` 的 cell 是 `uint16`：高 2 位操作码、低 14 位长度
（`op << 14 | len`，长度上限 16383）。所有窗口共用一条 cell 缓冲，窗口只用
`cidx`/`clen` 指向自己那一段。

`get_interval` 把命中当前窗口的 overlap 分成两组，两个容易读错的地方都已对拍：

- **满窗组**从数组**头部**往后写（`overlapID[length++]`），遍历顺序即列表顺序；
- **部分窗组**从数组**尾部**往前写（`overlapID[size - lengthNT]`），而
  `verify_window` 用 `overlapID[size - i]` 取值，所以消费顺序同样是列表顺序。

`Window_Pool` 的最后一行窗口只发一次：当 `window_end == read_length - 1` 时
先置 `terminal = 1` 并照常返回该窗口，下一次调用才结束循环。真实 `Correct.cpp`
实测：`read=250/window=100 → [0,99] [100,199] [200,249]`，`read=300` 则是
`[0,99] [100,199] [200,299]` 三条而不是四条。

`determine_overlap_region` 的黄金值（`thr=15, window=405, target=2000`）：

```
start=5    -> extra_begin=10 extra_end=0  y_start=0    o_len=395
start=0    -> extra_begin=15 extra_end=0  y_start=0    o_len=390
start=50   -> extra_begin=0  extra_end=0  y_start=35   o_len=405
start=500  -> extra_begin=0  extra_end=0  y_start=485  o_len=405
start=1595 -> extra_begin=0  extra_end=0  y_start=1580 o_len=405
start=1640 -> extra_begin=0  extra_end=30 y_start=1625 o_len=375
```

注意 `extra_end` 在起点变负之前就已按 `window - o_len` 算好，所以
`extra_begin + o_len + extra_end` 恒等于窗口长度；`y-start` 越界的三类输入
（`< 0`、`>= 长度`、距离末端过近）都返回拒绝。

### P4 带状编辑距离（`alignment::banded_bpm`）

`Reserve_Banded_BPM` 是 Myers 位并行算法在 `2*threshold+1` 宽的 band 上滑动：
`Word = uint64_t`，每步 `D0 = ((VP + (X&VP)) ^ VP) | X`（加法按无符号回绕），
`D0` 最低位为 0 时错误数加一；band 每步右移四个碱基掩码并从
`pattern[2*threshold + i+1]` 补入最高位；最后沿 `VP`/`VN` 走
`p_length - t_length` 个候选终点取最小错误。

**最容易译错的一点**：等式表是 256 槽的，初始 band 会按 `pattern[i]` 写入任意
字节的槽，但紧接着 `memset(Peq, 0, ...)` 只恢复 A/C/G/T——所以初始 band 里的
`N` 位被抹掉；而**滑动阶段**的 `Peq[pattern[i_bd]] |= Mask` 写在 memset 之后，
`N` 位会留下且不参与右移。Rust 版保留完整 256 槽表并复刻这一区别，10 组
手工用例覆盖 `N` 在初始 band / 滑动 band / 仅 query 侧三种位置。

**四路 SSE 的等价性**：`verify_window` 实际调用
`Reserve_Banded_BPM_4_SSE_only`（`Word_32` × 4 lane）。用真实头文件做了
4000 组随机用例（含 `N`、阈值 1..8、长度 8..48）逐 lane 对比，
**SSE 与 4 次标量调用的 `(end_site, error)` 全部一致（mismatches=0）**，
因此 Rust 版用标量实现即可，无需引入 SIMD 也能保证输出一致。另取其中
40 组固定种子的随机序列作为常驻黄金测试。

### P4 窗口验证（`correction::verify`）

`y_start_offset` 是 f_cigar 上的阶跃函数：命中最后一个 gap site 时直接取它的
shift，否则取"最后一个 `pos <= x_start`"的 shift；落在记录范围外时 C 版会
`exit(0)`，Rust 版返回 `None` 并跳过该 overlap。对拍真实 `Hash_Table.h`
得到：sites `(0,0) (100,5) (200,-3) (300,8)` 下
`x=50→0`、`100→5`、`150→5`、`300→8`，`x=-1/350` 均为越界。

`verify_window` 把已验证的部件串起来：满窗 overlap 用固定 `THRESHOLD=15`
和 `Window_Len=375+30=405`，部分窗 overlap 按覆盖长度算
`threshold = floor(x_len * max_ov_diff_ec)` 再 `Adjust_Threshold` 并压到 31；
每成功一次就累加 `align_length` 并 `append_window_list`。C 版四人一组的
SSE 调用在 Rust 中改为四次标量调用（等价性见上），顺序与原组一致。

### P4 回溯与 cigar 生成（`alignment::banded_path` / `correction::cigar_generator`）

`Reserve_Banded_BPM_PATH` 是标量版的超集：除逐步保存每个 `text` 位置的
`D0/VP/VN/HP/HN` 五列外，还先走两条捷径——`old_error == 0` 直接返回（此时
`path_length` 在 C 里**不写**，`generate_cigar` 靠 `error==0` 特判），或
`try_cigar` 证明差异全为替换时复用旧结果。回溯时按
`delta/h_value/v_value` 三路取最小值决定方向（0 match / 1 mismatch / 2 多 y /
3 多 x），`current_value` 归零后剩余位补 0 并把 `start_site` 相应前移。

15 组黄金值覆盖：完全匹配、1/2 个错配、真实随机序列、`300000000`
（一个插入）、以及两个阈值内无解的失败用例。Rust 版逐项比对
`end_site / error / start_site / path` 的每一位。

`generate_cigar` 的三个易错点已锁定：

1. **path 是从后往前存的**：`path[0]` 是最后一步，游程编码同样倒序扫描；
2. **两端错配裁剪**：`path[0..]` 连续的 `1` 改成 `3` 并递减 `end`，
   `path[末尾..]` 连续的 `1` 改成 `3` 并递增 `start`；`trem_p` 记录裁剪边界，
   之后贪心挪 gap 的扫描只走到 `trem_p+1`；
3. **cell 打包**沿用 `op<<14 | len` 并写回窗口的 `cidx`/`clen`。

黄金值：

```
exact           -> cells=[(0,10)]                    path=0000000000
one_mismatch    -> cells=[(0,8),(1,1),(0,1)]         path=0100000000
two_leading     -> cells=[(0,8),(3,2)]  end=9        path=3300000000
middle_mismatch -> cells=[(0,12),(1,1),(0,9)]        path=0000000001000000000000
indel           -> cells=[(0,8),(3,1)]                path=300000000
```

`get_init_err_thres` / `double_error_threshold` 的黄金值
（`block_s=375, block_err=15, rate=0.04`）：`len=375→15/31`、
`100→4/8`、`10→1/2`、`3→0/0`、`4→1/2`、`300→12/31`、`1→0/0`。
注意 `double_error_threshold` 对 `x_len>=300` 有一个"先抬到 31"的下限。

`append_unmatched_wins` 按窗口 id 把已验证窗口散到对应槽，其余槽填哨兵
（`y_start=y_end=-1, error=error_threshold=-1`）；真实 C 探针的四种分布
（只有中间窗、中间+末尾、首尾、全部命中）已逐槽比对。C 版断言要求输入窗口
边界严格落在 `block_s` 网格上，越界输入本身非法。

### P4 边界修正与一致性回溯（`correction::boundary` / `correction::dag_con` / `correction::consensus_graph`）

`fix_boundary` 只在 `real_y_start == 0` 或 `end_site == Window_Len-1` 时尝试移窗，
且 `old_extra_begin/old_extra_end != 0` 时直接放弃；新结果必须
`new_error < old_error` 才接受。黄金值（`target=2000A, x_len=10, thr=2`）：

```
shift_left          -> total=3   start=2 end=11 error=0 path_len=10
right_boundary      -> total=502 start=2 end=11 error=0 path_len=10
shift_left_no_gain  -> 拒绝（old_error=0）
left_edge_no_shift  -> 拒绝（移窗后 y_start 未变）
with_extra_begin    -> 拒绝
away_from_boundary  -> 拒绝
```

`Merge_DAGCon` 把 `Edge.length` 当作 visited 标记复用（`Visit(E)` 宏），
这一点在 Rust 里保留同一个字段而不是另设 bool；`delete_Node_DAGCon` 的
边移除必须在**节点标记为 'D' 之前**完成，否则 `get_bi_Edge` 会因
`If_Node_Exist` 拒绝而留下悬空边——这是 `shared_prefix` 用例抓到的真实差异。

`build_DAGCon` 黄金值：

```
bubbles GG(3), TT(1)  -> support=3 seq=GG  nodes=6
bubbles GG(3)         -> support=3 seq=GG  nodes=4
bubbles G(5),T(2),C(1)-> support=5 seq=G   nodes=5
bubbles GG(4), GA(1)  -> support=5 seq=GG  nodes=6   ← 共享前缀合并后权重相加
```

`get_seq_from_Graph` 在 `max_count >= total_count*0.60` 时纠错，否则若位点
落在严格同聚物里用 `0.515` 阈值，再否则直接沿用 backbone 碱基。黄金值：

```
plain            -> ACGT   corrected=0  cigar=[(0,4)]
alternative(9)   -> ATGT   corrected=1  cigar=[(0,1),(1,1),(0,2)]
weak_alt(1)      -> ACGT   corrected=0  cigar=[(0,4)]
insertion(9)     -> AGGCGT corrected=2  cigar=[(0,1),(2,2),(0,3)]
```

`generate_consensus` 单轮主循环对 1/100/375/376/750/1000 长度、无 overlap
的读段输出原序列（与真实 `generate_consensus` 一致），窗口池最后一行只发一次、
尾部追加条件为 `window_start < read_length`。

**重要路径发现**：`ha_assemble` 里 `ha_overlap_and_correct(r)` 是**被注释掉的**
（`Assembly.cpp:2091`），实际调用的是 `ha_ec(r, ...)` → `cal_ec_r`
（`ecovlp.cpp:6268`）→ `cal_ec_multiple` → `worker_hap_ec`（`ecovlp.cpp:3234`）。
也就是说默认 HiFi 主流程的纠错走的是 `ecovlp.cpp` 的新链路
（`h_ec_lchain` / `gen_hc_r_alin_ea` / `rphase_hc` / `wcns_gen` /
`push_nec_re` / `sl_ec_r`），而 `Correct.cpp` 的 `worker_ovec` /
`correct_overlap` 属于旧路径（调试开关与 UL 分支仍会用到）。

因此 `--write-ec` 的逐字节一致性**必须以 `ecovlp.cpp` 链路为准**；
已完成的 `Correct.cpp` 模块是全量翻译的必要组成部分，也覆盖该链路调用的
底层原语（带状 BPM、POA、cigar），但不能单独满足门禁。下一阶段需按
`worker_hap_ec` 的调用顺序补齐 `h_ec_lchain` → `gen_hc_r_alin_ea` →
`rphase_hc` → `dedup_chains` → `wcns_gen` → `push_nec_re` →
`sl_ec_r`（后者把纠错序列写回 `R_INF.read_sperate`）。

**已知缺口**：`process_boundary` 已移植并接入单轮循环，但 `ecovlp.cpp`
新链路尚未移植；`--write-ec` 门禁还不能开启。

### P4 端到端门禁实测（早期 `pipeline` 路径，已被取代）

```
rshifiasm -t1 --write-ec --write-paf -o /tmp/rs_gate chr11-2M.fa.gz
```

早期实现在 `assembly::pipeline` 里直接跑单轮纠错，约 19 秒写出 `*.ec.fa`
（81 MB）与 `*.ovlp.paf`，但纠错结果与 C++ 有差异（首个差异在 `*.ec.fa`
第 5160 字节）。该路径现已由 `correction::ec_driver` 的多轮
`ha_ec` 链路取代，最终门禁结果见后文
“P4 `--write-ec` 门禁实测（2026-11-06，三轮逐字节一致）”。
输出器格式仍由单测锁定（`>name\nseq\n` 与 PAF 的 12 列 + `255`）。

### P4 二轮纠错与读段写回

`merge_cigars` 的真实黄金值（50 碱基窗口，安全中段 `[10, 40)`）：

```
tmp cells (0,10)(1,30)(0,10) -> tmp=[(1,30)] obtained=40 result=[(0,10),(1,30)]
tmp cells (0,50)            -> 纯 match，不需合并
tmp cells (0,10)(2,2)(0,38) -> tmp=[(0,30)] obtained=40 result=[(0,40)]
```

`get_corrected_read_from_cigar` 的黄金值（原序列 `ACGTACGT`）：

```
match 8                       -> ACGTACGT
match 3, mismatch(10), match 4 -> ACGCACGT   （10 = new C<<3 | old G）
match 2, ins TT, match 2, del AA, match 2 -> ACTTGTGT
```

`CigarRecord` 的默认构造必须用 `current_operation = -1` 哨兵；用 `derive(Default)`
得到 0 会让第一次 `add` 误判为"接在 match 后面"并越界——这个坑由二轮纠错的
测试抓到，已改为手写 `Default`。

### P4 新链路参数（`overlap::lchain_parameters`）

`set_lchain_dp_op` 用 `expf`（单精度）与 `0.0005f` 单精度字面量计算惩罚值，
不能用 `f64::exp` 或 `0.0005` 代替。HiFi（`accurate=1, k=51`）黄金值
`gap=0.30024778842926025, skip=0.00030024780269026496`，10 组 `(accurate, k)`
组合已逐项比对。

`minimizers_qgen0` 的种子按 `rid<<33 | rev<<32 | query_position` 排序（旧链路
按 target 位置），反向链目标坐标做 `UINT32_MAX-1-(pos+1-span)` 镜像再按
`target_length - (UINT32_MAX - offset)` 还原；权重打包为 `weight<<8 | span`，
出现次数介于 `low_occ`/`high_occ` 之间时取 1、不高于 `low_occ` 取 2、高于
`high_occ` 取 `pow(1 + ceil(occ/(2*high_occ)), 1.1)`。

`comput_sc_ch_ec` 的两个易错点已被 8 组黄金用例锁定：

- **gap 惩罚的 min/max 是反的**：`dd < 4` 时取 `min(lin_pen, a_pen)`，
  `dd >= 4` 时取 `max(...)`；写成统一 `max` 会让短 gap 的分数完全不对。
- **对角差超过 16 还要过带宽闸门**：`cal_bw` 用两条读段的剩余长度算上界，
  `dd > 16 && dd > cal_bw` 直接返回 `INT32_MIN`。`gap20`（dd=20, bw=19）
  与 `gap100`（dd=100, bw=18）就是靠这条拒绝的；`gap5` 虽通过闸门但
  因 `max` 分支得 `-100`。

黄金值（`rate=0.02, gap=0.30024778842926025, skip=0.00030024780269026496`）：

```
collinear    -> bw=20 score=25
gap2         -> bw=19 score=25
gap5         -> bw=19 score=-100
gap20        -> bw=19 score=INT32_MIN
gap100       -> bw=18 score=INT32_MIN
near_end     -> bw=20 score=25
heavy_weight -> bw=19 score=5
short_span   -> bw=19 score=5
```

### P4 qgen 链化包装（`overlap::qgen_chain::qgen_mcopy_fast`）

C++ 的 `lchain_qgen_mcopy_fast`（anchor.cpp:1920）是 `h_ec_lchain` 实际调用的入口，
它负责分组、去自环、按类型限链、r485 低支持链过滤。黄金值来自链接真实
`anchor.o` 的探针（`/tmp/qgen_probe.cpp`，W1–W5 数据集）：

```
cd /tmp && g++ -O2 -I/home/keith/projects/hifiasm qgen_probe.cpp \
  $(ls /home/keith/projects/hifiasm/*.o | grep -v '/main.o$') -lz -lpthread -o qgen_probe
./qgen_probe > /tmp/qgen_probe.out
```

五个用例锁定的行为（均为 `apend_be=1, gen_off=1, mcopy_num=1, mcopy_rate=0.7,
mcopy_khit_cut=32`）：

| 用例 | 场景 | 黄金结果 |
| --- | --- | --- |
| W1 | 3 组 anchor，自身组跳过 | 2 条 overlap；目标 2 `x=[0,799] sc=101 nh=4`，目标 1 `x=[0,999] sc=126 nh=0` |
| W2 | `max_n_chain=2` + 类型 3 覆盖率窗 | 保留 A/B/C，D 被覆盖丢弃；C 的 `x=[750,839] sc=61` 过 0.7 规则 |
| W3 | r485 低支持链被 16 倍强链包含 | 只留 32-anchor 链 `sc=826 x=[0,1999]`，单 anchor 链消失 |
| W4 | 无类型达到上限 | `s[w]` 全 0，不设阈值，3 条 overlap 全保留 |
| W5 | 类型 0 超限、类型 1 未超限 | 类型 0 保留前 2 条（151/126），类型 1 的 `sc=76` 保留 |

移植时踩到的点：

- 阈值过滤块的 `for (i = 0, k = 0, lch = 0; ...)` 会**重置 `lch`**，Rust 起初
  没重置；虽因 r485 只作用于 `al < cutoff` 的保留链而不可观测，仍按 C 修正。
- 类型 3 覆盖率窗按低 32 位 `cc[i]` 累加，`cc[i] = (cwe-cws)*(max_n_chain>>1)`，
  只有整窗被覆盖才算 covered；比窗口短的链靠 `cw0 >= (cw0+cw1)*0.7` 幸存。
- r485 的门槛是 `ss_stronger >= ss_low*16` 且 `al_stronger >= al_low*16`，
  且要求低支持链的 anchor 全部落在强链 `[ob, oe)` 区间内（锚点自身跨度
  `cnt&0xff` 参与判定）。

`cargo test --release qgen_chain`：15 passed。

### P4 空洞重算与边界 cigar（`correction::non_trim` / `boundary_cigar` / `advance`）

窗口验证只处理有 overlap 覆盖的窗口，剩下的空洞由 `recalcate_window_advance`
从两侧补偿。这一段的黄金值全部来自链接真实 `Correct.o` 的探针：

```
non_trim_error_rate  A 三窗连续        0.000000000000
                     B 目标在窗口内结束 0.250000000000
                     D 中间缺一窗+20 错配 0.020000000000
                     E block=100 缺第六窗  0.000000000000
                     F 同上但 60 个错配    0.030000000000
verify_sub_window    正向 err=20 y_end=405 x_end=374 aligned=375
                     反向 err=20 y_end=32  x_end=1   aligned=375
check_if_fully_covered 6 窗全覆盖=1/0   5 条=0/0   无 overlap=0/2（`-2` 提前返回）
                       仅前两窗=0/2      中间缺口=0/1   250bp 单窗=1/0
calculate_boundary_cigars 完美  err=0  y=[399,598] cells=[(0,200)]
                          边界错配 err=10 y=[399,598] cells=[(0,96),(1,10),(0,94)]
                          邻居未比对 err=-1
```

三个容易译错的点已被用例锁定：

- `get_available_fully_covered_interval` 返回 `-2` 时 **`real_length` 不会被清零**，
  而 `check_if_fully_covered` 的 `while` 条件同时看 `flag != -2`，所以循环会立刻结束，
  无 overlap 的读段只处理第一窗并把 `abnormal` 置 2。
- `calculate_boundary_cigars` 里 `determine_overlap_region` 会**覆盖** `y_start` 与
  `extra_begin`，最终坐标必须用覆盖后的值（本次实现最初用了覆盖前的 y_start，
  导致 200bp 的边界窗口整体偏移 16bp，被黄金用例抓到）。
- `recalcate_window_advance` 的 `z->align_length` 是**累加**：`verify_window` 已经
  记过一次，重算只在此基础上增加，不要重置。

### P4 新链路比对（`correction::hc_align`）

0.25.0-r726 的默认 HiFi 纠错不再走旧 `ha_overlap_and_correct` 的窗口验证，
而是 `ecovlp.cpp` 的 `worker_hap_ec` → `cal_ec_multiple` → `align_hc_ed_post_extz`
（Correct.cpp 12951）。这一层的语义与旧路径不同：

- `bit_extz_t` 是状态机：`is_align` 只看 `err <= thre`；`clear_align` 只把
  `err` 置 `INT32_MAX`，`ps/pe/ts/te` 留给调用方；trace 宏在 `err > thre` 时
  重置后跑 DP，在 `err == 0` 时不跑 DP 直接推全 match cigar，其余情况保留
  调用方给的边界再跑 DP。
- 宏的参数顺序是 `(t_string, q_string)`，即 **pattern 是目标、text 是查询**，
  与 `ed_band_cal_semi_64_w_absent_diag` 的命名相反，不能照抄参数名翻译。
- 窗口分块 `WINDOW_HC = 775`；`init_waln` 的拒绝条件是
  `s < 0 || s >= l || l - s + 2*err + THRESHOLD_MAX_SIZE < w_l`。`+31` 的余量
  正是右边界窗口（E 用例 `thre=5`、`ee=5`）仍能通过的原因。
- 只有**反向补洞**（`aln_wlst_adv_exz(..., is_cigar=1)`）和
  `gen_backtrace_adv_exz` 往共享 cell 缓冲写 CIGAR；正向直通窗口 `clen = 0`。

黄金探针（只放 `/tmp`，链接真实 `.o`）：

```
cd /tmp && g++ -O2 -I/home/keith/projects/hifiasm ec_align_probe.cpp \
  $(ls /home/keith/projects/hifiasm/*.o | grep -v '/main.o$') -lz -lpthread -o ec_align_probe
./ec_align_probe > /tmp/ec_align_probe.out
```

探针必须先调 `init_aux_table()`：否则 `recover_UC_Read_sub_region` 从
`bit_t_seq_table` 取到全 0，所有用例都 `ret=0`。反向链用例要把
`reverse_complement(prefix+query+suffix)` 存进 `R_INF`，overlap 的 y 坐标按
RC 读给出（`F_reverse` 与 `A_exact` 逐字段相同）。

黄金表（`W x=[..] y=[..] eb ee err thre`，cell 为 `op<<14|len`）：

```
A_exact      ret=1 al=1600  三窗 err=0
B_subs       ret=1 al=1600  err=6/8/2（每 100bp 一个替换）
C_ins3       ret=1          W1 err=3（目标多 3bp）
D_left_edge  ret=1          W0 y=[0,584] eb=13
E_right_edge ret=1 al=700   W1 y=[1170,1299] ee=5 thre=5
F_reverse    ret=1 al=1600  与 A 相同
G_del60      ret=0          60bp 缺失，无窗口保留
H_diverged   ret=0          qovlp 拒绝
I_forced     ret=1 wn=0     force_aln 跳过 qovlp，但无窗通过
K_bt         ret=0 al=775   W1 cells = 225E 3M 1I 2E 1D 544E
L_midfail    ret=0 al=825   W1 cells = 250E（反向补洞）
J_shift      ret=1          W2 y=[2143,2402]（f_cigar 位移在第三窗生效）
```

四个容易译错的点已被用例锁定：

- **直通窗口的 `y_start` 是 `init_waln` 之后的 `t_s`**，不是 overlap 的
  `y_pos_s + offset`；`y_end = y_start + pe`，非回溯 DP 不会给出 `ps`。
- **`push_hc_wlst_exz` 在推入当前窗口前会反转 `windows[a_n..]` 的尾部**，
  反向补洞是按从右到左的顺序生成的；K/L 各只有一个反向窗口，多窗顺序
  按源码实现，若后续出现回归应补一组双反向窗口的探针。
- **`recal_boundary_exz` 左支（`ps == 0`）永远返回空**：它要设的 `ts = toff`
  与调用方传入的 `t_pri_l` 相同，`determine_overlap_region` 结果不变；只有
  右支（`pe + 1 == target.length`）会真正改变映射，且只接受 `err < err0`。
- **`pass_qovlp(o, a, r) = (a > 0) && (o * r <= a)`**，`o` 是查询跨度
  `x_end + 1 - x_start`，不是对齐长度。

Rust 侧单测：`cargo test --release hc_align`，用例
`matches_the_reference_window_records` 复刻探针的 LCG（同一组常数）与全部
13 个构造，逐窗口断言 `x/y/eb/ee/err/thre` 与 cell 序列，并断言两轮延伸误差
（见下一节）。

### P4 延伸误差（`correction::hc_align::gen_extend_err_exz`）

`gen_hc_r_alin` 在 `align_hc_ed_post_extz` 成功之后调用
`gen_extend_err_exz`（Correct.cpp 13400），把未对齐块的误差也计入整条
overlap 的错误率；`ga_hc_r_alin` 之后用这个 `rr` 与 `e_rate` 比较。这一层
包含四个函数：`ed_cut`（13063）、`gen_extend_err_0_exz`（13224）和
`gen_extend_err_exz`，以及被 `ed_cut` 直接调用的反向延伸内核。

**反向端点约定（最容易译错）**：`Reserve_Banded_BPM_Extension_REV`
（`Levenshtein_distance.h` 216）内部把两条序列反转后调用正向内核，返回
`p_end = p_length - return_site - 1`、`t_end = t_length - i - 1`；而
`alignment_extension` 包装（374）把结果换算成 `p_length - p_end` /
`t_length - t_end`。`ed_cut` 用的是 **REV 函数本身**，因此
`aln_qlen = ql - q_end` 中的 `q_end` 是“距左端”的坐标（0 表示反向比对
从最右端开始）。Rust 侧为此单独实现了
`reserve_banded_bpm_extension_reverse`，把两端映射回 `length - 1 - end`
后复用正向内核，再原样返回 REV 约定，**不能**直接复用
`alignment_extension`。

`gen_extend_err_0_exz` 的邻窗锚定：

- `tb[0]` 只在 `pk > 0 && qs == w[pk].x_end + 1` 时取
  `w[pk].y_end + 1`；若该窗 `clen == 0`，先 `gen_backtrace_adv_exz`
  原地重算（会改写 `y_start/y_end` 并补上 CIGAR），**再用更新后的
  `y_end`**。
- `tb[1]` 只在 `pk + 1 < an && qe + 1 == w[pk+1].x_start` 时取
  `w[pk+1].y_start - ql`；同样在 `clen == 0` 时先回填 CIGAR。
- 两端都缺失时 `tb[0] = tb[1] = ts`，
  `ts = (qs - x_pos_s) + y_pos_s + y_start_offset(qs, f_cigar)`。
- 两端都成功且 `al[0] + al[1] > ql` 时按 `ql / (al0 + al1)` 缩放误差，
  C 的 `(int64_t)` 转换是向零截断；Rust 的 `as i64` 等价（结果为有限值）。

`gen_extend_err_exz` 从右向左遍历窗口，每两个窗口之间用
`get_win_id_by_s` 按块长补齐；`e_max > 0 && tot_e > ovl * e_max` 时立即
返回 `DBL_MAX` / `INT64_MAX`（调用方据此跳过该 overlap）。所有块补齐后
`assert(tot_l == ovl)`，所以探针只有在完整覆盖时才跑第二轮。

黄金探针 `/tmp/ec_extend_probe.cpp` 由比对探针复制而来，增加
`gen_extend_err_exz` 的手写原型（该函数没有头文件声明），在
`ret || w_list.n > 0` 时调用两轮：第一轮 `e_max = e_rate*1.5 + 1e-6`，
第二轮 `e_max = 0`（关闭提前中止，走完整条 overlap），每轮后打印窗口表
以观察补洞 backtrace 的原地改写：

```
cd /tmp && g++ -O2 -I/home/keith/projects/hifiasm ec_extend_probe.cpp \
  $(ls /home/keith/projects/hifiasm/*.o | grep -v '/main.o$') -lz -lpthread -o ec_extend_probe
./ec_extend_probe > /tmp/ec_extend_probe.out
```

黄金表（`rr` 为 `re/1600` 的 IEEE double 除法，E 用例的分母是 700）：

```
A_exact      early 0 / 0           full 0 / 0
B_subs       early 0.01 / 16       full 0.01 / 16
C_ins3       early 3/1600 / 3      full 3/1600 / 3
D_left_edge  early 0 / 0           full 0 / 0      （eb=13 不计入）
E_right_edge early 0 / 0           full 0 / 0
F_reverse    early 0 / 0           full 0 / 0
I_forced     early DBL_MAX/INT64_MAX  full 0.3125 / 500（wn=0，整条 overlap 都是尾巴）
K_bt         early DBL_MAX/INT64_MAX  full 0.21875 / 350
L_midfail    early DBL_MAX/INT64_MAX  full 527/1600 / 527
M_bt_gap     early DBL_MAX/INT64_MAX  full 163/1600 / 163
J_shift      early 3/1600 / 3      full 3/1600 / 3
```

新增的 `M_bt_gap` 用例专门覆盖“补洞时回溯直通窗口”分支：前两窗直通对齐、
第三窗因目标 2200..2399 被打乱而失败；第一轮补洞对 `w[1]` 做 backtrace，
把它从 `y_start=1344, clen=0` 改写成 `y_start=1375, clen=1, 775E`（与
`A_exact` 的直通窗口不同，`y_start` 是实际对齐起点）。Rust 侧用
`backtrace_window_in_place` 先克隆再写回，避免与 `region.window_cigar`
的可变借用冲突；单测在 `final_windows` 上锁定这次改写。

`gen_extend_err_exz` 的 Rust 单测用第二轮 `e_max = 0` 的完整数据校验
`ed_cut` 的正/反向、`gen_extend_err_0_exz` 的双端/单端/全缺失三种组合以及
向零截断；`e_rate * 1.5 + 0.000001` 保持与 C 相同的两次浮点运算。

### P4 exz 比对层（`correction::cal_exz`）

`hc_aln_exz_adv_hc` 依赖的坐标映射与阈值尝试层。每个函数都有对应的 C
黄金值，探针 `/tmp/ec_calexz_probe.cpp` 直接链接未修改的 `Correct.o`：

```
cd /tmp && g++ -O2 -I/home/keith/projects/hifiasm ec_calexz_probe.cpp \
  $(ls /home/keith/projects/hifiasm/*.o | grep -v '/main.o$') -lz -lpthread -o ec_calexz_probe
./ec_calexz_probe > /home/keith/projects/rshifiasm/reference/ec_calexz_probe.out
cargo test --release probe_goldens -- --ignored
```

黄金覆盖 `INITWALN`（288 组坐标）、`ADJUST`（8 组）、以及四个场景
（精确目标、每 97bp 替换、6bp 插入、反向链）的 `SEMI`/`EST`/`ESTHC`/
`EXACT`/`INFIA` 记录，共 381 行；Rust 回放逐字节一致。

容易译错的细节：

- `init_waln` 的 `aux_end` 在 `aux_beg` 缩短区间**之前**计算，因此 `r_s < 0`
  时 `aux_end != w_l - r_l`；`WindowAlignment::auxiliary_end` 保留了这个值。
- `scale_ed_thre` 的 `err << 1` 是 **32 位**运算，先回绕再转 64 位；带宽
  向上取整到 64 的倍数后 `(w-1)>>1`，最后钳到 `maximum_error`。
- `cal_exact_exz` 在 `memcmp` **之前**就写入 `*pts/*pte`，失败时也不回滚；
  成功时 `ps/pe/ts/te` 先写成区间内相对坐标再加偏移。它不检查 `q_tot_l`。
- `cal_exz_infi_adv` 先 `clear_align`（只清 `err`），坐标映射用
  `min(thre, max(ql,tl))` 预钳制，映射后再把 `thre` 钳到 `max(ql,tl)`；
  `thre <= *pthre` 时直接返回，**`exz.thre` 与 CIGAR 保持上次调用的值**。
  调用成功时 `*pthre` 已被更新，失败也保留更新值。
- 两个误差估计器把完整窗口按整数累加，部分窗口按
  `(int64_t)((double)total + error*overlap/window_len)` 逐步截断，未覆盖的
  碱基再按 `e_rate` 计费；`_hc` 版本额外要求窗口误差全零、左端
  `q[0]==t[0]`、目标连续拼接且首尾坐标完全吻合。
- `push_replace_alnw_adv` 依赖 `x_start == qs` / `x_end == qe+1` 断言
  （参考构建未定义 `NDEBUG`），Rust 侧同样用 `assert_eq!` 保留。

### P4 比对驱动（`correction::aln_adv`）

`hc_aln_exz_adv` / `hc_aln_exz_adv_hc` 用固定的阈值阶梯尝试一个链块：
估计值 → `ql*e_rate` → 上一档的两倍 → `ql*0.51` → 短块（`ql<=512`）直接
用 `maxe`。探针与回放：

```
cd /tmp && g++ -O2 -I/home/keith/projects/hifiasm ec_alnadv_probe.cpp \
  $(ls /home/keith/projects/hifiasm/*.o | grep -v '/main.o$') -lz -lpthread -o ec_alnadv_probe
./ec_alnadv_probe > /home/keith/projects/rshifiasm/reference/ec_alnadv_probe.out
cargo test --release aln_adv -- --ignored
```

黄金覆盖四个场景 × 16 组参数（两条驱动、`est` 为 -1/0/5/100、`maxl`
限长、`maxe` 限错、前/后向延伸、半全局、`ql<=16` 精确小块），182 行逐字节
一致。两个容易忽略的点：

- 入口的 `(te-ts) <= 0` 检查发生在 `mode` 改成 3 **之后**，所以
  `ts == te == -1` 的“半全局”调用总是直接返回 0，只执行 `clear_align`
  （`err` 置 `INT32_MAX`、`thre` 置 0），CIGAR/其他字段保持上一次调用的值。
  Rust 回放的 `A_exact hc i=6`、`B_subs hc i=6`、`pl i=15` 锁定这一点。
- `push_alnw` 不写对齐窗口的 `extra_begin`/`extra_end`（参考里是 malloc
  垃圾值），探针在打印前把它们归零；这两个字段只对未对齐窗口有意义。

### P4 整段比对（`correction::ovlp_base`）

`hc_ovlp_base_direct` 把一条链按相邻种子切成块，逐块调用
`hc_aln_exz_adv_hc`；失败块落成未比对窗口。探针与回放：

```
cd /tmp && g++ -O2 -I/home/keith/projects/hifiasm ec_ovlpbase_probe.cpp \
  $(ls /home/keith/projects/hifiasm/*.o | grep -v '/main.o$') -lz -lpthread -o ec_ovlpbase_probe
./ec_ovlpbase_probe > /home/keith/projects/rshifiasm/reference/ec_ovlpbase_probe.out
cargo test --release ovlp_base -- --ignored
```

十组场景共 25 次调用、70 行逐字节一致：精确窗口捷径（含平移、非零 error、
`non_homopolymer_errors!=0`、目标长度不足）、`pre_mode` 为
-1/0/1/2/其它值时首尾块的取舍、空链的 mode 3、超长块失败、反向链、以及
零迭代调用下 `extz` 的透传。三个容易忽略的点：

- 捷径要求 `pre_mode < 0 && non_homopolymer_errors == 0`，且窗口列表首尾
  外扩后 `q[0]`/`q[1]` 必须仍是原窗口端点；平铺判定比较的是相邻窗口的
  `y_end` 差值而非 `y_start`，所以跨度不等的窗口也可能通过。
- `l < 0` 的块沿用 `q[0]=0` 但 `t[0]` 保持 `-1`，随后 `adjust_ext_offset`
  以 mode 2 向左扩展；`i == ch_n` 的块 `q[1]=ql`、`t[1]` 保持 `-1`，以
  mode 1 向右扩展。两者失败时 `push_unmap_alnw` 记录的是 `t[0]` 与
  `t[1]-1`（不是 `t[0]-1`）。
- 当 `pre_mode` 为 0/1/2 且链过短时循环一次都不执行，函数不改 `extz`，
  `aux_o` 也不追加窗口；回放的 `H_no_iteration` 锁定这一透传行为。

### P4 窗口重链（`correction::window_chain::rechain_alignment_hc`）

`gen_hc_fast_cigar0`（Correct.cpp:17813）在首轮整段比对后，对每个未比对窗口调用
`rechain_aln_hc` 再试一次。黄金值来自链接真实 `Correct.o` 的探针：

```
cd /tmp && g++ -O2 -I/home/keith/projects/hifiasm ec_chain_probe.cpp \
  $(ls /home/keith/projects/hifiasm/*.o | grep -v '/main.o$') -lz -lpthread -o ec_chain_probe
./ec_chain_probe > /home/keith/projects/rshifiasm/reference/ec_chain_probe.out
cargo test --release window_chain -- --ignored
```

黄金文件 691 行，sha256
`28c7ca0540d84248d58471676da3e8da38347399aca370ae08a901d79f29a398`。
覆盖 `generate_window_chain`（GW/CH 段，11 组）与三组合成读段
（精确平铺、错配、反向）× 四组重试规格（mode 0/1/2，含目标整体平移 -7）。
被锁定的语义：

- **`z` 是原始 overlap 区域，不是首轮输出。** 首轮比对把窗口写进 `aux_o`，
  `z->w_list` 仍是 `clen == 0` 的平铺窗口；重试时 `gen_win_chain` 对这些窗口
  走 `gen_backtrace_adv_exz`，并从 `z->y_id` / `y_pos_strand` 取目标读段。
  把首轮输出当作 `z` 会让 C 读出目标读段末尾（未初始化/越界内存），Rust 的
  `sub_region` 则截断补 `A`，两者必然分叉。
- `rechain_aln_hc` 从 `aux_o->w_list.a[aux_i]` 取 `qs/qe/ts/te`（`qe/te` 加 1），
  任一段短于 `FORCE_SIN_L`（512）直接返回；`mode` 取自该窗口的
  `error_threshold`，不是 -1。
- 链化后按 mode 断言首/尾锚点，并且 `mode 0` 需要至少 3 个锚点、
  `mode 1/2` 至少 2 个才真正重比对。
- 只有重比对恰好新增一个且非未比对窗口，且新旧窗口在 mode 要求的端点上完全
  相同时，才用新窗口覆盖 `aux_o->w_list.a[aux_i]` 并弹出末尾窗口。
- `cl->length` 必须在返回前恢复为进入时的值（`must reset!!!!`），否则
  `gen_hc_fast_cigar0` 的后续窗口会从错误偏移取链。

已知差异：`recover_UC_Read_sub_region` 不检查读段边界，越界区间在 C 里读到的是
压缩堆的遗留字节；Rust 的 `PackedSequence::sub_region` 会截断并用 `A` 补齐。
本组黄金用例的区间都在目标读段内，不触发该差异；真实数据若出现越界窗口，
两者的输出无法逐字节对齐（C 侧本身依赖未定义行为）。

### P4 快速比对驱动（`correction::fast_cigar::generate_fast_cigar_direct`）

`gen_hc_fast_cigar0`（Correct.cpp:17813）是 `worker_hap_ec` 里每个 overlap 的
实际比对入口：按链索引取同一条目标读的种子组，整段比对，再对未比对窗口逐个
重试。黄金值来自链接真实 `.o` 的独立探针：

```
cd /tmp && g++ -O2 -I/home/keith/projects/hifiasm ec_fast_cigar_probe.cpp \
  $(ls /home/keith/projects/hifiasm/*.o | grep -v '/main.o$') -lz -lpthread -o ec_fast_cigar_probe
./ec_fast_cigar_probe > /home/keith/projects/rshifiasm/reference/ec_fast_cigar_probe.out
cargo test --release fast_cigar -- --ignored
```

黄金文件 216 行，sha256
`d25d76fdc1882cd967e536132052e8f5dc02cce810c503fd02c9feb69e2c558f`。
覆盖六组合成读段（精确平铺、错配平移、反向链、长读精确/错配、长读大段损坏）
× 四组规格（plain / 前缀链索引 / `h_khit=1000` / `nhpe`）。

被锁定的语义：

- `ch_idx = z->shared_seed`，`ch_n` 是 `cl->list` 中从 `ch_idx` 起 `readID`
  相同的一段；`ch_n <= 0` 直接返回 0，不写 `*re`。
- 首轮 `hc_ovlp_base_direct(..., pre_mode = -1)` 写入 `aux_o->w_list`；
  `aux_o` 从 `z` 复制 `y_id` / `y_pos_strand` 与四个端点坐标，窗口清空。
- 重试循环只处理 `is_ualn_win`（`error == INT16_MAX && clen == 0 &&
  extra_end < 0`）的窗口；`rechain_aln_hc` 会从 `z` 取平铺窗口、把新窗口
  追加到 `aux_o->w_list`，并恢复 `cl->length`。
- 只有 `aux_o->w_list.n > aux_n` 时才整理：丢弃仍未比对的初始窗口、压缩其余
  窗口，然后 `radix_sort_window_list_xs_srt` 按 `x_start` 稳定排序。
- `update_overlap_region(z, aux_o, ql, tl)` 之后 `z` 的窗口列表即最终结果；
  `*re` 是最终窗口的误差总和，未比对窗口按 `max(x 长度, y 长度)` 计。

已走到的关键路径：首轮精确捷径（`W1_exact`）、重试成功恢复（`L_mixed_*`）、
首轮与重试都失败后窗口被丢弃（`X_mixed_corrupt` 的 `[4800,8944]`）、前缀链
索引（`prefix` 变体打印 `ch_idx=3 length=8`）。

### P4 空洞重排（`correction::gap_reassign`）

`reassign_gaps`（Correct.cpp:25409）在 `gen_hc_fast_cigar` 之后运行：同聚物与
串联重复里 indel 的位置不唯一，比对器可能把 gap 放在任意等价位置并留下
"错配" cell。该函数把每个 gap 尽量左移，沿途把碱基其实相同的 mismatch cell
改回 match，并相应减少 `non_homopolymer_errors` 与窗口 `error`。

```
cd /tmp && g++ -O2 -I/home/keith/projects/hifiasm ec_reassign_probe.cpp \
  $(ls /home/keith/projects/hifiasm/*.o | grep -v '/main.o$') -lz -lpthread -o ec_reassign_probe
./ec_reassign_probe > /home/keith/projects/rshifiasm/reference/ec_reassign_probe.out
cargo test --release gap_reassign -- --ignored
```

黄金文件 62 行，sha256
`c7836fedd263e103d192ce459c18b4c0c8b4758091111eea5be18095b46dfba2`。
15 组场景：9 组真实 overlap（先经 `generate_fast_cigar` 对齐，再重排；含精确、
置换、indel、反向链、低复杂度同聚物/串联重复插入删除、`nhpe == 0` 提前返回）
+ 6 组手写 window cigar（`manual_match_recovery` / `manual_mismatch_break` /
`manual_insertion` / `manual_no_move` / `manual_unaligned` /
`manual_unaligned_move`），逐窗口打印 `cidx`/`clen` 与每个 cell。

被锁定的语义：

- `non_homopolymer_errors == 0` 直接返回，`aux` 保持原样（`no_errors` 场景）。
- 未比对窗口（`error == INT16_MAX && clen == 0 && extra_end < 0`）整结构复制到
  `aux`，不触碰 cell 缓冲；`error == 0` 的窗口整段复制并 `push_wcigar`。
- 没有"match 之后的 gap"时整段复制并 `ajust_end_cigar`：首尾连续的 mismatch
  cell 改成 `more query`（op 3），同时 `y_start += len` / `y_end -= len`。
- `adjust_gap` 从窗口尾部向前扫描：遇到 gap cell 就把新 gap 放在它前面；
  `op0 == 2`（目标多碱基）先 `query_position--`，`op0 == 3` 先
  `target_position--`；mismatch 分支对整段逐个碱基比较，只有
  `l[1] > 0`（匹配位之后还有未匹配位）才把 `ff` 置 1 并计入 `rd_err`。
- `ff == 0` 时 gap 原样追加（`manual_no_move`）；`ff != 0` 时按扫描停下的
  cell 截断后倒序回放缓冲；扫描走完全部 cell 则从窗口开头重新拼。
- 任一窗口移动（`rr`）后 `update_overlap_region` 用 `aux` 替换 `z` 的窗口并
  重算四个端点；`z->non_homopolymer_errors -= re`（uint32 回绕减法）。

已知限制：C 的 `kv_pushp` 不初始化新窗口，`move_wins` 又没有写
`extra_begin`/`extra_end`，因此这两个字段在 C 里是未初始化内存。本移植写 0
（与 C 实际看到的全新缓冲区值一致），黄金探针在打印前也会把非未比对窗口的
这两个字段清零。

### P4 整段精化（`correction::hc_r_alin`）

`gen_hc_r_alin`（Correct.cpp:25617）是 `worker_hap_ec` 里一条读段的 overlap
列表的总入口：逐个 overlap 依次调用 `align_hc_ed_post_extz` →
`gen_extend_err_exz` → `gen_hc_fast_cigar`，再经过两个末端过滤器
`ff_lunalign`（Correct.cpp:25554）与 `ff_tend`（Correct.cpp:25466），最后
`reassign_gaps`。黄金值来自链接真实 `.o` 的独立探针：

```
cd /tmp && g++ -O2 -I/home/keith/projects/hifiasm ec_hc_ralin_probe.cpp \
  $(ls /home/keith/projects/hifiasm/*.o | grep -v '/main.o$') -lz -lpthread -o ec_hc_ralin_probe
./ec_hc_ralin_probe > /home/keith/projects/rshifiasm/reference/ec_hc_ralin_probe.out
cargo test --release whole_overlap_refinement -- --ignored
```

黄金文件 111 行，sha256
`8a843cea62b253266d3251e330569d1a4a8295f874bb61450029a4e69dd97854`。
24 组单读段批次覆盖：精确、周期错配、插入、缺失、大段缺失、块状损坏（100bp /
150bp）、同聚物插入/删除、高分歧、中部损坏、反向链、`wl=400`、链索引前缀
（`filler_hits=3`）、ONT 参数（`wl=375`、`chem_drop=1`、`gap_rate=0.006`、
`gap_max=64`）、`align_gap_max=0`，以及 5 条混合 overlap 的压缩；另有 6 组
直接调用 `ff_tend` 的手写窗口场景（`clean` / `left_bad` / `left_bad_low` /
`right_bad` / `middle_bad` / `unaligned_first`）。

被锁定的语义：

- `z->shared_seed = z->non_homopolymer_errors` 在进入对齐前写入，作为
  `return_t_chain` 的链索引；`non_homopolymer_errors` 随后被
  `gen_extend_err_exz` 的 `re` 覆盖。
- `gen_extend_err_exz` 的 `maximum_error_rate = e_rate * 1.5 + 0.000001`，
  `rr > e_rate` 时丢弃；`is_match` 在调用之后无条件清 0。
- `gen_hc_fast_cigar` 恒返回 1，其 `re` 输出被调用方丢弃；存活 overlap 的
  `non_homopolymer_errors` 仍来自上一步。
- `align_gap_max >= 0` 才运行 `ff_lunalign`；`chem_drop` 为真才运行
  `ff_tend`，其错误率上界是 `min(e_rate*10, 0.36)`。
- 存活者按 `swap(list[k], list[i])` 就地压缩到列表前部，随后
  `is_match = 1`、`strong = 0`、`without_large_indel = 0`；失败的 overlap
  留在尾部，顺序就是交换后的排列（`mix` 场景：保留 y1/y3/y5，尾部为 y4/y2）。
  Rust 版返回新的长度、不截断切片，调用方再收缩列表。
- `ff_tend` 的前向/后向扫描共用 `mm`/`err` 累加器，`pop_trace` 会把相邻同
  操作 cell 合并；`op == 2`（目标多碱基）不推进 query 坐标。

已知限制：`ff_lunalign` 的深路径（非单窗口、存在未比对空洞）在合成数据上难以
到达——同一处损坏会先被 `align_hc_ed_post_extz` 的 `qovlp` 检查拒绝，探针的
`gapblock` 系列因此停在 `alen=775`。单窗口精确平铺捷径与常规窗口链已覆盖；
深路径留待 ONT 数据（P8）实测时用真实读段锁定。`ff_tend` 是 `inline` 函数，
`Correct.o` 里仍有弱符号，探针直接调用该符号取得黄金值。

### P4 精确重叠回收（`gen_hc_r_alin_ea` / `gen_hc_r_alin_nec`）

`gen_hc_r_alin_ea`（ecovlp.cpp:2810）是 `worker_hap_ec` 中 `gen_hc_r_alin` 的
替代入口：它先尝试用上一轮已经存下的**精确子区间**（`ma_hit_t.el == 1`）直接
恢复 overlap，不跑任何比对；剩下的交给 `gen_hc_r_alin_nec`
（Correct.cpp:25677）——即只在 `is_match != 1` 时才做整段精化的版本。黄金值
来自链接真实 `.o` 的独立探针：

```
cd /tmp && g++ -O2 -I/home/keith/projects/hifiasm ec_hc_ralin_ea_probe.cpp \
  $(ls /home/keith/projects/hifiasm/*.o | grep -v '/main.o$') -lz -lpthread -o ec_hc_ralin_ea_probe
./ec_hc_ralin_ea_probe > /home/keith/projects/rshifiasm/reference/ec_hc_ralin_ea_probe.out
cargo test --release exact_overlap_recovery -- --ignored
```

黄金文件 106 行，sha256
`17732ec5f83caf31a34eca174df75e16478aa86577c252e9c22de79cb97a70de`。
18 组场景：`ea_exact_prefix` / `ea_fallback_prefix`（各带 3 条 filler 命中，
验证链索引转存）、`ea_empty`（空 overlap 表）、`ea_no_el`（无 `el` 命中）、
`ea_exact`、`ea_exact_wrong_seq`、`ea_coord_mismatch`、`ea_rev_exact`、
`ea_strand_mismatch`、`ea_two_hits_first_mismatch` /
`ea_two_hits_first_match`、`ea_two_targets`、`ea_mixed`（精确 +
存活 + 淘汰）、`ea_nec_filters`（ONT 参数）、`ea_exact_merge` /
`ea_exact_no_merge`（窗口列表已有的旧窗口），以及直接调用 nec 的
`nec_mixed` / `nec_all_prematched`。

被锁定的语义：

- `srt` 只收 `el == 1` 的命中，键为 `(tn<<1 | rev)<<32 | k`；当前 overlap
  的 `oi` 键为 `(y_id<<1 | y_strand)<<32 | k`，两者都用
  `radix_sort_ec64`（MSD、每趟 8 位，与 `ksort::radix_sort_u64` 一致）排序。
  `oi` 的构建同时把每个 overlap 的 `is_match` / `strong` /
  `without_large_indel` 清零。
- 查询游标只在新键**小于**当前目标键时前进；命中后**不前进**，因此同一个
  目标键的多个 stored 命中只有排序后的第一条会被查询，且重复的 overlap 会
  复用同一条命中（`ea_two_hits_first_*` 两组正反顺序）。查不到键时不做任何事。
- 坐标相等按存储侧的**半开区间**比较：`x_pos_s == qns`、
  `x_pos_e + 1 == qe`、`y_pos_s == ts`、`y_pos_e + 1 == te`。反向命中
  （`rev = 1`）的目标子区间在 RC 坐标系中取
  （`recover_UC_Read_sub_region(tu->seq, ts, te-ts, trev, ...)`），与链上的
  `y` 坐标同一框架；`ea_rev_exact` 用 `query` 与 `reverse(exact)` 锁定。
- `exact_ec_check` 只比较长度与 `memcmp`（`ql`/`tl` 参数未用）；Rust 版在
  切片越界时判为不匹配（C 会越界读，但真实 `paf` 坐标不会产生这种输入）。
- 命中成功后：`is_match = 1`，`shared_seed = non_homopolymer_errors`（链索引
  转存），`non_homopolymer_errors = 0`，`strong = without_large_indel = 0`；
  `set_exact_exz` 用半开端点建一条 exact 记录，`push_alnw` 把它追加到窗口
  列表，`nec++`。该路径下 **`push_alnw` 的合并分支不可达**：exact 窗口从
  `x_pos_s` 开始，而能与之合并的旧窗口必须结束于 `x_pos_s - 1`；探针的
  `ea_exact_merge` / `ea_exact_no_merge` 两组都得到 `n=2` 即是证据，合并规则
  本身已由 `cal_exz` 的黄金用例覆盖。
- 全部命中时 `ol->length` 保持不变（返回 `on`）；`on > nec` 时调用
  `gen_hc_r_alin_nec`，后者把列表压缩到 `k` 并返回。`srt->n` 在函数入口先
  清零再填充，空 overlap 表则在清零**之前**返回（保留调用方原值）。
- `gen_hc_r_alin_nec` 与 `gen_hc_r_alin` 的唯一差别是整段对齐/
  `gen_extend_err_exz` / `gen_hc_fast_cigar` / `ff_lunalign` / `ff_tend` /
  `reassign_gaps` 全部位于 `if (z->is_match != 1)` 内；`swap` 压缩、
  `is_match = 1`、`strong = 0`、`without_large_indel = 0` 和 `k++` 对每个
  `i` 无条件执行，已匹配的条目因此原地通过（`nec_mixed` 的 y=4 保留原窗口，
  `nec_all_prematched` 全部保留）。`el` 列表为空时走的是 `gen_hc_r_alin`
  而非 nec。

### P4 已存 overlap 辅助（`dedup_chains` / `push_nec_re` / `extract_max_exact` / `push_ne_ovlp` / `check_well_cal`）

`worker_hap_ec` 在共识之后用这组函数把结果写回两条已存 overlap 表：
`push_nec_re` 取出本轮共识 cigar，`push_ne_ovlp` 把 `is_match == flag` 的
overlap 写进 `ma_hit_t`（前向表带 `el` 精确子区间，反向表不带），
`check_well_cal` 判断读段是否完全纠正、覆盖是否异常，`dedup_chains` 在整段
精化后按目标读去重。黄金值来自链接真实 `.o` 的独立探针：

```
cd /tmp && g++ -O2 -I/home/keith/projects/hifiasm ec_stored_overlap_probe.cpp \
  $(ls /home/keith/projects/hifiasm/*.o | grep -v '/main.o$') -lz -lpthread -o ec_stored_overlap_probe
./ec_stored_overlap_probe > /home/keith/projects/rshifiasm/reference/ec_stored_overlap_probe.out
cargo test --release stored_overlap_helpers -- --ignored
```

黄金文件 227 行，sha256
`4aabcd7258766ad713f5ae21fff26bf11b9a81e55204afa51aba7bef5e8f0a2d`；
其中 78 行是 `dedup_chains` 压缩后参考仍打印的**旧缓冲区尾部**（见下），
回放只比对余下 149 行。22 组场景：

- `dedup_mixed`（同目标多条、`is_match`/得分/跨度优先级）、`dedup_ties`
  （得分与跨度都相等，`>` 严格比较保留先遇到者）、`dedup_equal`（4 条全等，
  走插入排序）、`dedup_many`（70 条全等，越过 `RS_MIN_SIZE` 走 MSD radix，
  验证等键稳定）、`dedup_single`、`dedup_empty`。
- `exact_full`、`exact_two_windows`（游标跨窗口续扫）、`exact_with_indel`、
  `exact_none`（全 `(u32)-1` 哨兵）、`exact_short_second`、`exact_offset`。
- `push_forward`（`flag=1`、带共识 cigar）、`push_reverse`（`flag=2`、
  `ec == NULL`）。
- `well_full_double`、`well_interior_hole`（`abnormal=1`）、`well_left_gap` /
  `well_right_gap`（`abnormal=2`）、`well_low_depth`（`f_ec=0`）、
  `well_empty`、`well_non_match_cell`（consensus 含非 match cell）、
  `well_short_read`。

被锁定的语义：

- `dedup_chains` 用 `overlap_region_sort_y_id`（4 字节键、从 bit 24 起的
  MSD radix，低于 `RS_MIN_SIZE` 时插入排序；两条路径都稳定），随后
  `k = 1, l = m = 0` 扫描分组；组内多于一条时按优先级选一条：`is_match`
  更小、`sc = (x_pos_e + 1 - x_pos_s) - non_homopolymer_errors * 12` 更大
  （`plus`/`minus` 按 `u32` 回绕后再转 `i64`）、跨度更大；`mm_sc` 初值
  `INT32_MIN`。单条组保留自身。选中者与位置 `m` 交换后 `m++`；
  `ol->length = m`。Rust 版排序的是索引排列再用同一键重建列表，结果与
  直接排序记录一致；列表用 `Vec::truncate` 收缩，因此参考缓冲区里被交换到
  尾部的旧条目无法复现——回放测试从黄金中过滤 `TAIL` 行。
- `push_nec_re` 断言辅助 overlap 恰有一个窗口，按 `cidx/clen` 原样复制
  cell。
- `extract_max_exact` 逐窗口工作：窗口切片上的 `pop_trace`（0x3fff 掩码，
  相邻同操作 cell 合并）；`wx[1] <= wx[0]` 或跨度不大于当前最优时跳过；
  只有 `op == 0` 的 run 才调用 `extract_max_exact_sub`。sub 先按游标
  `exk/eyk/eck` 回退（`ck > 0 && xk >= xs0`），回退时按 op 取长度掩码
  （D/I 为 0xfff、M 为 0x3ff、E 为 0x3fff），`op == 2` 不移动 query、
  `op == 3` 不移动 target；再前向扫描到 `xk >= xe0`，用 `pop_trace_bp_f`
  （D/I 0xfff、M 0x3ff 且 bq/bt 分别在 bit 10-11/12-13、E 0x3fff），按
  op+两个碱基合并；仅 `op == 0` 参与，求 `[wx0,wx1)` 与 `[xs0,xe0)` 的
  交集，严格长于当前最优才更新 query（`wy0 + os - wx0`）与 target
  （`ys0 + os - xs0`）。`ye0` 参数从未被读取。游标跨窗口、跨调用持续；
  每次调用内部的最优区从零开始，外层跨窗口保留最优。无结果时四个输出
  `(u32)-1`，Rust 用 `Option` 表达。
- `push_ne_ovlp` 只收 `is_match == flag`：`qns = (x_id << 32) | x_pos_s`，
  `qe = x_pos_e + 1`，`ts = y_pos_s`，`te = y_pos_e + 1`，`rev` 取
  `y_pos_strand`，`bl` 取目标读长度，`ml = strong`（C 里是 1 bit 字段），
  `no_l_indel` 原样。`ec` 非空时调用 `extract_max_exact`，先置 `el = 0`，
  找到子区间则重写四个坐标并置 `el = 1`。`ec == NULL`（反向表）时 C
  **从不写 `el`**，读到的是堆上的旧值（探针首次运行得到 184）；Rust 版在
  该路径置 `exact_region = false`，且 `reverse_paf.el` 在参考里从未被读取，
  探针因此只在带 `ec` 的场景打印 `el`。
- `check_well_cal` 取 `s = (uint32_t)qns`、`e = qe`，压入 `s << 1` 与
  `(e << 1) | 1`，用 `radix_sort_ec64`（即 `radix_sort_u64`）排序后扫描
  深度（奇数键 `--`、偶数键 `++`）；每个 `ed > st` 区间：`old_dp < min_dp`
  则 `f_ec = 0`；`old_dp == 0` 时，`st > 0 && ed < len` 置 `abnormal = 1`，
  否则仅在 `abnormal == 0` 时置 2。最后一段 `[st, len)` 重复同样的判定。
  若仍 `f_ec`，再要求 consensus cigar 每个 cell 的 op 都为 0。

已知差异：`dedup_chains` 压缩后的尾部旧条目（参考打印、Rust 截断）；
`ec == NULL` 时 `el` 的不确定值（参考从未消费）；Rust 用 `PafHit` 承载
`ma_hit_t` 的已消费字段。

### P4 单读段单倍型划分（`correction::rphase` / `correction::haplotype_evidence`）

`worker_hap_ec` 的相位环节：`rphase_hc` 把查询读段切成滑窗，逐窗口解码相交
overlap 的窗口 cigar，统计每个查询位置上与其不一致的 overlap 数；过阈值的
位置成为候选位点。`extract_sub_cigar_hc` 两遍扫描（`set_f=1` 统计 flag、
`set_f=0` 产出证据），`push_info` 把同一位点的证据折叠成每个碱基的等位统计，
`generate_haplotypes_naive_HiFi` 决定哪些 overlap 支持第二单倍型并更新
`strong` 与 `mapped_overlaps_length`。黄金值来自链接真实 `.o` 的独立探针：

```
cd /tmp && g++ -O2 -I/home/keith/projects/hifiasm ec_rphase_hc_probe.cpp \
  $(ls /home/keith/projects/hifiasm/*.o | grep -v '/main.o$') -lz -lpthread -o ec_rphase_hc_probe
./ec_rphase_hc_probe > /home/keith/projects/rshifiasm/reference/ec_rphase_hc_probe.out
cargo test --release correction::rphase::probe_goldens -- --ignored
```

黄金文件 153 行、7819 字节，sha256
`bc1abd60076f7ccb392d6353429407954729b2963c0d1cb19a9b8b16ddca5a5c`；
确定性输出，Rust 回放逐行一致。六组场景，均为 HiFi 配置
（`bd=0, wl=WINDOW_MAX_SIZE, occ_thres=1, hpc_len=0, std_bs=0, lindel=0`）：

- `strong_and_weak`：强等位（occ_0/occ_1 都达标）与弱等位混排，一条反向链
  overlap；验证 `is_match 1→2`、`strong=1`、`mapped -= 360`、match 统计的
  `occ_0` 递减（7→4）。
- `multi_check_weak`：四个弱位点，`multi_check` 的 32bp 间距过滤与
  `snp_srt.n -= o` 只弹 scratch 段；验证 `score=-1` 且 `aim=1`、`mapped=0`。
- `two_windows`：两个滑窗、候选压缩、只在第二窗出现的候选；验证
  `snp_srt` 为空、`mapped -= 1500`。
- `indel_shift`：`D`/`I` cell 改变目标坐标；验证 `t-xk+yk` 投影与
  `overlapSite` 原值存入 `cov`。
- `adjacent_sites`：site 20/21 相邻被整组丢弃，保留组的
  `overlapSite -= m_off`（`m_off=2`）；验证 `id=2` 保留原压入下标。
- `strong_reset`：未晋级的 overlap 其 `score=1` 被清回 -1，最终不再晋级；
  验证 `mapped -= 240`。

被锁定的语义：

- `WINDOW_MAX_SIZE = WINDOW + (int)(1.0 / HA_MIN_OV_DIFF) + 3`
  `= 375 + 50 + 3 = 428`。
- 候选键 `(x_start+bd) << 32 | c_idx 下标`，一次 64 位 radix；同起点分组内
  改写为 `(x_end+1-bd) << 32 | 下标` 后组内再排。滑窗 `s=0..ql` 步长 `wl`；
  上一轮 `hc_phase_robust_rr` 返回 1（有候选 `q1 <= e`）时压缩候选表，只留
  与本窗相交者；新增候选指针 `i` 跨窗口前进、`q0 >= e` 即停。
- 第一遍只统计 `op==1` cell 覆盖的查询位置，flag 饱和加 1 到 127；阈值
  `> occ_thres` 者 `nn_snp++`、flag 置 1；首尾清零区间 `si..ei` 在第二遍
  **之后**才 memset，因此第二遍仍能看到 1。
- 第二遍对 flag 非零位置产证据：`op==0` 取查询碱基；`op==1` 延迟
  `recover_UC_Read_sub_region`，缓冲区起点由**首个**错配的
  `t-xk+yk-hpc_len`（截 0）确定、长度 `yk1-yk0`，同一批证据共用该起点
  （参考的 `yk0` 只赋值一次）；`om = (flag==3)`，`type=(om<<1)+op`。
- 游标两次调用之间必须一致（参考 assert）；`set_f` 结束时把
  `xk0/yk0/ck0` 复位到回退后的位置，`cur_ylen` 记录
  `min(yk1+hpc_len, 读长) - yk0`。
- `push_info`：4 字节 radix 按 `overlapID`；同 ID 折叠取该段排序后的第一条
  （依赖 radix 置换含等键顺序）；`occ_2` 按 u32 回绕；`occ_0==0 || diff<=1`
  或（有 overlap 表时）`rev_n==occ_0` 直接返回 0；每个 `occ_1[i] >= 2` 的
  碱基压一条 `SnpStats`（`occ_0 = 1+occ_0` 截 u32、
  `occ_2 - occ_0 - occ_1` u64 回绕后截 u32、`score=-1`、`id` 为压入时
  下标）；压缩证据时 match 指向该位点最后一条统计，mismatch 把原
  `overlapSite` 存入 `cov`（32 位）再指向等位 id。
- `generate_haplotypes_naive_HiFi`：先丢弃与相邻位点（`site-1` 或
  `site+1`）同组的统计，`overlapSite -= (l - m_snp_stat)`；`snp_srt` 键
  `((u32::MAX - o) << 32) | l`，`o` 为强等位 mismatch 数（`occ_0>=2`、
  `occ_1>=2`、非 `is_st_bs`、`occ_0>=s_hap_cov`、`occ_1>=infor_cov`）；
  第一遍只处理有强等位的 overlap：`is_match 1→2`、mismatch 统计
  `score=1`、match 统计沿同 site 段向下把 `occ_0`（以及 `std_bs` 下的
  `overlap_num`）减去 `cov`（断言 >=1）；第二遍按 `score==1` 再晋级一批；
  仅在 `snp_srt` 非空时把未晋级 overlap 的 score 清回 -1。
- multi_check：跳过 `is_match==2`，收集弱位点（非强等位、`score!=1`）的
  `overlapSite`；`snp_srt.n -= o` 只弹掉 scratch 段；`o >= align_length*up`
  时对 scratch 段排序做 32bp 间距过滤（`i==0` 用前一轮遗留的 `s`、末位用
  遗留的 `t`，复现参考的 stale 指针）；`z>=2` 才保留。随后整表排序做重复
  计数，但 `l = k` 在 if 外，重复计数恒不触发——保持死代码。
- 最终分类：`is_match==2` 置 `strong=1` 并扣 `x_pos_e+1-x_pos_s`（u32
  回绕后转 i64）；`is_match==1` 扫证据，`score==1 && occ_0>=2 && occ_1>=2
  && !is_st_bs` 置 `strong=1`，mismatch 再升为 2、扣跨度并 break。
- `hpc_mask_ff` 在 HiFi 路径从不被调用（`hpc_len=0`）；仍按参考实现覆盖
  `r=1..4` 的四种扫描（含 `> round` 与 `>= round*cutoff` 双条件）。

已知差异 / 未接线：

- `rphase_lidel`、`generate_haplotypes_sv` 未移植；`lindel != 0` 显式
  assert 拒绝（HiFi 调用方传 0）。
- `hpc_mask_ff` 的 flag 清零分支省略：该路径所有调用传 NULL。
- `haplotype_evdience_alloc` 的 `snp_matrix` / `r_snp` / `dp` 未移植
  （HiFi 路径不用）。
- 尚未接入 `pipeline.rs`，`--write-ec` 门禁保持不变；下一步 `wcns_gen`
  后再由 `worker_hap_ec`/`cal_ec_r` 驱动接线。

### P4 窗口一致序列（`correction::consensus_gfa`）

`worker_hap_ec` 的 `wcns_gen`：按 `wl=512` 的滑窗对查询逐位投票，覆盖不足的
区间回退到小 POA 图（`cns_gfa`）的一致序列；结果作为一条 base-aware 的
`aux_o` 窗口 cigar 交给下一步。投票、短空洞 `cns_gen0`、长空洞
`cns_gen_full`（按 `cns_g_wl=256` 切窗）、`refine_cns_g` 图化简、
`gseq_cns_g` 最重路径、`push_correct0_fhc` / `push_correct1_fhc`（含
`indel_exz` 重写）全部按参考逐位对齐。黄金值来自链接真实 `.o` 的探针：

```
cd /tmp && g++ -O2 -I/home/keith/projects/hifiasm ec_wcns_probe.cpp \
  $(ls /home/keith/projects/hifiasm/*.o | grep -v '/main.o$') -lz -lpthread -o ec_wcns_probe
./ec_wcns_probe > /home/keith/projects/rshifiasm/reference/ec_wcns_probe.out
cargo test --release correction::consensus_gfa::probe_goldens -- --ignored
```

黄金文件 700 行、96178 字节，sha256
`f922bd27e1f34a0548be5c4aac683495478c5132a15c635cbc56faec654b3c37`；
20 组场景逐行一致（`NEC`、每个 overlap 的 `without_large_indel`、`aux_o`
窗口边界与 cell 序列）。参数固定为 HiFi 调用方实参
`bd=0, wl=512, occ_tot=3, occ_exact=0.500001, cns_g_wl=256, max_trace=ql,
rid=0`。

场景覆盖：多窗口全匹配 / 短查询 / 无候选；少数与多数错配串（`cns_gen0`
成功与失败）；空锚点插入（`s==e`）与 `cns_gen0` 失败后的空图路径；长错配
空洞与多子窗图；所有位点失败时的退化尾窗口；unaligned 窗口
（`without_large_indel=0`）；跳过 `is_match!=1` 与零窗口 overlap；部分窗口；
缺失多数；尾部短/长空洞；反向链读段；插入后缺失的图路径。

被锁定的语义：

- 候选键 `(x_start+bd) << 32 | c_idx 下标`，一次 64 位 radix；同起点分组内
  改写为 `(x_end+1-bd) << 32 | 下标` 后组内再排。`ii_a` 与 `ii_b` 共享同一
  `c_idx` 数组，只是 `idx` 不同；Rust 用单一候选表即可。
- 投票计数器：每站点 `ct[k<<1]` 高 32 位记 op0，低 32 位记覆盖
  （op0/op1/op3）；`ct[(k<<1)+1]` 插入计数对 op0 内部碱基加高+低、对
  op1/op3 只加低、对 op2 的 `ws` 只加低；op0 串首碱基只加覆盖。通过条件
  `oc0 > oc1*occ_exact && oc0 > oc1-oc0 && oc1 >= occ_tot && oc0 > 1`。
- 低覆盖（`oc1 < occ_tot`）进入锚点分支并默认参考碱基；插入位点先推前一个
  锚点，再以当前位置开新锚点。
- `cns_gen0`：`e > s+6` 直接失败；`extract_sub_cigar_ii` 要求
  `(s0<iws)||(s0==0)` 与 `(e0>iwe)||(e0==ql)`；阈值与多数表排序后各查一次，
  失败置 `idx->ru=1`。
- `cns_gen_full` 以 `cns_g_wl` 切窗，`s_end=(s==s0)`；`cns_gen_full0` 在
  `iter_cc_idx_t` 之后马上把 `rr` 清 0，循环内按 `q1 <= e` 重新置 1。
- `push_correct1_fhc_indel_exz`：回溯跳过非 op0 cell 累计 `ok/nk/e0/f`；
  `!f`/`!ok`/`!nk` 时直接用图碱基；否则 `diff=min(31,ml)` 调
  `cal_exz_global`，失败再用 `diff=min(63,ml)` 重试；`exz->err < e0` 时从
  `cn0=ck+1` 起按 exz cigar 重写并截断 `sc->n`（**不刷新 `clen`**），否则按
  `ml < e0` 走较短者回退或图碱基。
- 尾锚点 `(ql, ql, is_tail=1)` 即使无任何锚点也会产生半开窗口
  （如 `xs=600, xe=599`）；截断只缩短共享 cell 缓冲，`clen` 保持旧值，下一次
  push 仍按旧 `clen>0` 尝试合并。

黄金探针发现并修正的两处图分区错误（此前单测未覆盖）：

- `trace_cns_bp` 沿节点**出边** `[0, nou)` 找同碱基后继；Rust 曾误扫入边。
- `merge_cns_g_in` 统计 `w` 的反向弧时扫 `w` 的**出边**；Rust 曾误扫入边，
  在空图插入场景触发 `target == current` 断言。

已知差异 / 未接线：

- `get_sub_seq` 越界请求在参考中返回未初始化 scratch 字节，Rust 补 `A`；
  仅在参考会读未初始化内存的输入上不同。
- `cns_gen_full0` 的 `occ_tot`/`occ_max` 参数在参考中未使用，Rust 保留形参
  但不起作用。
- 尚未接入 `pipeline.rs`；下一步把 `recover_UC_Read`、`h_ec_lchain`、
  `fetch_aux_ovlp`、`gen_hc_r_alin_ea`、`rphase_hc`、`dedup_chains`、
  `wcns_gen`、`push_nec_re`、`push_ne_ovlp`、`check_well_cal` 与
  `refresh_ec_ovec_buf_t0` 串成真实 `worker_hap_ec`，再由
  `cal_ec_r`/`cal_ec_multiple`/`cal_sec_ec_multiple`/`cal_update_ec_multiple`
  驱动到 `--write-ec` 逐字节门禁。

### P4 `--write-ec` 门禁实测（2026-11-06，三轮逐字节一致）

```
cargo build --release
./target/release/rshifiasm -t1 -o /tmp/p4fast --write-paf --write-ec \
    /home/keith/Downloads/chr11-2M.fa.gz
cmp /tmp/p4fast.ec.fa reference/write_paf_ec.ec.fa
# 84,247,384 B，无差异
```

正式门禁同样通过：

```
cargo test --release --test golden_pipeline \
    corrected_reads_are_byte_identical_to_the_reference -- --ignored
# test result: ok. 1 passed; 0 failed ... finished in 153.83s
```

逐轮日志与 C++ 基准完全一致（`# bases` / `# corrected bases` /
`# exact o` / `# non-exact o` / counting 三类计数 / peak 值）：

| 轮次 | bases / corrected | exact / non-exact | collected / indexed / distinct | peak_hom / peak_het |
| --- | --- | --- | --- | --- |
| 1 | 84042542 / 180149 | 231097 / 18161 | 2264901 / 2224445 / 65289 | 36 / 18 |
| 2 | 84042541 / 2444 | 248616 / 6741 | 2264594 / 2259905 / 64530 | 53 / 36 |
| 3 | 84042686 / 275 | 250558 / 4824 | 2264589 / 2260311 / 64524 | 36 / 18 |

耗时（单线程，同一数据）：Rust `Real time: 162.873 sec`（pec 70.201 / 34.554 /
32.834）对 C++ `109.106 sec`（pec 42.923 / 23.694 / 22.493），约 1.5×；
Peak RSS 7.531 GB（C++ 16.48 GB，bloom filter 惰性分配所致，不影响结果）。

关键修复：`PackedSequence::sub_region` 原先每次调用都 `unpack()` 整条读段再切片，
而 `hc_ovlp_base_direct` 的每个窗口、每次重试都会调用它；快照读段长度下它是
O(read_length)，纠正一轮从 ~70 s 被放大到 ~994 s。现按
`recover_UC_Read_sub_region` 直接解码所请求区间（O(length)，反向链按
`bit_t_seq_table_rc` 向下走，N 位点事后覆盖，越界位置解码为零字节的 `A`），
新单测 `sub_region_matches_the_whole_read_unpack_at_every_offset` 覆盖
全部起点/长度/两条链与整读解包切片对拍；修复前后两次全量运行的 `.ec.fa`
互相比对也逐字节一致。

调试钩子（仅定位用，默认关闭）：

```
RSHIFIASM_EC_PROFILE=1     每轮每读段打印 seeds/chain/refine/phase/dedup/consensus 耗时
RSHIFIASM_EC_STOP_AT=<n>   只纠前 n 条读段
```

### P4 已知差异与后续

- `--write-paf` 产物已逐字节一致：`chr11-2M.fa.gz` 上
  `.ec.fa`（84,247,384 B）与 `.ovlp.paf`（28,934,587 B）均与 C++ 完全相同
  （见 `tests/golden_pipeline.rs` 第二条测试）。
- ONT/超长读段模式（`-l`）与部分纠错（`--write-partially-corrected-reads`）
  未接入，当前只接受 HiFi 默认链路。
- 日志中的时间/CPU/内存字段不参与逐字节比对（仅格式对齐），GFA/FASTA/PAF
  属于强一致范围。

### P5 门禁实测（2026-11-06，r_utg/p_utg/p_ctg 逐字节一致）

- **探针级**（`tests/bp_graph_golden.rs`，黄金值来自 `/tmp/bp_probe.cpp`
  链接真实 `.o` 的转储）：
  - `output_bp_graph` 主图半边：`probe.bp.p_ctg.gfa` / `.p_ctg.noseq.gfa` /
    `.p_ctg.lowQ.bed` 逐字节一致（`adjust_utg_by_primary` 的 purge/tip/bubble
    清理 + `print_utg` 写盘）；
  - trio 转移表：`clean_u_trans_t_idx` 前后（`f15`/`f16`）逐条目一致，
    `set_trio_flag_by_cov` 产出的 trio 标志（`f17_trio_flag.bin`）逐字节一致；
  - 相位图：`ma_ug_gen_phase`（`/tmp/trioprobe`）逐节点一致。
- **端到端**：`./target/release/rshifiasm -t1 -o reference/candidate_default
  chr11-2M.fa.gz`（171.6 s，峰值 7.5 GB），`scripts/compare_outputs.sh`
  对 `p_ctg`/`p_utg`/`r_utg` 的 gfa、noseq.gfa、lowQ.bed 共 9 个文件
  `cmp` 无差异。
- 运行在预期断点停下：`output_bp_graph` 写完 `bp.p_ctg.*` 后返回
  `output_trio_graph_joint is not implemented yet`（P6 待移植）。
- Rust 版当时不写 `*.ovlp.source.bin` / `*.ovlp.reverse.bin` / `*.ec.bin`
  重跑缓存索引；重启缓存现已在选项空间轮次实现（见下文），`compare_outputs.sh`
  对 `ovlp` 两个 bin 报 identical，对 `ec.bin` 报 DIFFERS 属预期
  （C++ 侧非确定性 padding，以 `read_cache_golden` 分节比对为准）。

### P6 门禁实测（2026-11-07，hap1/hap2 逐字节一致）

移植范围（`Overlaps.cpp`）：

- 清理函数族（`src/assembly/trio_clean.rs`）：`if_primary_unitig`、
  `set_drop_trio_flag`、`update_unitig_graph`、`magic_trio_phasing`、
  `get_candidate_uids`、`unitig_simi`、`all_to_all_deduplicate`、
  `is_useful_node`/`recover_chain_nodes`/`rescue_useless_trio_nodes`、
  `delete_useless_trio_nodes`、`update_hap_label`、`force_trio_clean`、
  `clean_trio_untig_graph`（redo 循环 + `FlySimp` 相位缓存跨轮复用）。
- 驱动层（`src/assembly/bp_graph.rs`）：`adjust_utg_by_trio`、
  `output_trio_unitig_graph`（is_bench 路径）、`DedupIndex`
  （`gen_dedup_idx_t` 的侵入式链表布局）、`cal_exact_ug_o`、
  `dedup_exact_ug`、`push_ma_utg_t`、`append_miss_nid`、
  `output_hap_graph`、`output_trio_graph_joint`。
- 默认链路关键分支：`self_scaf=0` 走 `output_hap_graph`；
  `purge_level_trio=0` 跳过 purge_dups；`HA_F_BAN_POST_JOIN` 默认关，
  两个 rescue 都执行；`break_ug_contig`/`update_dump_trio` 默认不触发。

门禁结果：

- **探针级**（`tests/bp_graph_golden.rs`）：`output_bp_graph` 全流程（含
  hap 输出）对 `/tmp/bp_probe.cpp`（f 系列转储）与
  `reference/default_repeat.bp.*` 逐字节一致——p_ctg/hap1/hap2 的
  gfa、noseq.gfa、lowQ.bed 共 9 个文件。
- **端到端**：`rshifiasm -t1`（159.4 s，峰值 7.5 GB，exit 0），
  `scripts/compare_outputs.sh` 对 r_utg/p_utg/p_ctg/hap1/hap2 的
  gfa、noseq.gfa、lowQ.bed 共 15 个文件 `cmp` 无差异。
- 日志核对：`[M::output_trio_graph_joint] dedup_base::0, miss_base::0`
  与 C++ flag 探针输出一致（该语料无精确重复、无遗漏 unitig）。

移植中踩过的坑（对后续阶段有用）：

- **HG 标签来自活数组**：`ma_ug_print2` 的 `HG:A:` 取自
  `"apmaaa"[trio_flag]`（FATHER→`p`、MOTHER→`m`、其余→`a`），而 C++ 的
  `R_INF.trio_flag` 是被各阶段持续改写的同一个数组。Rust 侧打印原先读
  `ReadStore::trio_labels`（初始副本），而变异写在独立的 `trio_flags`
  克隆上；p_ctg 在 `set_trio_flag_by_cov` **之前**打印（全 `a`）掩盖了
  这一问题，hap 输出才暴露。修复：打印链显式传入活的 `trio_labels`
  切片（`ma_ug_print2`/`ma_ug_print`/`ma_ug_print_simple` 及其调用方）。
- `f17_trio_flag.bin` 带 4 字节长度头，逐字节对拍时不要按文件偏移直接
  索引读段 id（本项目测试读取已正确跳过）。
- 定位工具：`/tmp/flag_probe.cpp`（加载 d20 转储，逐步打印
  `R_INF.trio_flag` 分布），其 hap 输出与正式参考逐字节一致。

### 选项空间门禁实测（2026-11-07，--primary / -e / 重启缓存）

#### `--primary` 分支（`output_contig_graph_primary` + `output_contig_graph_alternative`）

- 非 partition 分发（`clean_graph` 末尾 else，`Overlaps.cpp` 39659）：
  primary 图在**活的读图**上跑 `adjust_utg_by_primary`（读段类别被改写、
  未认领读段被 drop），alternative 图再从改写后的读图生成
  `ma_ug_gen_primary(ALTER_LABLE)`——两个驱动的先后依赖与 C++ 相同。
- 门禁：`hifiasm --primary -t1` vs Rust `--primary -t1`——
  `p_ctg`/`a_ctg`/`r_utg`/`p_utg` 各 gfa/noseq/lowQ.bed 全部逐字节一致
  （165.9 s，exit 0）。输出前缀无 `.bp` 段（`get_outfile_name` 在
  `HA_F_PARTITION` 清除后返回原名），与 C++ 一致。
- `break_ug_contig`（`--b-cov/--h-cov`）仍为显式 Err（默认不触发）。

#### `-e`（`HA_F_BAN_ASSEMBLY`）

- `build_string_graph_without_clean` 在缓存写盘后检查
  `HA_F_BAN_ASSEMBLY`（`Overlaps.cpp` 39714），命中即返回 Ok、exit 0、
  不产出任何 GFA。
- 门禁：归一化时间/CPU 字段后与 `hifiasm -e -t1` 日志逐行一致，
  包括写盘日志 7 行（`Writing reads to disk...` …
  `bin files have been written.`）。
- 顺带补齐了默认路径缺失的
  `Inconsistency threshold for low-quality regions in BED files: 70%`
  收尾日志行（`Overlaps.cpp` 39673）。

#### 重启缓存（`*.ec.bin` / `*.ovlp.source.bin` / `*.ovlp.reverse.bin`）

新模块 `src/assembly/read_cache.rs`，对应 `write_all_data_to_disk` /
`load_all_data_from_disk`（`Overlaps.cpp` 23567/23590）与
`write_All_reads` / `write_ma_hit_ts` / `write_ma` 的逐字段格式
（`Process_Read.cpp` 69/127，`Overlaps.cpp` 23442/23420）。

- **格式要点**：`ec.bin` 头部携带 `index_size`/`name_index_size`——这是
  C 的分配容量（`READ_INIT_NUMBER=1000` 起，`< total_reads+2` 时
  `*2+2`，5626 读 → 8014），不是条目数；`name_index` 有效条目为
  `[0..=total_reads]`（哨兵 = `total_name_length`）；每读 packed 序列
  写 `len/4+1` 字节。`ovlp.bin` 为 `i64 n_read` + 每读
  (u8 fully_corrected, u8 abnormal, u32 length) + 每 hit 42 字节流
  （u64 qns、u32 qe/tn/ts/te、u8 el/no_l_indel、u32 ml/rev/bl/del；
  `cc` 位字段不序列化）。
- **写盘门禁**：`ovlp.source/reverse.bin` 与 C++ 逐字节一致（默认、
  `--primary`、`-i`、`-e` 四个 profile 全部）；`ec.bin` 分节一致
  （`tests/read_cache_golden.rs`，掩码见下）。
- **加载门禁**：同目录第二遍运行打印
  `[M::ha_assemble::0.041*0.96] ==> loaded corrected reads and overlaps
  from disk`，跳过全部纠错（0 行 pec 日志），17 个组装产物仍与 C++
  逐字节一致，且缓存不被重写（`write = !ovlp_loaded = 0`）。
  `het_cov == -1024` 时按 `hap_recalculate_peaks`
  （`Assembly.cpp` 1881）重建峰值（本语料不触发：写入值为 36/18）。
  `-i` 时完整重算，产物一致。
- **`ec.bin` 的两处 C++ 非确定性**（字节门禁必须掩码，见
  `read_cache` 模块文档与 `tests/read_cache_golden.rs`）：
  1. 每读最后一个 packed 字节的 padding 位（`len%4==0` 时整字节）——
     C 的堆残渣，同命令重跑也不一致（实测 7–24 字节差异，全部位于
     padding）；`.ec.fa` 不受影响（只读有效碱基）。
  2. `name_index` 哨兵之后的容量槽——malloc 残渣（实测跨命令恰好
     相同，但语义上未初始化）。
  `compare_outputs.sh` 对 `ec.bin` 报 DIFFERS 属预期，以
  `read_cache_golden` 的分节比对为准；`ovlp` 两个 bin 可裸 `cmp`。
- **`total_reads_bases` 语义修复**：C++ 在载入时冻结该计数
  （`ha_insert_read_len` 之后无任何更新），而 Rust 的
  `write_corrected_read_back` 曾同步更新 `total_bases`——导致缓存头部
  写出修正后总长（84042711 vs 84042542）。已改为冻结语义（长度数组
  仍随纠错更新，与 C++ 相同）。
- 写盘时机在 `try_rescue_overlaps` **之前**（`Overlaps.cpp` 39707），
  缓存内容是 pre-rescue overlap。

#### 选项空间剩余缺口（2026-11-07 收尾后，见文末记录）

- `--ul` 现为 parse 时显式错误；`--ul-*` 参数组（`--ul-rate/--ul-tip/
  --ul-cut/--ul-round/--path-max/--path-min`）单独设置时 C++ 同样无效
  （仅 UL 路径消费），保持静默以对齐参照语义。
- `--chem-c/--chem-f`、`--rl-cut/--sc-cut`、`--kpt-rate`、
  `--integer-correct` 同上：仅 ONT/UL/trio 路径消费，单独设置无效，静默。
  Hi-C 参数组（`--seed/--s-base/--n-weight/...`）为静默 no-op（Hi-C 前端
  `--h1/--h2` 本身已是显式错误）。
- trio 前端（`Trio.cpp`，`-1/-2/-3/-4/-5/--trio-dual`）与 Hi-C 分发
  （`--h1/--h2/--enzyme`）为显式 gate 错误（需测试数据）。
- 调试项：`ha_assemble_pair`（`--sec-in`，读 `.ec0/.ec1` 旧分片缓存，
  无产出方）为桩；`--dbg-gfa/--prt-raw/-b/--dbg-het-cnt/--ex-list/
  --ex-iter` 与 `--ont/--fast` 现均为 parse 时显式错误（`check_supported`）。
- 双方一致地死掉的选项：`--skip-triobin`、`--low-het`、`--ul-m`、
  `--bin-only`（2026-11-07 核实：C++ 全库亦无消费点）。

#### 输出层收尾（2026-11-07，--b-cov/--h-cov、--telo-*、--dbg-ovec）

- **`break_ug_contig`**（`src/assembly/break_contig.rs`，`Overlaps.cpp`
  11503–12565）：覆盖深度扫描（`push_cov_interval_{direct,advance}` 的
  beg/end 键打包 + 回绕处理）、扫描线区间合并（low/high 两套合并规则）、
  谷/峰跃迁检测 + 圆形回绕、`get_break_point_idx` 的方向感知
  exact/inexact 拆分与无符号 `c_break` 取位、
  `push_sub_unitig`/`break_all_contigs` 的旧→新顶点重映射。
  `print_utg`/`output_trio_unitig_graph`/`output_contig_graph_primary`
  三处 "not ported" 守卫替换为实际调用（low/high 与 C++ 一样单开关）。
  门禁：`hifiasm --b-cov 20 -t1` 17/18 逐字节一致（p_ctg 相对默认运行
  真实断裂）。
- **端粒**（`gfa_clean.rs::gen_telo_ends`，`Overlaps.cpp` 191–347）：
  motif 轮转集合 + 每读双端两遍扫描（检测 + 轮转连续性精化），单线程化；
  `ul_clean_gfa` 穿线到三个切割函数（消费侧 `TeloEnd` 早已移植）。
  门禁：`--telo-m CCCTAA -t1` 的 `gen_telo_end_t` 摘要行逐字一致
  （0/0——chr11-2M 无端粒读段，本语料为 no-op，默认奇偶性由既有门禁
  保证）。
- **`--dbg-ovec`**（`ec_driver.rs::run_debug_overlap`，`Assembly.cpp`
  1053 / `ecovlp.cpp` 6311–6452）：无过滤索引 + 每读链化精化后直接
  落盘对齐窗口（带 cigar）——把 C++/Rust 分歧隔离在共识之前。
  门禁：125,745,896 B 的 `.ovlp.paf` 与 C++ 逐字节一致。

**introsort 差一修复（影响全局）**：klib `ks_introsort` 的深度消耗是
`if (--d == 0)`——进入分区**先减后判**；Rust 移植原为"先判后减"，comb
兜底晚一层才触发。oloist 超过深度上限（≥2^5=32 元素）时相等键元素的
排列与 C++ 不同。纠错语料的 per-read olist 从未越过上限，故 P4–P6 门禁
未暴露；dbg 的未过滤种子让 read 25 的 61 区域 olist 越线。修复后
`--dbg-ovec` 全量一致，且默认/`--write-paf --write-ec`/`--b-cov 20`
四条门禁复验全部保持逐字节一致。`lchain_qgen_mcopy_fast` 的两个
`ks_introsort_or_ss/or_xs` 调用点改为逐行字面转录的
`introsort_klib_by`（含 klib comb 排序的 gap=9/10→11 与收尾插入序）。

### P2 门禁实测（与 C++ 日志逐字节一致）

```
[M::ha_analyze_count] lowest: count[5] = 1208
[M::ha_analyze_count] highest: count[36] = 103835
[M::ha_ft_gen] peak_hom: 36; peak_het: 18
[M::ha_ct_shrink::...] ==> counted 114 distinct minimizer k-mers
[M::ha_ft_gen::...] ==> filtered out 114 k-mers occurring 180 or more times
[M::ha_opt_update_cov] updated max_n_chain to 180
```

Rust 版耗时约 11 秒（C++ 约 8.6 秒），峰值 RSS 7.5 GB（C++ 16.5 GB，
差异来自 bloom filter 的惰性分配，不影响结果）。

### P3 minimizer 门禁实测（与 C++ 日志逐字节一致）

```
[M::yak_count] collected 2264901 minimizers
[M::ha_pt_gen::...] ==> counted 105745 distinct minimizer k-mers
[M::ha_pt_gen] count[4095] = 0 (for sanity check)
[M::ha_pt_gen] peak_hom: 36; peak_het: 18
[M::ha_ct_shrink::...] ==> counted 65289 distinct minimizer k-mers
```

`mz1` sketch 的窗口环形缓冲、同分 minimizer 的重复写出、高占用 streak 的
窗口最小值标记与 `ks_heapmake`/`ks_heapdown` 采样都已逐条对移植，因此
minimizer 总数与位置分布与 C++ 完全一致。

### 位置索引门禁实测（与 C++ 日志一致）

```
[M::ha_pt_gen::] counting in normal mode
[M::ha_pt_gen::...] ==> indexed 2224445 positions, counted 65289 distinct minimizer k-mers
```

索引的分区起点按桶序累加，occ 列表按 sketch 顺序填充，两者都与
`worker_pt_gen` / `ha_pt_insert_list` 一致。

## 比对方法

```bash
# 生成基准（单线程）
scripts/reference_run.sh 1

# 运行 Rust 版
cargo run --release -- -t1 -o reference/candidate chr11-2M.fa.gz

# 逐字节比对，首个差异会打印字节偏移
scripts/compare_outputs.sh reference/reference reference/candidate
```

## 编辑距离内核（`alignment::edit_distance`）

`Levenshtein_distance.h` 是 overlap/纠错主链路真正使用的比对内核（位并行 Myers
带状编辑距离）。已完成 64 位路径的三个最常用函数：

| 移植函数 | C++ 对应 | 主链路调用次数 |
| --- | --- | --- |
| `global_64` | `ed_band_cal_global_64_w` | 2 |
| `semi_global_64` | `ed_band_cal_semi_64_w` | 1 |
| `semi_global_64_absent_diagonal` | `ed_band_cal_semi_64_w_absent_diag` | 7 |
| `extension_forward_64` | `ed_band_cal_extension_64_0_w` | 2 |
| `extension_backward_64` | `ed_band_cal_extension_64_1_w` | 2 |
| `global_64_trace` | `ed_band_cal_global_64_w_trace` | 2 |
| `semi_global_64_absent_diagonal_trace` | `ed_band_cal_semi_64_w_absent_diag_trace` | 8 |

黄金值由真实 C++ 头文件编译出的小程序产生（`gcc -msse4.2 -fgnu89-inline`），
17 组用例覆盖：完全匹配、单碱基替换、超出阈值、pattern 比 text 长/短、
`abs_diag` = 0/1/2、前向/后向延伸的 `err`/`pe`/`te`/`ps`/`ts`，结果全部一致。

### 宽度无关内核与 `_trace` 全族（2026-11-06）

`Correct.cpp:15481-15600` 的 `cal_exz_*` 按 `nword` 分派：`bd = (thre<<1)+1`，
`nword = (bd>>6) + !!(bd&63)`；`nword <= 1` 走 64 位内核，`== 2` 走 128 位，
其余走 `infi`（运行时 `nword`）。纠错新链路的阈值上界是 `MAX_SIN_E = 2047`，
因此这三档都必须存在。现已移植四个族的全部 `_trace` 内核：

| Rust 函数 | C++ 对应 | 分派档 |
| --- | --- | --- |
| `global_trace` | `ed_band_cal_global_64_w_trace` / `_128_w_trace` / `_infi_w_trace` | 全部 |
| `extension_forward_trace` | `ed_band_cal_extension_64_0_w_trace` / `extension_infi_0_w_trace` 等 | 全部 |
| `extension_backward_trace` | `ed_band_cal_extension_64_1_w_trace` / `extension_infi_1_w_trace` 等 | 全部 |
| `semi_global_absent_diagonal_trace` | `ed_band_cal_semi_64_w_absent_diag_trace` / `semi_128_w_absent_diag_trace` / `semi_infi_w_absent_diag_trace` | 全部 |

C++ 的 64/128 位内核是宏 `HA_ED_INIT` 按字宽展开的实例，`infi` 内核用运行时
`nword` 循环；Rust 用一个 `BandedVectors` 表示 5 个 pattern 向量与
`VP/VN/D0/HN/HP`，`core()` 对齐 `ed_infi_core`，`shift_pattern()` 对齐
`ed_infi_post_Peq`。逐字展开的 128/256 位实例与运行时循环在字级完全等价，
因此 664 条黄金记录（含 128 位与多字 `infi` 档）同时验证了两条路径。

对拍探针 `/tmp/ed_trace_probe.cpp` 链接**除 `main.o` 外的全部真实 `.o`**，
直接调用 `cal_exz_global/extension_0/extension_1/semi` 及其 `_simi` 变体：

```bash
cd /tmp && g++ -O2 -I/home/keith/projects/hifiasm ed_trace_probe.cpp \
    $(ls /home/keith/projects/hifiasm/*.o | grep -v '/main.o$') -lz -lpthread -o ed_trace_probe
./ed_trace_probe > /home/keith/projects/rshifiasm/reference/ed_trace_probe.out
cargo test --release trace_kernels_match_the_reference_dispatch -- --ignored
```

黄金值覆盖 10 组用例（exact、subs50、subs7、ins10、del12、shift_plus、
shift_minus、longer、shorter、tail_only）× 17 个阈值
`{1,3,16,31,32,63,64,65,100,128,255,300,400,511,512,1000,2047}` × 4 种模式，
记录 `nword err ps pe ts te pl tl cigar` 与逐 cell 的 packed CIGAR；
共 664 条 TRACE 记录，Rust 回放全部逐字节一致。

**可达性约束（semi-global）**：真实链路中 `cal_exz_infi_adv` 只在
`update_semi_coord → init_waln(thre, s, l, aln_l = ql + 2*thre, ...)` 之后调用
semi 内核，几何关系保证 `ai = pn - tn + aux_beg <= 2*thre`。超出该范围时 C 内核
会读到 `ez->nword` 之外的 `VP/VN` 字（`resize_bit_extz_t` 只增长不清理，这些字是
未初始化的堆数据），不属于可复现行为。探针与 Rust 回放测试都显式跳过
`pn - tn + aux > 2*thre` 的合成调用，并附注释说明。

`gen_trace` 的回溯状态机：每列保存 `D0/VP/VN/HP/HN` 五个位向量，回溯时按
`D`（对角）、`H`（水平）、`V`（垂直）三者比较并用 `<=` 保证「prefer indels」，
再按 `push_trace` 的 14 位长度打包规则输出 CIGAR。

尚未移植：非 `_trace` 的 `_simi` 内核（探针参数 `2` 的两阶段/复用分支）。
它们与 `_trace` 内核共享同一递推，差别只是不记录回溯列；P4 新链路实际调用的
全部是 `_trace` 版本。

## 链化（`overlap::chaining`）

`overlap::chaining::chain_dp` 移植自 `Hash_Table.cpp` 的 `chain_DP`，配套
`chain_check`（`ha_chain_check`，共线快速路径）、`chain_length`
（`get_chainLen`）与 `normal_weight`（`normal_w`）。数据结构
`KmerHit` / `ChainData` / `FakeCigar` / `OverlapRegion` 对应
`k_mer_hit` / `Chain_Data` / `Fake_Cigar` / `overlap_region`。

对拍方式与其他模块不同：直接链接**原始 `Hash_Table.o`**（连同 `kalloc.o`、
`CommandLines.o`）后调用真实 `chain_DP`，因此黄金值来自真正的目标文件而非
重新抄写的代码。4 组用例覆盖共线快路径、带 indel、锯齿、反向链：

```
collinear chain=5 xs=100 xe=500 ys=100 ye=500 ovlp=1000 seed=255 fclen=1
indels    chain=5 xs=100 xe=500 ys=100 ye=502 ovlp= 998 seed=211 fclen=4
zigzag    chain=3 xs=200 xe=400 ys=200 ye=400 ovlp=1000 seed=119 fclen=1
reverse   chain=3 xs=100 xe=560 ys=100 ye=550 ovlp= 790 seed= 57 fclen=3
```

`chain`、四个坐标、`overlapLen`、`shared_seed`、`f_cigar` 长度全部一致。
`Fake_Cigar` 的 `(site << 32) | (|shift| << 1 | sign)` 打包方式也已对齐。

### 候选种子（`overlap::candidates`）

`collect_seeds` 移植自 `anchor.cpp` 的 `ha_get_new_candidates` 前半段：

- 每个查询 minimizer 用 `ha_pt_get` 取位置列表，展开成 `anchor1_t`
- `reverse = (z->rev != y->rev)`；反向时 `self_off = rlen - 1 - (z->pos + 1 - z->span)`
- `srt = rid << 33 | rev << 32 | other_off`，先按 `srt` 做 radix sort，
  再在每个相同 `srt` 的区段内按 `self_off` 二次 radix sort
- 计数权重 `low_occ = hom_cov * 0.333`、`high_occ = hom_cov * 1.667`：
  稀有种子记 2，正常记 1，过高频按 `1 + ceil(cnt / (2*high_occ))` 再取 `^1.1`

`chain_candidate_groups` 移植自 `calculate_overlap_region_by_chaining`：
按 `(readID, strand)` 分组、跳过自比对、每组 `max_skip = 25` 调用 `chain_DP`。

种子收集与分组目前只有单测（镜像位置、计数权重、分组顺序），端到端一致性
要等 PAF 门禁；`chain_DP` 本身已与真实目标文件逐字段对拍。

### overlap 追加（`overlap::append`）

`append_inexact_overlap` 移植自 `append_inexact_overlap_region_alloc`，黄金值
同样来自链接真实 `Hash_Table.o`：4 组用例（正向/反向 × 是否带首尾 padding），
比对的字段包括镜像后的坐标、双向 strand、以及 **f_cigar 的每一个
`(pos, shift)` 条目**：

```
forward_be0 x=0..999 y=0..999      cigar=[(500,-2),(100,0)]
forward_be1 x=0..999 y=0..999      cigar=[(0,0),(500,-2),(100,0),(999,0)]
reverse_be0 x=0..599 y=200..999    cigar=[(499,0),(99,2)]
reverse_be1 x=0..599 y=200..999    cigar=[(0,0),(499,0),(99,2),(599,2)]
```

两个容易译错的点已被用例锁定：镜像时 start 取 `xLen - x_pos_e - 1`
（start/end 互换），反向链的 f_cigar 按 tmp 顺序镜像而 shift 原样保留，
正向链则倒序遍历且 shift 做 `initial_shift - gap` 变换。同 `y_id` 的链按
`shared_seed` 高者优先保留（相等时取更大 `overlapLen`）。

## 已知限制：ksw2 的内核等价性

`alignment::ksw2::extz2` 是对 `ksw_extz2_sse` 标量语义的移植：8 位回绕算术、
`max_sc` 夹紧、有符号比较与 tie-break、Z-drop、traceback 状态机都已对齐，
评分与独立的全矩阵仿射 DP 逐例一致（单测 `score_agrees_with_a_full_matrix_alignment`）。

但它**尚未**做到与 SSE 版本逐字节等价，原因是 SSE 实现按 16 条反斜线为一组
计算，会写入带外的 padding 单元，这些单元的 `s[]` 是上一轮遗留值、`qrr[]`
甚至越界读取；带内单元会读到这些 padding 值，因此要完全复现需要模拟整个
16 lane 的斜带布局（含越界读取的未定义行为）。

优先级说明：hifiasm 的 overlap/纠错主链路使用的是自带的
`Levenshtein_distance.h`（`ed_band_cal_*`），`ksw_extz2_sse` 仅在
`Correct.cpp:11449` 一处调用。因此 P3/P4 的下一步是把
`Levenshtein_distance.h` 作为主目标移植，ksw2 的 SSE 级等价性留到该路径需要时
再补。

## 已知差异与注意事项

- 日志中的时间、CPU、峰值内存数值不参与逐字节比对，仅格式对齐。
- 默认参数下 hifiasm 会写 `*.ovlp.source.bin`、`*.ovlp.reverse.bin`、`*.ec.bin`
  索引并在下次运行时加载；比对时必须使用相同的 `-i` 状态（默认加载，
  基准脚本会重新生成）。
- `--primary` 未开启时只输出 `bp.p_ctg.gfa` 系列；默认（`HA_F_PARTITION`）还输出
  `bp.hap1.p_ctg.gfa`、`bp.hap2.p_ctg.gfa`，两者都在比对范围内。

## 收尾：选项拒绝、introsort 黄金用例与死代码（2026-11-07）

本轮不触碰任何默认执行路径，只为消除"选项被解析但静默无效"的误导与
遗留杂物。门禁基线（`--write-ec`/`--write-paf`/`--primary`/`-e`/`-i` 全组
SHA-256）在收尾后复跑不变。

### `check_supported`（`command_line/parser.rs`）

`check_options` 通过后追加一轮 Rust 侧支持面检查，把"参照实现会进入未实现
分支"的开关在 parse 时显式拒绝（此前 `--ont/--fast` 在纠错驱动里才拒绝，
其余完全静默）。消息沿用现有 gate 风格（`[ERROR] ... is not implemented yet`）。

- 拒绝：`--ul`、`--dbg-gfa`、`--prt-raw`、`-b`、`--dbg-het-cnt`、
  `--ex-list`/`--ex-iter`、`--ont`（连带 `single_cell_mode`）、`--fast`。
  逐项核实过 C++ 消费点：`--dbg-gfa` 触发 `.pt_flt` 索引检查点与调试图
  dump（`htab.cpp:1144`/`Assembly.cpp:2078`），`--prt-raw` 在 `gen_init_sg`
  dump raw GFA（`Overlaps.cpp:39248`），`-b` 单读调试（`Assembly.cpp:1003`），
  `--dbg-het-cnt` 末轮 het 计数（`Assembly.cpp:1014`），`--ex-list/
  --ex-iter` checkpoint 抽取（`extract.cpp`）——这些在默认 HiFi 模式下就
  生效，Rust 字段此前赋值后无任何消费点，属真实静默分歧。
- 有意不拒绝（C++ 单独设置同样无效，报错反而偏离参照语义）：
  `--ul-rate/--ul-tip/--ul-cut/--ul-round/--path-max/--path-min`（仅 UL 路径）、
  `--chem-c/--chem-f/--rl-cut/--sc-cut`（仅 ONT/is_sc 路径）、
  `--kpt-rate`（仅 `output_trio_graph`）、`--integer-correct`（仅
  `gfa_ut.cpp:17648` 的 `ul_re_correct`）。
- 双方一致死掉的选项维持不动：`--skip-triobin`、`--low-het`、`--unskew`、
  `--ul-m`、`--bin-only`（本轮核实 C++ 全库亦无 `binning_only` 消费点）。
- `ec_driver.rs` 的运行时 ONT/单细胞/partial/`--fast` 防御检查保留。
- parse 失败退出码仍为 0，忠实 C++ `main` 行为（报错走 stderr）。

### introsort 回归测试固化

HEAD（d5f3ce7）introsort 深度 off-by-one 修复的临时对拍脚本
`tests/klib_sort_tmp.rs` 依赖 `/tmp` 黄金数据，已将两组 61 条
`(i32,u32,u32)` 记录内嵌为 `container/ksort.rs` 的黄金用例
（`introsort_klib_matches_the_c_reference_permutation`，比较器
`(second<<32)|third`，锁定相等键的参照排列），删除临时文件。

### 死代码与注释

- 删除 `src/correction/pipeline.rs`：早期单轮 `ha_ec` 路径，无生产调用点
  （仅 mod.rs 导出），git 历史保留。
- 修正 `bp_graph.rs` `output_contig_graph_primary` 文档注释：`break_ug_contig`
  自 f685eec 起已移植并在该函数内调用，旧注释仍称"not ported"。
- P3 阶段表中 "overlap 精化" 行由 ⏳ 改注：`append_window_list/
  ovlp_chain_gen/worker_ovec` 属旧 `ha_overlap_cal` 路径（仅 `--sec-in` 与
  C++ 死代码经过），默认纠错与 `--dbg-ovec` 均不经此，非缺口。

## Trio Stage A 进行中（2026-11-07，门 A 未过，定位到 rescue_bubble_by_chain）

### 已落地并验证正确的部分

- 语料：`scripts/gen_trio_corpus.sh`（读名 cksum 奇偶拆两个伪单倍型），名单
  `reference/trio_l1.txt`（2752, FATHER）/`trio_l2.txt`（2874, MOTHER）；
  C++ 基线 `reference/reference_trio.*`（`-t1 -3 l1 -4 l2`）。
- 前端：`assembly/trio_bin.rs`（`ha_triobin`/`_list`/`_set_list`，Trio.cpp
  450/328/294；khashl 换 std HashMap，首现优先语义保持；gz 嗅探；日志三行
  逐字一致——flagged 2752/2874 与 C++ 相同）。调用点在 `assemble()` 重算
  分支 `write_requested_outputs` 之后（Assembly.cpp:2113 语义）。
- `hit.rs::drop_edges_by_trio`（Overlaps.cpp:1590），插在两次
  `normalize_ma_hit_t_single_side_advance` 之后——gdb 取
  `cxx_predrop/cxx_postdrop2` 与 Rust d00a/d00b 掩码比对**完全一致**
  （del=127414 相同）。
- `gfa_clean.rs::get_tip_trio_infor`（Overlaps.cpp:6890，本地
  TRIO_THRES/DOUBLE_CHECK_THRES 用 f64 对齐 Overlaps.h 宏）+
  `asg_arc_cut_inexact`/`asg_arc_cut_length` 的 `is_trio` 三处分支
  （gfa_ut.cpp:1007-1041/1326-1366），assert 拆除，`ul_clean_gfa` 传
  `is_trio`（轮内两次传 trio，轮后 3157/3163 按 C++ 传 0）。
- 分发：pipeline.rs 按 Overlaps.cpp 39605-39651 顺序——`-5` poly Err →
  trio+Hi-C benchmark Err → trio → `output_trio_graph_joint`
  （清 HA_F_PARTITION、`(max_short_tip*2,0.15,3,…)`、打印收尾
  Inconsistency 行）→ Hi-C Err。
- **清洗全程逐轮验证一致**（gdb 断点 + write_debug_graph，掩码后）：
  round 0-3（cxx_r0-r3 vs d04-d07）、large_indel（d08）、
  cut_length_final（d09）、semi_circ（p10 vs d10）全部 **0 差异**。

### 剩余缺口（门 A 的唯一阻塞）

`post_rescue` 末尾的 `rescue_bubble_by_chain`（gfa_ut.cpp:3024 →
Overlaps.cpp:36832）。先前会话在默认语料上探测为 no-op 而跳过；trio
语料下 C++ 2387→2396 条活弧（p10→p14），Rust 停在 2387→2390。需移植：

- `bubble_type`（hic.h:37-53）、`get_bubbles`（hic.cpp:2217）、
  `dfs_bubble_broken`（hic.cpp:2232）、`identify_bubbles`
  （hic.cpp:2354-2607，Hi-C Phase 3 的机械提前入场）
- `reset_bub`（Overlaps.cpp:36534）
- `rescue_bubbles_by_contained_reads`（35620）、
  `rescue_bubbles_by_missing_ovlp`（36092）、
  `rescue_bubbles_by_missing_ovlp_backward`（36247）
- `rescue_bubble_by_chain` 本体（36832：ma_ug_gen_primary → 副本上调
  `adjust_utg_by_primary`（第三参 TRIO_THRES）→ reset_bub ×3 → 三个
  rescue；trio_recover 块在 0.25 已被注释，`no_trio_recover` 实参未用）
- Rust 已有可复用：`copy_read_graph`/`copy_untig_graph`/
  `ma_ug_gen_primary`/`adjust_utg_by_primary`（bp_graph.rs/unitig.rs）、
  `asg_bub_pop1_primary_trio`（bubble_identify.rs）、
  `asg_bub_pop1_primary_trio_switch_check`（hamming.rs:957）、
  `TransChain`（hap_cov.rs）。

### 侦察方法论备忘（下次直接用）

1. gdb `finish` 不可靠（cxx_postdrop 首次为假象）；应在**调用后的行**下
   断点（如 Overlaps.cpp:39375）再 dump。
2. `write_debug_ma_hit_ts`/`write_asg_t` dump 的掩码比对：ma_hit 的
   `cc/bl`（cigar 缓冲偏移，重算运行留陈旧值、装载运行为零、输出不可见、
   缓存不序列化）与 asg 弧的 `ou`（=min(bl,0x3fff)，Overlaps.cpp:3012，
   仅 UL 模式消费，Rust 有意置零——graph_build.rs:146）必须掩码。
3. `--dbg-gfa` 运行在 build_string_graph_without_clean 头部走
   `load_debug_graph` 分支且 HA_F_VERBOSE_GFA 改变纠错轮的 r_out 行为，
   其 39436 dump（1514 弧）与普通运行（2397 弧）不可比；一律用 gdb 在
   普通运行上取 dump。
4. pipeline.rs 现有临时探针 d00a_after_normalize/d00b_after_drop_edges
   （待门 A 过后决定去留）。

## Trio Stage A 续（2026-11-08，气泡子系统移植完成，门 A 推进到 purge 深链）

- `src/assembly/bubble_chain.rs`（~4200 行 Rust）：簇 1-4 全部落地——
  BubbleType/Pdq/identify_bubbles/build_bub_graph/get·detect_bub_graph/
  update·resolve_bubble_chain/tangle 簇/三 rescue/reset_bub/
  rescue_bubble_by_chain 驱动；hamming.rs 的 switch_check 转 pub；
  rescue.rs 的 get_best_contained_read/update_unitig 复用。
  cargo test 294 全绿。
- 接入：pipeline post_rescue 的 d15 之后（C++ post_rescue 末行），
  **暂限 trio_binning_requested()** 进入——默认语料 C++ 证明 no-op
  （d15==d16），跳过等价；默认门禁已复验通过（17 identical + ec.bin 豁免）。
- 已修 bug：驱动内 beg_idx 误用 s_bub，应为 f_bub（C++ rescue_bubble_by_chain
  的 beg_idx = bub.f_bub，只处理链/tangle 新增气泡）。
- 门 A 现行阻塞：rescue 内 adjust_utg_by_primary（trio 标签下）触发
  Purge_Dups.cpp:2964-3179 的 hap_alignment_advance_worker 候选扫描
  （purge_level_primary=3 在 trio 模式照样生效）→ 其 DP 链化助手簇
  （get_candidate_hap_alignment 等，~2000 行 C++）未移植。这是此前
  purge 深链 stub 清单的既知项，详见 bubble_rescue_map.md 的 ⚠️ 节。
- 遗留：簇 3 shortest_path 有一个未定位的越界 panic（默认语料触发过一次，
  trio 门控后当前不可达；hap-alignment 移植后 trio 全链路会暴露）。

## Trio Stage A 续二（2026-11-08，purge 深链打通 + 气泡子系统五处修复，门 A 差异收敛到 3 个读段）

### adjust_utg_by_primary 的 skip_purge 语义（根因修复之一）

C++ `adjust_utg_by_primary`（Overlaps.cpp:32425）在 `i_cov != NULL && collect_p_trans == 0`
时 `goto skip_purge`——跳过 purge、post-join rescue、tip 切割、读段重标记与
coverage-range 默认值，直接进 `recover_utg_by_coverage`。此前 Rust 移植漏掉该短路，
导致 trio 流程里 rescue 内的 adjust 白跑了整套 purge（产出 6 条 hap overlap，
后续 normalize 桩 Err）。`i_cov` 的有无在 Rust 中以新增 `has_cov_out: bool`
形参表达：rescue/ug_ext_gfa 调用点（`&cov, 0, 0` → skip）、`output_bp_graph`
（17642，`&cov, 1, 1` → 不跳）、`output_contig_graph_primary`（32654，`NULL`
→ 即使 collect_p_trans=0 也不跳）。**修改后 trio 运行的 purge 阈值打印次数
与 C++ 一致（1 次）**。

### hap-alignment 候选扫描簇（新模块 src/assembly/hap_alignment.rs）

`Purge_Dups.cpp` 611-3216 的活代码全量移植：prefilter/get_xy_pos/
deduplicate_edge/比较器/get_hap_overlapLen/classify_hap_overlap/
get_pair_hap_similarity(_by_base)/check_hap_match/determin_..._single_side/
bi_direction_.../get_hap_alignment_boundary/get_idx_by_base/get_base_boundary_
chain/advance(_back)/get_chain_score/determine_hap_overlap_type_advance(_back)/
calculate_pair_hap_similarity_advance/hap_chaining/get_candidate_hap_alignment/
filter_secondary_chain/init_pip_hh + `hap_alignment_scan_worker`（投票半部 +
候选扫描，替代 purge_dups.rs 的 Err 桩）。浮点混合精度（f32/f64 各比较点）、
无符号回绕、`push_hap_can` 不写 index_beg/score 的未初始化行为（Rust 取 0 并
注释）均在模块头记录。**接线后 trio 语料的 all_ovlp 现在能真实产出**（6 条），
但 skip_purge 修复后 C++ 在该语料不触发扫描——扫描代码为后续 purge 级语料
与 Hi-C 复用保留，当前默认/trio 路径不可达。

### 气泡子系统五处移植 bug（gdb 逐检查点对拍 C++ 后修复）

gdb 检查点（每次 `reset_bub` 后）：HIT1 14/74/13、HIT2 16/96/15、rescue 入口
18/106/17（num.n/list.n/b_g 节点数）——修复后 Rust 三处全部逐项一致，num 数组
逐字节相同，b_bub=2/b_end=2/tangle=0 一致：

1. **Pdq::id_mode 移位方向**：C `((uint64_t)-1)>>（64-i）`（低 i 位），Rust 误写
   `>>i`（2^56-1），unpack 出 `(dis<<8)|v` 的垃圾顶点（文档记载的
   shortest_path panic，本次终获定位）。
2. **get_bubbles 哨兵**：C 只在请求 body 时读 `num[id+1]`；update_bubble_graph
   新推气泡无哨兵，Rust 无条件读越界。Rust 改为缺哨兵时回退 list.len()
   （参考在该处读未初始化内存且 body 从不被消费）。
3. **扫描接受条件反相**：C `i == b.b.n` 为"校验循环完整跑完=接受"，Rust 写成
   反相——所有候选都被标 bit31 拒绝，气泡排序退化、b_s_idx/update 全盘失配。
4. **m_het/p_het 定义颠倒**：hic.h 实际 `P_het=num.n、M_het=num.n+1`
   （bubble_rescue_map 旧记录写反），Rust 跟随错误文档导致 index 分类错位、
   update_bubble_graph 的 hom 过滤吞掉 broken 集合。
5. **update_bubble_chain 三处**：end 分支弧计数须用链条图（bub_ug.graph）的
   `(i<<1)+1`/`(i<<1)`（Rust 误用 bsg）；middle 的 p_t 须为 bsg 弧槽位
   （C 持指针原地置 del，影响链条再生成；Rust 改传槽位下标）；`arc_first`
   无 count 守卫（按 offset 取槽位）。

另有 `update_unitig` revert 打包修复：C 的 `mode<<31` 后整体 `<<32`，mode 落在
bit63；Rust 直接 `mode<<31` 侵入 uId 字段（panic 0x80000012）。

### 门 A 现状

trio 运行完整跑通（exit 0），默认 9-profile 门禁无回归，294 单测全绿。
12 个 GFA/bed 产物仍有差异；d12 掩码对拍（write_debug_graph 转储 + ruIndex
解析）显示唯一分歧为 **3 个读段（1361/1455/1460）的容器指派**：C++ 未恢复
（0xffffffff），Rust 保留了陈旧的 contained 映射（0x549 等）——即
rescue_bubbles_by_contained_reads 的 chain_edges 收集在 Rust 侧漏了这 3 条。
下一层对拍点：rescue 内 chain 链收集（get_edge2existing_node_advance_by_
broken_bub / bub_complex 一带）对这 3 个读段的决策。

调试钩子（仅定位用，默认关闭）：`RSHIFIASM_PURGE_DEBUG`（打印 purge 后
all_ovlp 计数）；C++ 侧对拍脚本 /tmp/gdb_bub*.gdb（gdb 检查点转储）。

## Trio Stage A 完成（2026-11-08，门 A 通过：12 个产物逐字节一致）

### 收尾修复（bubble_chain.rs）

1. **`finish_bubble_rescue` 缺索引失效**：C++ 在 `asg_cleanup`/`asg_symm` 前执行
   `free(r_g->idx); r_g->idx = 0; r_g->is_srt = 0;`。Rust 的 `push_arc` 只追加，
   索引还指向旧区间——新增弧对 `arcs_of` 不可见，`symmetrize` 又把它当作
   非对称弧删掉。修法：push 之后、cleanup 之前 `index.clear(); is_sorted = false;`。
2. **backward rescue 的 `beg_idx` 传错**：C++ 传 `bub.f_bub`（full bubble 数），
   Rust 传了 `simple_bubble_count`，导致遍历的 m_bub 落在 [12,15) 而不是
   [15,18)（本语料 f_bub=15、s_bub=12）。
3. **backward rescue 的最终 recover_edges 回填同样缺索引失效**：修复同 (1)；
   修后 d16 与 C++ b3 转储在 ruIndex/trio 标签/coverage_cut/source/reverse/
   4786 条弧/序列上全部一致。

### 缺失的收尾阶段：`flat_soma_v`（关键）

`ul_clean_gfa` 末尾在 `output_unitig_graph` 之后还调用
`flat_soma_v(sg, src, rI)` → `flat_bubbles_advance`（Overlaps.cpp:38992）。
它在**读图**上反复弹出杂合气泡（`asg_bub_pop1_primary_trio` 先试探后真弹），
最后用 `clean_sg_by_utg`（Overlaps.cpp:38733）把读图重写成展开后的 unitig 图。
本语料 r_utg 有 95 条 unitig，展开后 p_utg 只剩 20 条——这正是此前
"Rust p_utg == r_utg（95 条）而 C++ p_utg 20 条"的根因；不是 ALTER_LABLE
（trio 模式下 `hic_clean_adv` 根本不执行：调用点 gfa_ut.cpp:39452 被
`if(!ha_opt_triobin(&asm_opt))` 门控）。

新增移植（bubble_identify.rs）：

- `flat_soma_v`：`dip_thre_max` 的三种来源（`somatic_cov >= 0` /
  `hom_global_coverage_set` / `hom_global_coverage / HOM_PEAK_RATE`），
  再 `*1.15 / polyploidy` 并截断为 u64（与 C++ 两次 double→u64 截断一致）。
- `flat_bubbles_advance`：`tLen = get_bub_pop_max_dist_advance`；每轮按
  `bs_flag` 状态（0 未试→1、1 跳过、2 已弹可重试）扫描；对每个成泡先用
  `is_pop=0` 试探并统计 r_flag/C_bases，`C_bases/R_bases <= het_thres` 时
  以 `is_pop=1` 真弹；轮次由 `n_pop` 驱动。
- `get_primary_path_len`：沿 `BubbleInfo.parent` 从 sink 回溯到（不含）v0。
- `clean_sg_by_utg`：先全删读图弧，再恢复 primary unitig 内部弧与端到端弧；
  ALTER_LABLE unitig 的读段整条删除；最后 `cleanup`。

接入点：`pipeline.rs` 的 `ul_clean_gfa` 等价段，`output_unitig_graph` 之后、
`clean_graph` 的 `output_contig_graph_primary_pre` 之前（新增 dump
`d19_after_flat_soma`）。

### 门禁结果

- **门 A（trio）**：`rshifiasm -t1 -3 reference/trio_l1.txt -4
  reference/trio_l2.txt -o /tmp/trio_wip/out chr11-2M.fa.gz`，
  `scripts/compare_outputs.sh reference/reference_trio.dip
  /tmp/trio_wip/out.dip` → 12 个产物（r_utg/p_utg 的 gfa/noseq/lowQ.bed +
  hap1/hap2 的 p_ctg gfa/noseq/lowQ.bed）**全部逐字节一致**。
- **9-profile 默认回归**：default / write_paf_ec / primary / bcov / telo /
  dbg / e / i / load 共 102 个 gfa/bed 的 SHA-256 与改动前完全相同；
  default 与 C++ `reference.bp.*`、write_paf_ec 与 C++ 同名产物（除已记录的
  `.ec.bin` padding 豁免）逐字节一致。`flat_soma_v` 在默认语料是 no-op
  （`hom_global_coverage` 路径下阈值不触发弹泡）。
- `cargo test`：326 项全绿。

调试钩子已全部移除（`RSHIFIASM_CHAIN_DEBUG`/`RSHIFIASM_BSG_DEBUG`、
BW-*/IR_HET/BUB/reset2/reset3 打印、b1/b2/ug_before_backward 临时转储）。

### 门 A 之后的 trio 剩余分支勘察（2026-11-08）

- **`--kpt-rate`（335）在 0.25.0-r726 的活路径上是死选项**：`asm_opt.kpt_rate`
  只在 `output_trio_graph`（Overlaps.cpp:21329-21363）里被读取，而该函数在
  `clean_graph` 的 trio 分发处已被注释（Overlaps.cpp:39636 起走
  `output_trio_graph_joint`）。Rust 现在不拒绝也不使用该参数，与 C++ 等价；
  不需要再实现。
- **`--trio-dual`（`output_bp_trio_graph`，Overlaps.cpp:17680）** 的缺口：
  `hic_clean`（已移植）→ `adjust_utg_by_primary`（已移植，需
  `collect_trans=1`）→ `trio_phasing_refine`（hic.cpp:17741）→
  `output_trio_graph_joint`（已移植）。`trio_phasing_refine` 还缺：
  `clean_u_trans_t_idx_filter_adv`（Overlaps.cpp:15503）、`asg_seq_set`、
  `ug_occ_w`/`ug_occ_hap_w`、`init_ps_t`/`destory_ps_t`（小件），以及
  `mc_solve`（rcut.cpp:3710，含 `mc_opt_init`/`init_mc_g_t`/`update_mc_edges`/
  `mb_solve_core`/`mc_solve_core_adv`/`p_nodes` 的完整 maxcut 簇，约 2k 行）。
  即 `--trio-dual` 与 Hi-C Phase 4 共用同一个 maxcut 依赖。
- **`-5` 多倍体（`output_poly_trio`，Overlaps.cpp:14682）** 额外需要
  `output_trio_unitig_graph` 的多倍体分支与 `mc_solve_general`（同一
  maxcut 模块），且需要多倍体测试语料（`-5` 三亲本名单），后置。
- **`--dual-scaf`（`self_scaffolding`）** 走
  `gen_self_scaf`/`gen_contig_path`/`double_scaffold`/`output_hap_sc_graph`
  与 `scaf_res_t` 簇（Overlaps.cpp:21585-23240，约 1.5k 行），独立于 maxcut，
  后置。

因此 trio 剩余分支的下一步是 maxcut（rcut.cpp）整体移植，它同时解锁
`--trio-dual`、`-5` 与 Hi-C Phase 4；`--dual-scaf` 的 scaffolding 簇次之；
`--kpt-rate` 已可视为完成。

## `--ex-list` / `--ex-iter`（2026-11-08，`extract.cpp` 全量移植）

语义（`ha_extract_print`，extract.cpp:82）：

1. 名单来自逗号串或 `@file`（gzip 透明，取每行首个空白前的 token）。
2. 命名的读段标为第 `-1` 层；随后 `--ex-iter` 轮逐层扩展：对两侧 overlap
   表逐读逐命中，恰有一端 `< round` 时把另一端登记为 `round` 层。
3. **轮次循环只扩展、不打印**；循环结束后一次性转储所有「至少一端已被发现」
   的命中，字段为
   `qname qlen qs qe ± tname tlen ts te ml bl !l`（`!l` 对 forward 表为 1、
   reverse 表为 0）。
4. 仅在重启缓存加载成功时触发（`Assembly.cpp:2065`），打印后 `exit(0)`；
   `-i`（不读缓存）时与 C++ 一样被忽略。

门禁（`reference/reference` 的 C++ 缓存 vs `reference/candidate_default` 的
Rust 缓存，stdout 逐字节）：

| 名单 | 逗号 1 名 | 逗号 3 名 | `@file` | `@file.gz` |
|---|---|---|---|---|
| `--ex-iter 0/1/2/3` | 114/5358/8748/12256 行，全部 identical | 364/5878/8748/12256 行，全部 identical | identical | identical |

`cargo test` 326 全绿；默认 9-profile（102 个 gfa/bed 的 SHA-256）与 trio 门 A
均无回归。实现：`src/assembly/extract.rs`（新），接线在
`assembly::assemble` 的缓存加载分支、`hap_recalculate_peaks` 之前。

## `--prt-raw`（2026-11-08，`gen_init_sg` 的原始 unitig 转储）

`asm_opt.prt_dbg_gfa` 在 `gen_init_sg` 的 `ma_sg_gen` 之后、`asg_arc_del_trans`
之前触发一次 `prt_dbg_gfa(sg, "raw", ...)`：
`get_outfile_name()` 前缀 + `.raw`，再经 `print_debug_gfa(ug=NULL, is_seq=1)`
执行 `ma_ug_gen` → `ma_ug_seq(..., is_polish=0)` → `ma_ug_print`，写出
`{o}.{bp|dip|hic}.raw.r_utg.gfa`（带序列、A 行与 HG 标签）。

门禁：default `-t1 --prt-raw` 与 trio `-3/-4 --prt-raw` 的
`*.bp.raw.r_utg.gfa`、`*.dip.raw.r_utg.gfa` 均与 C++ 逐字节一致
（29,881,620 / 42,477,043 字节）；关闭该选项时行为不变（默认 9-profile 与
trio 门 A 复验全绿，`cargo test` 326 项）。实现：`graph_build::write_raw_debug_gfa`。

## `--dbg-het-cnt`（2026-11-08，`het_cnt.log` 兼容）

`asm_opt.is_dbg_het_cnt` 在最后一轮纠错分配 `het_cnt`，`cal_ec_r` 结束后由
`print_het_cnt_log` 写 `{o}.het_cnt.log`（`>{name}\t{count}` 每读一行）。
**0.25.0-r726 的活路径里数组从不会被填充**——只有已死的 `worker_ovec` /
`worker_ovec_cal0` 会写 `het_cnt`，`cal_ec_r` 的 `worker_hap_ec` 不写。因此
C++ 的日志 5626 行全为 0；移植按同一行为写出全零，并在 `ec_driver` 注释中
记录原因（不移植死 worker）。

门禁：`-t1 -i --dbg-het-cnt` 的 C++ 与 Rust 运行
（`/tmp/hetcheck/cxx|rs`）`.het_cnt.log` 逐字节一致（5626 行全 0），
且同一次运行的全部 GFA/bed 产物也逐字节一致；`cargo test` 326 项；
默认/trio 门禁无回归。

## maxcut 的 hap-overlap 路径、purge_dups 深链与 `--trio-dual`（2026-11-08）

本节收尾 `--trio-dual`，同时补齐 `--trio-dual` 所依赖的 purge_dups 深链。此前
`purge_dups` 只移植到「归一化之前」；`adjust_utg_by_primary(collect_p_trans=1,
collect_p_trans_f=0)` 需要的后半段（purge maxcut、含 hap-overlap 的 maxcut 求解、
contained haplotig 删除、purge 图构弧、泡弹出、弱弧剪切、unitig 合并与重连）
全部缺失。

### maxcut 的 hap-overlap 分支（`src/assembly/maxcut.rs`）

C++ 的 `mc_solve(NULL, NULL, &in, ug, sg, 0.8, NULL, 0, s->s, 1, bub, ref,
0, 0)`（hic.cpp:17741 `trio_phasing_refine`）走的是「ovlp 建边」而不是
`mc_match_t` 建边，且 `bub`/`ref` 均为非空，因此 `mb_solve_core` 与
`mc_solve_core_adv` 的带泡、带 ref 分支都会真实执行。移植内容：

- `update_mc_edges_ovlp`（maxcut.rs:469）：按 `hap_cov` 的 `HapOverlap`
  列表给 `mc_g_t` 建边。过滤条件与 C++ 一致：两端 `u_id`、`hs` 相位、
  `mc_interval_filter` 的 het/hom 长度比例 `f_rate=0.8`；反向相位的边标
  `FLIP`。
- `filter_ovlp_by_mc` / `clean_ovlp_by_mc`（maxcut.rs:533,559）：求解结束后
  按 spin 一致性镜像保留、再 `retain(status == FLIP)`，只留下应删除的弱边。
- `set_p_flag` / `p_nodes`（maxcut.rs:591,610）：仅当
  `(flags & FLAG_PARTITION) && t_ch` 时标记读段：`SET_TRIO(8)`，再按 spin
  置 `FATHER(1)`/`MOTHER(2)`；`spin == 0 && P_HET` 也归 `FATHER`。这是
  trio 门禁里 `trio_flags` 被改写的唯一来源。
- `mc_solve_purge`（maxcut.rs:3042）：purge 专用入口，序列数取
  `unitig_graph.graph.sequences.len()`，`is_sys=1`、`renew_s=1`，调用
  `mc_solve_core_adv` 时不传 `unitig_graph`/`bub`（与 C++ 传 `NULL` 的
  `ug`/`bub` 对齐，`ref` 仍非空以触发 `mb_solve_core`）。
- 修掉一处移植期引入的截断：`mc.cc`/`mb.cc` 的元素是
  `(group << 32) | vertex`，Rust 里所有下标必须 `as u32 as usize`
  （C++ 侧是 `uint32_t k`）；否则 `b.s[k]` 会在多组时越界。

### purge_dups 深链（`src/assembly/purge_dups.rs`）

按 C++ 顺序补齐，函数与 Purge_Dups.cpp 的对应关系：

| Rust | C++ |
|---|---|
| `normalize_hap_overlaps_advance_by_p_g_t`（:903） | 3190-3260 两端口过滤 + DELETE 压缩 + 顺序镜像 |
| `PurgeAlignContext` / `chain_origin_trans_uid_by_purge`（:960） | `chain_origin_trans_uid_by_distance` 的 purge 专用包装 |
| `collect_purge_trans_cov`（:1121） | 3354-3500 收集 RC_2 转换 |
| `remove_contained_haplotig`（:1208） | 4040-4460 删除被包含的 haplotig |
| `get_hap_arch`（:1243） | 4590-4680 由 `ma_hit_to_arc` 生成 `pg_h_lev` 弧，score 编入 `overlap`/`ul` |
| `clean_purge_graph`（:1599） | 4870-4990 弱弧剪切 |
| `fill_unitig` / `get_node_boundary_advance` / `purge_merge`（:1613/:1664/:1855） | 4990-5730 链式合并 |
| `link_unitigs`（:2046） | 合并后重连 |

关键移植细节：

- `new_rtg_edges` 在 `adjust_utg_by_primary` 中由 purge 追加连接弧，随后被
  `fill_unitig` 与 `renew_utg` 消费，因此参数从 `&[AsgArc]` 改为
  `&mut Vec<AsgArc>`（`bp_graph.rs`、`trio_dual.rs`、`bubble_chain.rs` 与
  `tests/bp_graph_golden.rs` 同步）。
- score 编码沿用 C++：`encode_score` 把符号放最高位，`decode_score` 还原；
  `purge_g_arc_del_short_diploid_by_score` 用 `f32` 比较
  `(ovlp as f32) > (ovlp_max as f32) * drop_ratio`，与 C++ 的浮点语义一致。
- `purge_g->seq[v > 1]` 这类 C++ 边界行为按原样复刻（v 为 0/1 时不取）。
- 唯一的剩余回退是 `if_ploid_sample`（既无 coverage 峰也无 k-mer 峰时按
  样本定倍性），在 0.25.0-r726 活路径中不会被触发。

### `worker_for_trans_clean_re` 的 `get_u_trans_spec` 修复

`clean_u_trans_t_idx_filter_adv` 第二遍按「query 侧配对」选代表。原实现直接
拿 `target` 的完整 transition 区间，C++ 实际走
`get_u_trans_spec(ta, a[l].tn, a[l].qn, &r_a, &r_n)`（Overlaps.cpp:17000），
只取该 target 中**恰好指向本 query** 的子区间。差异表现为：Rust 的 purge
maxcut 边数 62，C++ 64；gdb 抓到的 C++ `in` 列表里多出
`(qn=7, qs=1215651, qe=1469660, tn=10, ts=0, te=254028, nw=127, f=2, rev=1)`
及镜像对。`trans_clean.rs:1084` 已按 `get_u_trans_spec` 的子区间扫描并
复现该行为。

### 门禁结果（全部 `-t1`，逐字节）

| 门禁 | 结果 |
|---|---|
| `--trio-dual`（12 个产物：r_utg/p_utg vs `reference_trio.dip.*`；hap1/hap2 vs `reference_trio_dual.dip.*`） | 12/12 identical |
| triobin（6 个产物 vs `reference_trio.dip.*`） | 6/6 identical |
| default 配置（17 个 gfa/bed vs `reference.bp.*`） | 17/17 identical；仅 `.ec.bin` 因 padding 不同（既有豁免） |
| `--primary`（vs `reference_primary.*`） | 全部 gfa/bed identical；仅 `.ec.bin` 例外 |
| write_paf_ec（vs `reference/write_paf_ec.*`） | gfa/bed/ec.fa/ovlp.paf/ovlp bin 全部 identical；仅 `.ec.bin` 例外 |
| `cargo test --release` | 全绿（含 3 项 `bp_graph_golden`） |

至此 `--trio-dual` 在功能与门禁上完成；trio 仅剩 `-5`（`output_poly_trio`
+ `mc_solve_general`）与 `--dual-scaf`（`gen_self_scaf` 簇）。

## `-5` 多倍体 trio 输出（2026-11-08，`output_poly_trio`）

`-5 <file>` 是 `clean_graph` 输出分派的第一支（Overlaps.cpp:39613）：
`flag` 去掉 `HA_F_PARTITION` 后调用
`output_poly_trio(..., hapN = asm_opt.polyploidy)`。它从外部名单构造每条读段的
单倍型位图，然后对每个单倍型重写读段标签并输出一整套 `hapN.p_ctg.*`。

移植内容：

- `ha_polybin_list`（`src/assembly/trio_bin.rs`）：先按第一次出现建立读名到
  下标表（与 `ha_triobin_list` 相同），再逐行解析名单。`phrase_hstatus`
  （Trio.cpp:373）的格式语义完整复刻：第一字段截到首个空白，`hid` 取自
  第一字段最后一个 `_` 之后（`atoi` 语义），第二字段是读名；缺空白、缺
  `_` 或 `hid == (uint32_t)-1` 的行打印 `ERROR: wrong hap status` 后跳过。
  文件按 gzip magic 透明解压（`gzopen` 行为）；打开失败打印
  `ERROR: failed to open file '<path>'`。位图用
  `1u32.wrapping_shl(hid-1 & 31)`，对越界 `hid` 复现 x86 的移位数掩码。
- `update_poly_trio`（Overlaps.cpp:14646）：`DROP(5)` 保持不变；其余先置
  `AMBIGU(0)`，位图为 0 的读段留 `AMBIGU`，否则按当前掩码位选
  `FATHER(1)`/`MOTHER(2)`。即不在任何 bin 里的读段是 `AMBIGU`，而在其他
  bin、不在当前 bin 的读段是 `MOTHER`。
- `output_trio_unitig_graph` 的 `is_bench=0` 分支（`src/assembly/bp_graph.rs`）：
  `ma_ug_gen_phase` → `adjust_utg_by_trio` → 可选的两个 `break_ug_contig` →
  `ma_ug_seq` → `ma_ug_print`/`ma_ug_print_simple`/`ma_ug_print_bed`，文件名
  前缀为 `hap1..hapN`，GFA 标签固定 `h1tg`（C++ 传 `FATHER`）。原实现只在
  `is_bench` 路径下返回图，现在两条路径都可用；joint 路径仍传 `None` 前缀。
- `output_poly_trio`（bp_graph.rs）按 C++ 顺序逐单倍型调用
  `update_poly_trio(1<<i, ...)` 与上面的输出，`bubble_labels` 在循环间共享；
  全部分支结束后统一打印 BED 不一致阈值。

门禁（`-t1`，C++ 16.5 GB / Rust 7.5 GB 峰值，故参考运行须单独执行）：

| 名单 | 产物 | 结果 |
|---|---|---|
| 2 bin，`--n-hap 2`（默认） | hap1/hap2 的 p_ctg gfa/noseq/bed + p_utg/r_utg 等 15 项 | 全部 identical；仅 `.ec.bin` 例外（既有 padding 豁免） |
| 3 bin，`--n-hap 3` | hap1/hap2/hap3 同上共 18 项 | 全部 identical；仅 `.ec.bin` 例外 |

名单由读名哈希确定性生成（`x_<1+hash%hapN>\t<read_name>`），C++ 与 Rust 使用
同一文件。trio 剩余分支只剩 `--dual-scaf`。

## `--dual-scaf` 的 UL 基础（2026-11-08，进行中）

`--dual-scaf` 是 trio 的最后一个分支（`output_trio_graph_joint` 内，
Overlaps.cpp:23282 `gen_self_scaf`）。它与 `--ul` 共用 `inter.cpp` 的
ultra-long 索引/链化引擎，规模远大于此前各分支。本轮先落地可独立验证的
基础层，新增 `src/assembly/ul.rs`：

| Rust | C++ |
|---|---|
| `UlOv` / `UcBlock` / `UlVec` / `ScafRes` | `ul_ov_t` / `uc_block_t` / `ul_vec_t` / `scaf_res_t` |
| `UtgRid` / `UtgRidDt` / `gen_r_ug_idx` | `utg_rid_t` / `utg_rid_dt` / `gen_r_ug_idx`（21278） |
| `get_r_ug_region` | 13268 |
| `rov2uov_ctg` | 16223 |
| `get_unique_rctg_aln` | 16269 |
| `ctg_rg2ug_gen` | 16323 |
| `linear_ctg_chain_dp_adv` | 16391 |
| `gen_linear_chains_ctg` | 16497 |

要点：

- `gen_r_ug_idx` 的 CSR 填充复刻参考的「末槽 `off` 兼作游标」写法；前缀和
  从 `k = 0` 开始（`idx[0]` 也必须置 0），`keep_offset` 时 `off = k`。
- `rov2uov_ctg` / `get_unique_rctg_aln` 的 `rev` 取反、`adjust_rev_tse`
  的坐标翻转、`over_bd` 的圆环裁剪均按原文；`IS_SCAF_READ` 在纯 HiFi 输入
  下恒为假，故省略。
- 线性链 DP 的 `f`/`p`（`i32`/`i64`）与 `t` 跳过计数、`max_ii` 距离窗口、
  链压缩后 `qn = k + chain_offset`、`sv` 的 `adjust_rev_tse` 全部逐行对应；
  排序使用稳定排序并按 `qe` → 组内 `qs` 二级键，等价于参考的两次 radix
  sort（同为稳定排序，平序置换一致）。
- `gen_linear_chains_ctg` 的 `tn >>= 1` 发生在 DP 之前，`is_available_ctg_aln`
  的 16 bp 容差与 `sec >= n*0.333333` 门槛保持不变。

单元测试覆盖：合成 unitig 上的 `r_ug` 索引偏移、连续读段合并、
`sec` 累加与全长对齐保留（`cargo test --release --lib assembly::ul`）。

### 已打通（2026-11-08 续）

- `gl_chain_linear`、`find_mg_lchain_max`、`set_ul_ov_t_by_mg_lchain_t`、
  `infer_rovlp`、`get_ecov_adv`（mode 0）、`select_max_ctg_chain`、
  `gen_gchain_track`、`hc_shortest_k`、`dump_linear_chain`、`push_ctg_res`、
  `gen_contig_path` 全部落地于 `ul.rs`；`gen_contig_path` 对 cp0/cp1
  与参考在 gdb 断点（`double_scaffold` 入口）处逐块比对，6+2 条 contig 的
  链头与数据块（hid/qs/qe/ts/te/aidx/rev）完全一致。
- 关键修正：`select_max_ctg_chain` 的链头 `g_item` 必须在 `gen_gchain_track`
  之前 `kv_pushp`（此前漏推，导致 `gchains[len-res_n-1]` 下溢 panic）；
  `push_ctg_res` 的数据块写到 `bb[header.ts + z]`，不是 `bb[k + z]`。
- `ul_vec_t.bb.n` 语义：参考只把**链头数**记入 `n`，数据块位于容量区
  （`bb[ts..te]`）。Rust 新增 `UlVec::n` 保存链头数，`bb` 保留全部数据。
- `gen_contig_path` 的 `diff_ec_ul` 是 `asm_opt.ul_error_rate`（默认 **0.2**，
  不是图链的 `N_GCHAIN_RATE=0.04`）；`bw=G_CHAIN_BW=16`。

### `double_scaffold` 前半（新增 `src/assembly/dual_scaf.rs`）

| Rust | C++ |
|---|---|
| `gen_clean_ug` / `gen_shallow_clean_ug` | Overlaps.cpp 21623 / 21685 |
| `gen_shallow_ref_trans` / `gen_deep_ref_trans` | 22850 / 22868 |
| `merge_utgs` / `merge_scaf_res` / `get_scp_len` / `gen_scaf_res_via` | 22910 / 22970 / 22992 / 23020 |
| `gen_double_scaffold_gfa` / `gen_scaf_path_lst` / `get_arcs` | 22689 / 22578 / gfa_ut.cpp 463 |
| `gen_scaffold_bases` / `ma_scg_print` | 23162 / 11025 |
| `gen_self_scaf` | 23209 |

要点：

- `gen_clean_ug` 先把**所有** seq/arc 置 `del=1`，再沿 cp0/cp1 的链把节点
  复活、把链上相邻弧 `asg_arc_del(..., 0)` **复活**；`gen_shallow_clean_ug`
  只对 `bubble.index < f_bub`（无需打断的小气泡）相关弧做删除。
- `gen_deep_ref_trans` 的 `ccref = copy_untig_graph(cref)` 必须在
  `adjust_utg_by_primary` **之前**取快照；`ma_ug_seq(ccref, sg, ..., NULL, 0)`
  是 `is_polish=0`、`edge=NULL`（仅物化序列）。
- 本语料（chr11-2M + 2752/2874 条 trio 读段）C++ 的
  `work_ctg_path_trans_self` 打印 `# try:8, # done:0`；`# done` 是
  `direct_ctg_trans_chain` 的恒定返回值（push 后即 `return 0`），
  不代表空表——实际产出 10 条过渡，见下文“过渡链接”一节。
- 门禁结果：`out.dip.hap1.p_ctg.gfa`（3,019,998 B）与
  `out.dip.hap2.p_ctg.gfa`（2,310,344 B）与 C++ **逐字节一致**，其余
  p_utg/r_utg 产物沿用既有门禁。

### `--dual-scaf` 过渡链接（2026-11-08 续二，已完成）

早期观察到的 `# try:8, # done:0` 曾被误读为“无过渡”：`direct_ctg_trans_chain`
在 push 完过渡后无条件 `return 0`，其后的代码被注释掉了，所以 `# done:0`
不代表空表。用 gdb 在 `inter.cpp:17362/17367/17372/17374` 打印 `ll->lo.n`
后确认参考实际产出 **10 条过渡**（各查询条数 1,1,3,1,1,1,1,1），Rust
`ctg_trans.rs` 的条数与逐条内容与之相同。

本轮补齐的部件：

| Rust | C++ |
|---|---|
| `order_contig_trans` / `rm_dup_aln` | inter.cpp 21955 / 17483 |
| `scaf_tangle::is_noisy_ov` / `set_fsk` / `scaf_hit2arc` | Overlaps.cpp 21888 / 22648 / 21876 |
| `scaf_tangle::get_ref_nid_flank` / `dedup_tang_flank` | 22040 / 22314 |
| `scaf_tangle::cal_self_close_tang_flank{0,1}` / `trace_back_tangle_dis` | 22348 / 22441 / 22497 |
| `scaf_tangle::is_small_tangle` | 22547 |
| `gen_double_scaffold_gfa` 的每查询过渡循环 + 分支裁剪 + 路径枚举 | 22689 |

要点（易错处）：

- `get_ref_nid_flank` 里 `va/vd/wa/wd` 四个指针在 `dedup_tang_flank`
  **之前**用原始 `vn0/wn0` 算好；dedup 只搬移 `va`/`vd` 内容并缩小
  `vn`，`wa` 始终位于 `2*vn0`。因此切分必须是
  `vd=flank[vn0..]`、`wa=flank[2*vn0..]`、`wd=flank[2*vn0+wn0..]`。
- `dedup_tang_flank` 的循环条件 `k<vn && k<wn` 每轮重新求值，
  不能用一次性的 `0..min(vn,wn)`。
- `cal_self_close_tang_flank0` 的**反向半**（`v^=1` 之后）遇到
  `fs[v^1]==ef` 只清标志、不 push；只有正向半会 push `v^1`。
- `get_ref_nid_flank` 的第二次 `rev` 反转作用在 `res->n - rn` 指针上，
  即最后 `rn` 个词（距离块），不是节点块；节点块在 push 距离前已反转一次。
- `rm_dup_aln` 的 contained 过滤按 `(nw<<32)|k` 稳定排序，重复项中
  `nw` 大者优先；`order_contig_trans` 的 mm=1 会**复活** del 状态不一致
  的一对并标 `occ=-1`（不可靠），mm=2 再把不可靠区按 `qs` 合并成区间。
- 本语料 10 条过渡在链接阶段没有一条通过 `el` 判据（输出中 0 个 N gap，
  hap1/hap2 仍是 6/2 条 singleton scaffold），但代码路径已完整执行且与
  参考逐字节一致；tangle 的数值正确性需更大语料才能被输出覆盖。

### `--dual-scaf` 的 `trans_base_infer`（2026-11-09，已完成）

`trans_base_infer`（inter.cpp 22574）是最后一个未移植的 `--dual-scaf`
部件。本轮补齐并接线：

| Rust | C++ |
|---|---|
| `trans_base::trans_base_infer` | inter.cpp 22574（`init_ug_trans_t` + `gen_trans_base_count_comp`） |
| `trans_base::gen_trans_adaptive_aln` | inter.cpp 10297 |
| `trans_base::push_rc_3_transitions` | inter.cpp 22234-22301（tail 合并） |
| `trans_base::build_trans_base_index` | `gen_trans_base_count_comp` 的 count/sort/remap 半 |
| `qgen_input::minimizers_qgen_input` / `lchain_qgen_mcopy_input` | `ug_map_lchain` 的种子与链化半 |
| `ug_lalign::split_ug_lalign` / `ovlp_base_direct` / `ul_raw_aln` / `rechain_alignment_ul` | inter.cpp 9967 / 17367 / 24913 / 17145 |
| `target::TargetSource` | `uref`/`rref` 的统一读取 |
| `trans_clean::clean_u_trans_t_idx_adv` | Overlaps.cpp 15384（入口先重建过渡表） |

要点：

- 参数：`k=trans_mer_length`（31）、`w=trans_win`（31）、
  `bw_thres=diff_ec_ul=1-trans_base_rate`（0.07）、
  `bw_thres_sec=diff_ec_ul_double=1-trans_base_rate_sec`（0.5）、
  `sec_cutoff=0.85`、`mini_cut=512`、`chain_cut=3`、`high_occ=polyploidy+1`、
  `mcopy_rate=0.2`/`mcopy_khit_cut=3`（`enable_mcopy=0`，不生效）。
  注意 `init_ug_trans_t` 的调用顺序在 `ha_opt_update_cov` 之后，所以
  `max_n_chain` 取更新后的值。
- `worker_for_trans_ovlp_adv` 单线程顺序执行即 `kt_for(n_thread=1)`；
  `ll[tid].tk` 的合并退化为按查询升序遍历 track 的连续 `qn` 段。
- 两轮 `gen_trans_adaptive_aln`：首轮 `fi=filter`、`err_high=0.5`/
  `err_low=0.07`，返回 `test_het_aln`；失败时第二轮 `fi=NULL`、
  两个误差率均为 0.5，并追加 `remove_trans_ovlp_connect` 弧过滤。
- `RC_3` 的 `occ` 字段在参考里未赋值，随后 `clean_u_trans_t_idx_filter_adv`
  会先把它重置为 0，因此填 0 不影响可观测结果。
- **`rechain_aln` 在 UL 路径上是空操作**：`clear_overlap_region_alloc`
  会清空所有槽位的 `w_list`，`gen_win_chain` 看到的 `wn` 恒为 0 并返回 0。
  Rust 移植版此前会因 `window_count - 1` 下溢 panic，本轮在
  `generate_window_chain` 增加空列表短路，与参考的可观测行为一致。

门禁（`-t1`，语料 chr11-2M + `-3/-4` 名单）：

| 运行 | `hap1.p_ctg.gfa` | `hap2.p_ctg.gfa` |
|---|---|---|
| C++ `hifiasm --dual-scaf` | 3,019,998 B | 2,310,344 B |
| Rust 同命令 | 3,019,998 B | 2,310,344 B |

`scripts/compare_outputs.sh`：`dip.hap1/hap2` 的 p_ctg/p_utg/r_utg
gfa/noseq/lowQ.bed 共 9 项全部 identical，`.ovlp.source.bin`/
`.ovlp.reverse.bin` 亦 identical；仅 `.ec.bin` 以既有的 padding 豁免
（同尺寸 21,361,791 B，差异为散落的 padding 字节）。

数值层验证（`RSHIFIASM_PIPELINE_DUMP_DIR` 打开探针）：

- `trans_base_trans.txt`：58 条过渡（38 `RC_3` / 14 `RC_2` / 2 `RC_1` /
  4 `RC_0`），track 共 64 条。
- `dual_trans.txt`（`gen_contig_self` 之后）与旧记录中的 C++ 参考值逐条吻合：
  `nw=44144/161115/253461`、`qs=1262502`（此前 Rust 为
  `132565/161513/254026`、`qs=1262523`）。文档原先记录的“深过渡表数值差异”
  已关闭。

UL 基础与 `--dual-scaf` 已接入 `output_trio_graph_joint`
（`options.self_scaffolding != 0`）；默认/trio-dual/-5 路径不受影响，
`cargo test --release` 全绿（333 passed）。
