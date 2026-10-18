# 验证记录（2026-10-18）

当前简化版已移除参数校验、输入标签检查和比对后检查。下方相关检查结果
属于早先版本的开发记录，不再是当前流程中的任务。

## 已执行

- 真实 Y5Y-WT-HA 下机 BAM + 用户提供的 GRCh38.p14 FASTA：任务依赖构建通过。
- modified + 未比对输入：抽取 24 条真实 reads，创建临时测试参考，实际执行
  索引、pbmm2 比对、注释、FIRE、提取、QC 和 BAM 索引；全部通过。
- 测试参考的一半序列做了反向互补，用于覆盖反向比对。比对前后 24 条
  reads 的 modified_bases_forward 修饰位置及置信度逐项一致（排序后比较）。
- modified + 已比对输入：小样本实际执行通过。
- annotated + 已比对输入：复用核小体/MSP 注释，小样本实际执行通过。
- hifi 输入：从已有 CCS 测试 BAM 取 4 条含 fi/fp/ri/rp 的 reads，
  实际执行 predict-m6a 及下游分析通过。
- FiberHMM 可选分支：真实小样本执行通过，输出 BAM 和 BAI。
- FiberFold 可选分支：四类模拟 bigWig + 一个 2,097,152 bp 窗口，
  使用作者模型，通过 Snakemake 生成预测文件。PyTorch 多进程通信需要
  在当前沙箱外执行；CPU 推理成功，未验证 GPU。
- FIRE peaks 输入导出：实际生成带 sample/bam 表头的 manifest 和配置。
- 错误模式检查：无 kinetics 却选择 hifi、未对齐却声明已对齐、缺少区间
  注释却选择 annotated，均被拒绝。
- 另一个配置文件指定 samples 时，不会意外保留默认的开发样本。

## 测试材料和记录

临时 fixture 和各模式配置：`/tmp/fiberseq-smoke-ic7xif21/`。
临时测试参考来自 reads 本身，仅用于验证软件衔接，不能解释为 hg38 比对
准确性、有效覆盖度或生物学结果验证。测试 BAM 和模拟轨道不加入 Git。

执行日志位于 `/tmp/fiberseq-*-run.txt`，fixture 下还保留每个规则的日志。
临时目录可能被系统清理。本文件记录验证结论，不要求这些目录长期存在。

## 初版运行边界（后续更新见文末）

- 没有启动 44 GB 真实样本的全量比对和分析，也未构建完整人参考索引。
- 真实数据的全面 QC、有效深度、核小体周期性与结果合理性仍需验证。
- 完整 FIRE peaks/FDR 校准在官方流程中运行，当前只自动导出其输入。
- FiberFold 的 CpG 聚合轨道、CTCF motif/footprint/direction 轨道未自动生成。
  当前可选分支接收准备好的四类标准化 bigWig。
- 5hmC 标签保留在 BAM 中，当前不单独导出 5hmC 轨道。
- `ft qc` 输出统计表；尚未加入 QC 图表/HTML 报告。

简化版复测：移除预检查后，真实配置 dry-run 和 24 条 reads 的小样本主流程
（含 FiberHMM 可选分支）均通过。


CTCF 更新：motif、足迹比例和方向 bigWig 已自动生成，并接入FiberFold 预测步骤。
旧记录中“CTCF 轨道未自动生成”的限制已解除；FIRE/CpG 聚合和窗口生成
仍未自动化。具体方法与测试见 ctcf.md。


FiberFold 输入准备更新：FiberFold 输入准备模块已实现 FIRE/CpG 标准化轨道和预测窗口。
旧记录中这三项“未自动化”的限制已解除。未比对的小样本现已通过
索引、比对、核小体/MSP、FIRE、四类轨道、窗口和模型预测的 16 步链路。
仅为技术验证；实际三样本全量分析未执行。详见 fiberfold_inputs.md。


目录统一更新：规则、结果和日志编号现已一致，参考索引并入
01_align/reference，可选甲基化输出为 optional/methylation。
在新的临时结果目录中跑完 55 个小样本任务（包含 FiberHMM 和 FiberFold），
并实际验证 05_fire_peaks_inputs 的导出与日志。可选甲基化分支仅验证
dry-run 路径，未安装或执行 methbat。真实大样本全量分析仍未启动。

## FIRE peaks 接入：原生命令验证与兼容处理（2026-10-18）

第 5 步调整为 `05_fire_peaks.smk`，增加真正的群体水平 peak calling：
对第 4 步 FIRE BAM 生成同染色体随机背景，使用 `ft pileup --shuffle`
计算真实/随机分数，以区间长度加权的分数直方图累计
`shuffled_bp / real_bp` 生成 FDR 表，再通过
`ft call-peaks --fdr-table` 调用峰。FDR 限制在 0–1，不使用会丢失阈值的
表格简化；比例筛选叠加于 FDR，不能替代 FDR。

真实样本预期输出为 `05_fire_peaks/<sample>.peaks.bed.gz`、`.tbi` 和
`<sample>.fdr.tsv`，并保留随机分子、排除区间及覆盖度摘要。日志也使用
`05_fire_peaks/`。默认启用，独立目标为 `fire_peaks`。
`prepare_fire_peaks` 仍仅导出运行完整官方流程所需的样本表和配置。
旧记录中的“当前只自动导出其输入”描述对应早先版本。

本轮已经完成的原生命令实测：

- 在 `/tmp/fiberseq-native-peaks-BASIiJ/` 构造 740 条 reads 的合成 BAM，
  包含强 FIRE 区域、较弱信号及背景分子，并带 m6A 和至少 10 个 MSP。
  `build_fixture.py` 可重建材料，`compare.py` 和 `summary.json` 保存对比。
- 本机 fibertools-rs 0.13.1 接收完整/部分/空/原位置的 `--shuffled` 文件，
  甚至不存在的路径或省略该参数，均产生内容哈希完全相同的 FDR 表和峰。
  日志显示它自行生成随机位置，未使用外部 BED。因此不直接依赖该参数
  实现排除区域和用户指定随机背景。
- `ft call-peaks` 在这份合成数据上以 `--max-fdr 0.05` 产生 11 个峰和
  32 行 FDR 阈值；首峰 `chrTest:5000–5200`，FDR 为 0.002。
  这是命令调用证据，不代表对生物学真实峰的验证。
- `--fdr-table` 能加载已有表、跳过内部 FDR 表构建，并输出同样的 11 个峰。
- `.bed.gz` 是 BGZF，`pysam.tabix_index(..., preset="bed")` 成功。
  FDR 表头为 `threshold / FDR / shuffled_bp / real_bp`；峰表包含坐标、
  FDR、覆盖数、FIRE 覆盖数、分数和 `pass_coverage`。
- 有效 FDR 表下，把附加 FIRE 比例阈值设为不可能满足的 1.1，命令返回 0，
  输出仅含表头的 BGZF。反之，无法找到满足目标 FDR 的分数阈值时，
  原生命令返回非零；这两种情况不能混为一谈。
- 仅有强峰的极端合成数据还暴露出原生 FDR 表简化会丢失高阈值的问题；
  集成方案自行累计并保留阈值，避免依赖该简化步骤。

新增分支联调也已完成，记录位于 `/tmp/fiberseq-peak-integration/`：

- `snakemake.yaml` 以同一合成 BAM 模拟三个样本，复用预置的第 04 步 BAM。
  `fire_peaks` 目标实际完成 10 个任务：每样本背景、FDR、call peaks 各一个，
  加上汇总目标。日志为 `snakemake-run.txt`。
- 每个测试样本成功 shuffle 740 个分子，覆盖中位数 249，范围 170–328；
  流式统计按 0.01 分数分箱，生成 390 个 FDR 阈值（另有一行表头）。
- 每样本输出 10 个通过覆盖筛选且 FDR ≤ 0.05 的峰，并用 tabix 逐条读回。
  首峰为 `chrTest:5000–5200`。与原生命令的 11 峰记录不同，本次外部 FDR
  及覆盖筛选生效：有一个原始候选峰因 `pass_coverage=false` 被剔除。
- 检查全部 shuffle 的染色体、read name、原始起点对应关系一致；同种子
  再次运行位置相同，更换种子结果改变；加入 `chrTest:5000–6000` 排除区域后，
  没有随机分子落入该区域。
- 有效 FDR 表全部 FDR=1 时，`call-peaks --fdr-table` 正常返回零峰结果。
  这与空背景报错属于不同情况。
- `tests/test_fire_peaks.py` 六项测试通过，覆盖按 bp 加权、共同覆盖域、
  同分数累计、空背景失败与加权覆盖中位数。真实配置默认目标 dry-run 通过，
  包含新增的三个样本 peak calling 分支，共 53 个任务。

尚未运行真实大样本的全基因组峰调用，也未证明与完整官方 FIRE workflow
的结果等价。
