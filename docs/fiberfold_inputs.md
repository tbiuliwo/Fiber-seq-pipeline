# FiberFold 自动输入准备

`09_fiberfold_inputs.smk` 已补齐 FIRE bigWig、CpG bigWig 和预测窗口。
FiberFold 预测步骤直接依赖这些输出，以及CTCF 模块生成的 两类轨道。

## 使用

```yaml
fiberfold:
  enabled: false
  prepare_inputs: true
  window_step: 2097152
```

- `prepare_inputs: true`：默认任务包含四类轨道和窗口，不执行预测。
- `enabled: true`：自动准备输入并执行预测，不需要手工填轨道路径。
- 染色体统一使用 `ctcf.chromosomes`；参考使用 `ctcf.reference_fasta`，
  必须与比对参考一致。不要再填写旧的 `fiberfold.inputs`。

```bash
conda activate fiberseq
snakemake -s workflow/fiberseq.smk prepare_fiberfold_inputs -n --cores 8
snakemake -s workflow/fiberseq.smk prepare_fiberfold_inputs --cores 8
# 将 enabled 改为 true 后，可以直接执行主流程或指定某个预测文件作为目标。
```

## 计算方法

### FIRE

使用 `ft pileup` 的 `fire_coverage / coverage`，不是 pileup 的加权 `score`
列，也不是只把单分子 FIRE BED 转换格式。保留零覆盖区间，输出范围 0–1。
已有覆盖但没有 FIRE 时为 0，无覆盖区域也填 0。输出按染色体流式写入，
合并相邻相同值，不构建整条人染色体的 Python 数组，也不保存巨大中间表。
当前采用 ft 默认的分子过滤设置，没有额外启用 `--fire-filter`。

### CpG

在 fiberseq 环境安装 `pb-cpg-tools=3.0.0`。命令使用：

```text
aligned_bam_to_cpg_scores --pileup-mode model --modsites-mode reference
                         --min-coverage 4 --bam ... --ref ...
```

参考模式与 FiberFold 论文相同。v3 使用内置模型，不再传入旧版的 --model。
使用已有 MM/ML 中的 C+m 标签，不把 C+h 的 5hmC 合并到 5mC。
该工具的输出 bigWig 是 0–100 的百分数；转换后除以 100 得到 0–1，
不能再除一次覆盖度。模型硬编码最低覆盖为 4，低覆盖/未报告位置填 0，
这与上游 FiberFold loader 将缺失值转换为 0 的处理相符。

本模块不替代可选模块的独立 5mC/5hmC 分析；两者可以分别开启。
参考 FASTA 通过结果目录内的软链接与 FAI 配对，不复制整个基因组，
也不向源参考目录写索引。

### 预测窗口

长度固定 2,097,152 bp，默认步长也为 2,097,152 bp，从每条染色体的 0 开始。
仅输出完整窗口，末端不足长度的区域跳过，输出 gzip BED 三列。
步长可以修改；窗口长度不能随意修改，因为预训练模型要求固定长度。
当前没有根据 blacklist、N 区域或覆盖度筛选窗口。低覆盖、组装空缺窗口的
预测需要后续质量评估；这不是已经过质量筛选的生物学结论。

## 输出

```text
results/09_fiberfold_inputs/
  windows.bed.gz
  chromosomes.bed
  reference/genome.fa                 # 软链接
  reference/genome.fa.fai
  fire/<sample>/<chrom>.bw
  cpg_scores/<sample>.combined.bed.gz
  cpg_scores/<sample>.combined.bed.gz.tbi
  cpg_scores/<sample>.combined.bw       # 工具原始百分数
  cpg/<sample>/<chrom>.bw              # 模型使用的 0–1 轨道
```

CTCF 输入仍在 `results/08_ctcf/`，最终预测在 `results/10_fiberfold/`。

## 验证

- 数值测试：FIRE 2/4 输出 0.5，零覆盖和背景为 0；CpG 0/50/100 转为
  0/0.5/1；缺少某条染色体的 CpG 记录时仍生成零背景轨道。
- 窗口测试：长度固定、末端不补短窗口。
- 真实 reads 派生的小样本（重复扩充覆盖，仅用于技术测试）：从 BAM 自动
  生成四类轨道和三个窗口，运行模型得到 3 个有限值的 (1,256,256) 矩阵。
- 逐位点核对了 pb-CpG-tools 的 10 个实际输出位点，百分数换算正确。
- 未执行三个大样本的全量分析；这些测试不验证预测的生物学准确性。

方法参考：

- https://github.com/altemoselab/fiberFold
- https://pmc.ncbi.nlm.nih.gov/articles/PMC11957061/
- https://github.com/PacificBiosciences/pb-CpG-tools
- 官方 FIRE 的 percent_accessible 规则：
  https://github.com/fiberseq/FIRE/blob/main/workflow/rules/track-hub.smk
