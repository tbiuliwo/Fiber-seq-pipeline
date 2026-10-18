# CTCF 分支（2026-10-18）

## 入口

```bash
conda activate fiberseq
cd /data/tenghaohan/work/Fiber_seq
snakemake -s workflow/fiberseq.smk prepare_ctcf -n --cores 8
snakemake -s workflow/fiberseq.smk prepare_ctcf --cores 8
```

当前 `ctcf.enabled: true`，所以普通的 `rule all` 也会包含 CTCF 输出。
参考为用户的 GRCh38.p14 FASTA，默认处理 chr1–chr22、chrX。
默认下载作者使用的 JASPAR 2022 hg38 MA0139.1 位点集；不是重新扫描 motif。
可以在 `ctcf.motifs` 填写本地 BED6 或 JASPAR 七列 TSV/gzip 文件。
更换基因组时必须同时更换参考、motif 坐标和染色体列表。

## 依赖和结果

```text
FASTA → results 内的 FAI 索引 ─┐
MA0139.1 位点集 ──────────────┴→ 排序、按分数处理重叠 motif
                                      ├→ direction/<chrom>.bw（所有样本共用）
每个样本 FIRE BAM + BAI ────────────────┴→ ft footprint
                                         → <sample>.footprints.tsv.gz
                                         → <sample>.scores.tsv.gz
                                         → footprint/<sample>/<chrom>.bw
```

所有路径位于 `results/08_ctcf/`。三个分析模块为 [0,6)、[6,13)、[13,19)。
足迹比例采用作者 encodeFP.py 的定义：先要求跨越该 motif 的 MSP，再要求
含 A/T 的模块满足 footprint bit，计数后除以所有跨越该 motif 的分子数。
沿用作者对两条链均使用参考序列顺序确定 A/T 模块资格的编码约定。
方向：正链 +1、负链 -1、其余位置 0。足迹比例已按覆盖分子数归一化，
无需再次除以深度。没有覆盖时 bigWig 填 0，分数表额外保留覆盖数以供区分。

当前 ft 表有 `n_spanning_fires`、`haplotype` 等新增列，脚本按表头读取，
避免作者旧版固定列数解析失败。修正了旧脚本跨染色体比较重叠、缺少 sys
导入和随意跳过最后三条染色体的问题。相邻但不重叠的 motif 不被合并。

## 与FiberFold 预测步骤的连接

FiberFold 的 `ctcf_footprint` 和 `ctcf_direction` 现在自动使用本分支输出。
不需要在 `fiberfold.inputs` 中提供这两个路径。Snakemake 直接跟踪各个
bigWig 文件，轨道不存在时先生成轨道，完成后才调用预测脚本。

后续的 `09_fiberfold_inputs.smk` 已补齐 FIRE/CpG 轨道和 windows 自动生成，
现在FiberFold 预测步骤不再要求 `fiberfold.inputs`。详见 fiberfold_inputs.md。
FiberHMM 分支不作为这里的前置条件；本分支直接调用 `ft footprint`。

## 验证

- `python -m unittest discover -s tests -p 'test_ctcf_tracks.py' -v` 通过。
  覆盖跨染色体和重叠处理、相邻区间、正负方向、零覆盖、无 motif 染色体、
  A/T 模块筛选和按覆盖数计算 0.5/1.0 比例。
- 24 条真实 reads 的临时 BAM 和测试参考：三个样本名称的 CTCF 分支实际运行
  通过，生成全部方向/足迹 bigWig 和表格。测试目录 `/tmp/fiberseq-ctcf-test/`。
- 加入模拟 FIRE/CpG 轨道后，通过 Snakemake 自动生成 CTCF 轨道，再调用
  FiberFold 预测一个窗口，输出有限值的 (1,256,256) 矩阵。
- 真实三样本配置 dry-run 通过。未执行三个大 BAM 的全量分析。

## 方法来源

- https://github.com/altemoselab/fiberFold （README 和 misc 编码脚本）
- https://pmc.ncbi.nlm.nih.gov/articles/PMC11957061/ （Calling CTCF Peaks and CTCF Direction）
- https://fiberseq.github.io/fibertools/extracting/footprint.html
