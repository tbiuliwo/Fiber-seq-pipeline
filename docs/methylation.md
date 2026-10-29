# 第 08 步：MethBat 5mC 甲基化统计

规则在 `workflow/rules/08_methylation.smk`，读取第 05 步 FIRE BAM 和索引。
不依赖 FiberFold；分相关闭时也可独立运行。本地配置启用，示例配置默认关闭。

## 输出和口径

```text
05_fire/<sample>.fire.bam + .bai
    -> methylation_pileup
    -> 08_methylation/<sample>.5mC.bed.gz + .tbi
    -> 08_methylation/<sample>.summary.json
    -> methylation_bedgraph -> <sample>.5mC.bedgraph.gz
    -> index_methylation_bedgraph -> <sample>.5mC.bedgraph.gz.tbi
```

`methylation` 聚合目标要求每个样本以上五个输出；开启时也加入 `all`。
日志位于 `logs/08_methylation/`。路径遵循主配置的 `results_dir`、`logs_dir`。

MethBat 使用 `--skip-6mA --skip-5hmC`，不要求生成 5hmC 文件。
默认最小 MAPQ 为 1、最小覆盖为 4、两端修剪 20 bp。
绘图轨道选择 `type=Total` 的行，按列名读取 chrom/start/end/mod_score，
保留 0–100 的百分比，不混入单倍型特异行；通过 bgzip 压缩并建立 tabix 索引。
即使以后上游启用分相，这条 bedGraph 仍表示 Total 甲基化。

## 配置

```yaml
threads:
  methylation: 4
methylation:
  enabled: true
  min_mapq: 1
  min_coverage: 4
  edge_trimming_size: 20
```

methbat、bgzip、tabix 默认来自 `tools.bin_dir`，也可以在 methylation 下分别指定路径。
已有科学参数未改变，原结果使用 MethBat 1.1.0-fbf4686 生成。

## 命令示例

以下命令本次未执行：

```bash
snakemake -s workflow/fiberseq.smk methylation -n --cores 4
snakemake -s workflow/fiberseq.smk \
  results/08_methylation/Y5Y-WT-HA.5mC.bedgraph.gz.tbi -n --cores 4
```

旧的 `scripts/run_methbat_5mc.sh` 已删除，直接使用 Snakemake 目标即可。

## 目录迁移

早先 optional/methylation 已改为第 14 步；本次连续编号调整将原来的
`results/14_methylation/`、`logs/14_methylation/` 分别迁移至第 08 步。
原始分析输出和历史日志内容保留，内部记录的旧命令路径作为历史信息不重写。
Figure A 的三个甲基化输入已更新为 `results/08_methylation/`。
此次只修改路径和目录，不执行 pileup、轨道转换、索引或绘图。
