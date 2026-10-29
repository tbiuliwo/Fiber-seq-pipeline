# 第 02 步：可选变异检测、分相与单倍型标记

规则在 `workflow/rules/02_phasing.smk`，位于比对之后、Fiber-seq 注释之前。
本地配置和示例配置均默认 `phasing.enabled: false`。

## 数据流

```text
01_align/<sample>.chromosomes.bam + .bai
    -> DeepVariant PACBIO -> <alias>.deepvariant.vcf.gz + .tbi
    -> WhatsHap phase -> <alias>.phased.vcf.gz + .tbi
        -> WhatsHap stats -> 统计 TSV、分相区块 GTF
        -> WhatsHap haplotag（输入仍为第 01 步 BAM）
        -> 02_phasing/<alias>/<alias>.haplotagged.bam + .bai
        -> 03_annotate -> 04_qc、05_fire 及其他下游
```

`downstream_aligned_bam()` 负责切换输入：关闭分相时返回第 01 步 BAM；
开启时，选中的样本返回第 02 步带单倍型标签的 BAM，未选中的样本仍走原路径。
`input_mode: annotated` 时跳过重新注释，也通过这个函数选择下游输入。
比对质量的 flagstat/stats 仍统计第 01 步 BAM，Fiber-seq QC 则读取下游注释 BAM。
分相不依赖 FIRE，避免形成循环依赖。

## 配置

```yaml
phasing:
  enabled: false
  container: /path/to/deepvariant_1.10.0.sif
  singularity: singularity
  bind_paths: [/data:/data]
  reference_fasta: null
  mem_mb: 64000
  whatshap: whatshap
  whatshap_mem_mb: 16000
  samtools: null
  tabix: null
  samples: null
  sample_aliases:
    Y5Y-WT-HA: WT
    Y5Y-R354L-HA: R354L
    Y5Y-4A-HA: 4A
```

`reference_fasta: null` 复用主参考；必须是有 `.fai` 的 FASTA，不能是 `.mmi`。
`samples: null` 选择所有样本，列表可使用样本全名或输出别名。
别名只影响输出目录，不会替换 BAM/VCF 中的真实样本名。
WhatsHap 根据 BAM read group 的 SM 字段匹配 VCF 样本。

DeepVariant 使用 `threads.deepvariant`；phase/stats 每样本单线程；
haplotag 使用 `threads.haplotag`，剩余线程用于压缩。
`samtools` 和 `tabix` 默认使用 `tools.bin_dir`。
WhatsHap 可安装在分析环境，或以后使用 `--sdm conda` 管理
`workflow/envs/phasing.yaml`。这次没有安装或运行任何分析软件。

## 输出

所有新结果保存在 `results/02_phasing/<alias>/`：

- `<alias>.deepvariant.vcf.gz` 和 `.tbi`：变异检测结果。
- `<alias>.phased.vcf.gz` 和 `.tbi`：分相结果。
- `<alias>.haplotagged.bam` 和 `.bai`：第 01 步 BAM 加上可分配的 HP/PS 标签。
- `<alias>.phasing.stats.tsv`、`<alias>.phasing.blocks.gtf`：分相统计和区块。
- `intermediate/`、`logs/`：保留的 DeepVariant 中间文件和工具日志。

主流程日志位于 `logs/02_phasing/`。
无法可靠分型的 reads 仍保留。H1/H2 是分相区块内的标识，不直接代表父源/母源，
也不能默认在不同区块或独立样本间对应。

## 以后启用时

先设置 `phasing.enabled: true`，再按需要选择目标。下列仅为命令示例，本次未执行：

```bash
snakemake -s workflow/fiberseq.smk phasing -n --cores 16
snakemake -s workflow/fiberseq.smk variant_calls -n --cores 16
snakemake -s workflow/fiberseq.smk \
  results/02_phasing/4A/4A.haplotagged.bam.bai -n --cores 16
```

给已有结果目录改名不会补上单倍型信息。以后启用分相时，应检查新 BAM 的
HP/PS、MM/ML 及后续注释/FIRE BAM 的标签保留情况，再生成按单倍型的下游结果。
原来的 `13_phasing` 已迁移至 `02_phasing`；现有结果只有 DeepVariant VCF、索引
和中间文件，本次没有生成 phased VCF 或带单倍型标签的 BAM。
旧的 `scripts/analysis/13_phasing.sh` 已删除，不再作为入口。
