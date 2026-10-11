# Fiber-seq Analysis Pipeline

A reproducible Fiber-seq analysis pipeline built with Snakemake.

## Overview

This project aims to develop a reproducible workflow
for processing and analyzing Fiber-seq data.

## Main Analysis Modules

1. Quality control
2. DNA m6A modification analysis
3. Nucleosome and MSP identification
4. FIRE regulatory element analysis
5. Protein footprint analysis (planned)
6. 3D chromatin prediction (planned)

## Tools

- fibertools-rs
- samtools
- Snakemake
- FIRE

## Project Structure

- `config/`: Configuration files
- `workflow/`: Snakemake workflows and rules
- `scripts/`: Auxiliary scripts
- `docs/`: Documentation
- `data/`: Sequencing data (not tracked by Git)
- `results/`: Analysis results (not tracked by Git)
- `logs/`: Execution logs (not tracked by Git)

## Status

Project initialization and environment setup.
