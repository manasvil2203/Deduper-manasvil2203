# PCR Duplicate Removal from SAM Files Using UMIs

## Overview

This program removes PCR duplicates from a sorted SAM file using Unique Molecular Identifiers (UMIs). Reads are considered duplicates if they map to the same genomic position on the same strand and share the same UMI sequence. Only one representative read is retained, and all subsequent duplicates are removed.

The script outputs a deduplicated SAM file and a summary report describing the filtering process.

---

## Methodology

Each read is uniquely identified using the following combination:

- Chromosome  
- 5′ mapping position  
- Strand  
- UMI sequence  

Reads with identical values for all four components are considered PCR duplicates. The first occurrence is kept, and later duplicates are discarded.

Duplicate tracking is reset when the chromosome changes to reduce memory usage.

---

## Input Requirements

### 1. Sorted SAM File

The input SAM file must be sorted by chromosome and genomic position (RNAME and POS fields). This is required because duplicate tracking is reset when the chromosome changes. If the file is not sorted, some duplicates may not be detected.

Both uncompressed (`.sam`) and compressed (`.sam.gz`) formats are supported.

### 2. UMI List File

The UMI file must contain one valid UMI per line. Blank lines are ignored.

Example:
ATCG
GCTA
TTAA
CGAT


Only reads with UMIs present in this list are retained.

---

## Output Files

### 1. Deduplicated SAM File

The output SAM file contains:

- All original header lines  
- Only unique, non-duplicate reads  

### 2. Summary Report

The summary file reports:

- Number of header lines  
- Total reads processed  
- Number of unique reads retained  
- Number of reads skipped due to invalid UMIs  
- Number of PCR duplicates removed  
- Number of retained reads per chromosome  

---

## Usage

```bash
python deduper.py \
  -f input.sam.gz \
  -o output.sam \
  -u umi_list.txt \
  -s summary.txt
```
