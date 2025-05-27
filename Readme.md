# ClinVar Variant Analysis for Base Editing

## Description

This repository contains scripts to analyze ClinVar database variants, focusing on:
1. Statistical analysis of pathogenic variant types
2. Distribution of point mutations
3. Evaluation of which point mutations can be corrected by different versions of base editors (ABE8e, ABE8e-YA, ABE9, pcABE)

## Usage

This script reads a VCF file from the ClinVar database and filters for single-nucleotide variants (SNVs) with pathogenic significance. It then designs saturated sgRNAs for these loci and evaluates each sgRNA based on the editing characteristics of different genome editors—such as editing window, editing type, and PAM sequence. The suitability of each sgRNA is assessed accordingly, and the results are exported to a summary table.

> The `vcf` download from URL: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2024/clinvar_20240716.vcf.gz

``` python
python stats.py
```

The output is a tab-delimited text file summarizing the pathogenic SNVs from the ClinVar database that are potentially targetable by various genome editors. Each row represents a candidate sgRNA targeting a specific variant site. The columns are described as follows:

- **CHROM**: Chromosome where the variant is located (e.g., chr1, chr2).

- **POS**: Genomic position (1-based coordinate) of the SNV.

- **Ref**: Reference allele at the variant site.

- **ALT**: Alternative (mutant) allele reported in ClinVar.

- **CLNHGVS**: HGVS notation from ClinVar describing the variant (e.g., NM_000546.5:c.215C>G).

- **GeneName**: Gene symbol associated with the variant .

- **FeatureType**: Type of genomic feature containing the variant (e.g., transcript, exon).

- **HGVS**: Additional HGVS notation of the variant, possibly from transcript annotation.

- **Strand**: DNA strand of the transcript or gene (+ or -).

- **Motif**: Sequence context surrounding the variant.

For the following columns, {Editor} and {PAM} are placeholders that represent specific base editors (e.g., ABE8e, CBE4max) and their corresponding PAM requirements (e.g., NGG, NG):

- **{Editor}_{PAM}**: If multiple sgRNAs are suitable to correct the pathogenic mutation using {Editor} with PAM {PAM}, this column lists the total count.

- **{Editor}\_{PAM}_label**: A comma-separated list indicating whether each sgRNA can correct the mutation within the editable window (`True` for suitable, `False` for not suitable).

- **Is{PAM}Match**: Indicates whether there is at least one sgRNA using the specified {PAM} sequence (TRUE or FALSE).

- **sgrna\_{PAM}_list**: A list of all saturated sgRNAs (whit PAM) designed for this PAM. Each sgRNA is shown with the target mutation in lowercase (e.g., TCTCCTGTGCCGAGTTCTCaCT).

