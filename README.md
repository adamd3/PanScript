# StrainQuant

![StrainQuant](https://github.com/adamd3/StrainQuant/actions/workflows/ci.yml/badge.svg)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**StrainQuant** is a Nextflow pipeline for performing strain-specific bacterial RNA-Seq analysis without a reference genome.

![StrainQuant](docs/images/StrainQuant_pipeline_lowres.png)

## Pipeline summary

The pipeline requires the output from a pan-genome analysis with [`Panaroo`](https://gtonkinhill.github.io/panaroo/) and will perform the following steps:

1. Trim adaptors from reads ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
2. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
3. Pseudo-alignment to strain-specific gene sets ([`kallisto`](https://pachterlab.github.io/kallisto/))
4. Merging and length-scaling of counts for orthologous genes
5. Size-factor scaling of merged counts using the core genome ([`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)), and calculation of normalised TPM (nTPM)
6. Visualisation of gene expression across strains ([`UMAP`](https://umap-learn.readthedocs.io/))

## Installation

You will need to install [`Nextflow`](https://www.nextflow.io/) (version 21.10.3+).

You can run the pipeline as follows:

    nextflow run /path/to/StrainQuant \
        --meta_file /path/to/metadata.txt \
        --gpa_file /path/to/gene_presence_absence.csv \
        --group majority_ST \
        -profile conda -resume

You can run with [`Docker`](https://www.docker.com/) or [`Singularity`](https://sylabs.io/guides/3.5/user-guide/introduction.html) by specifying ` -profile docker` or ` -profile singularity`, respectively.

The `-resume` parameter will re-start the pipeline if it has been previously run.

Explanation of parameters:

- `meta_file`: metadata file (see below).
- `gpa_file`: gene presence/absence file from Panaroo (see below).
- `strandedness`: Strandedness of the reads. 0 = unstranded, 1 = forward-stranded, 2 = reverse-stranded. Default = 2.
- `fragment_len`: Estimated average fragment length for kallisto transcript quantification (only required for single-end reads). Default = 150.
- `fragment_sd`: Estimated standard deviation of fragment length for kallisto transcript quantification (only required for single-end reads). Default = 20.
- `group`: group for plots - must be one of the columns in metadata file.
- `skip_trimming`: do not trim adaptors from reads.
- `outdir`: the output directory where the results will be saved (Default: `./results`).

## Required input

- **Metadata file**: tab-delimited (TSV) file, which must contain at least the following named columns:

  - `rna_sample_id`: name of the RNA-Seq data file
  - `dna_sample_id`: gene sequence sample identifier (must match the column name in the Panaroo gene presence/absence file)
  - `sample_name`: unique identifier for the sample
  - `fastq1`: path to fastq file for RNA-Seq data
  - `fastq2`: (optional) path to 2nd fastq file for paired-end RNA-Seq data - leave blank for single-end data
  - `fasta`: path to fasta file containing gene sequences
  - Additional columns are optional

  Example:

  ```console
  rna_sample_id	dna_sample_id	sample_name	majority_ST	level7000	acsA	aroE	guaA	mutL	nuoD	ppsA	trpE	patient	collection_date	cohort	infection_type	fastq1	fastq2	fasta
  SRX5123742	SRR8737281	ZG301975	313	2	47	8	7	6	8	11	40	ZG301975		hzi_amr		https://raw.githubusercontent.com/adamd3/StrainQuant/main/test_data/SRX5123742_T1_sub.fq.gz		https://raw.githubusercontent.com/adamd3/StrainQuant/main/test_data/SRR8737281.fna
  SRX5123741	SRR8737282	ZG205864	111	3	17	5	5	4	4	4	3	ZG205864		hzi_amr		https://raw.githubusercontent.com/adamd3/StrainQuant/main/test_data/SRX5123741_T1_sub.fq.gz		https://raw.githubusercontent.com/adamd3/StrainQuant/main/test_data/SRR8737282.fna
  SRX5123744	SRR8737283	ZG302367	274	4	23	5	11	7	1	12	7	ZG302367		hzi_amr	ear infection	https://raw.githubusercontent.com/adamd3/StrainQuant/main/test_data/SRX5123744_T1_sub.fq.gz		https://raw.githubusercontent.com/adamd3/StrainQuant/main/test_data/SRR8737283.fna
  ```

- **Gene presence-absence file**: CSV-format output produced by [`Panaroo`](https://gtonkinhill.github.io/panaroo/).
  See below example (truncated):

  ```console
  Gene,Non-unique Gene name,Annotation,Pseudomonas_aeruginosa_PAO1_107_converted,SRR8737281,SRR8737282
  group_25769,,Protein of unknown function (DUF2845),PGD112333,SRR8737281_00662,SRR8737282_00657
  group_25760,;hypothetical protein,hypothetical protein;topoisomerase IIProtein of unknown function (DUF2790);Protein of unknown function (DUF2790);acyl- n-acyltransferaseUncharacterized protein conserved in bacteriaProtein of unknown function DUF482;,PGD109272,SRR8737281_05540,SRR8737282_06520
  group_25759,;hypothetical protein,Protein of unknown function (DUF3309);,PGD109276,SRR8737281_05538,SRR8737282_03481
  rubA2~~~Rubredoxin2,rubA2;Rubredoxin 2,rubredoxinRubredoxin-2anaerobic nitric oxide reductase flavorubredoxinRubredoxinRubredoxin;,PGD113575,SRR8737281_00021,SRR8737282_05845
  bfr_2~~~bacterioferritin,bfr_2;bacterioferritin,bacterioferritinBacterioferritinbacterioferritinBacterioferritin (cytochrome b1)bacterioferritinFerritin-like domain;,PGD109882,SRR8737281_05955,SRR8737282_01475
  mtnB~~~probablesugaraldolase,mtnB;probable sugar aldolase,methylthioribulose-1-phosphate dehydrataseMethylthioribulose-1-phosphate dehydratasemethylthioribulose-1-phosphate dehydrataseRibulose-5-phosphate 4-epimerase and related epimerases and aldolasesmethylthioribulose-1-phosphate dehydrataseClass II Aldolase and Adducin N-terminal domain;,PGD106137,SRR8737281_02504,SRR8737282_00940
  epd_2~~~epd_1~~~epd,epd_2;epd_1;epd,D-erythrose-4-phosphate dehydrogenaseD-erythrose-4-phosphate dehydrogenaseerythrose 4-phosphate dehydrogenaseerythrose-4-phosphate dehydrogenaseGlyceraldehyde 3-phosphate dehydrogenase C-terminal domain;D-erythrose-4-phosphate dehydrogenaseD-erythrose-4-phosphate dehydrogenaseerythrose 4-phosphate dehydrogenaseTransketolaseerythrose-4-phosphate dehydrogenaseGlyceraldehyde 3-phosphate dehydrogenase C-terminal domain,PGD103840,SRR8737281_03107,SRR8737282_02473
  emrE,emrE,SMR multidrug efflux transporterMethyl viologen resistance protein Cmultidrug efflux proteinMembrane transporters of cations and cationic drugsphosphonate utilization associated putative membrane proteinSmall Multidrug Resistance protein;SMR multidrug efflux transporterMethyl viologen resistance protein Cmultidrug efflux proteinMembrane transporters of cations and cationic drugsSmall Multidrug Resistance protein,PGD112851,SRR8737281_00400,SRR8737282_00337
  purA_1~~~purA_2~~~purA,purA_1;purA_2;purA,adenylosuccinate synthetaseAdenylosuccinate synthetaseadenylosuccinate synthetaseadenylosuccinate synthaseAdenylosuccinate synthetase,PGD112747,SRR8737281_00453,SRR8737282_00390
  ```

## Output

1. **trim_galore** directory containing adaptor-trimmed RNA-Seq files and FastQC results.
2. **gene_counts** directory containing three gene × sample matrices (one row per Panaroo gene group). All values are on a **linear scale** (not log-transformed), and genes absent from an isolate are `NA`:
   1. `raw_counts.tsv`: read counts adjusted **only** for ortholog length differences between isolates (scaled to the gene group's median length). Not normalised for library size.
   2. `norm_counts.tsv`: `raw_counts` divided by DESeq2 size factors estimated on the core genome. Use these to compare **the same gene across samples**.
   3. `ntpm.tsv`: normalised transcripts per million (nTPM). Counts per unit of each isolate's own gene length, divided by the core-genome size factors and scaled so per-sample totals average 10<sup>6</sup>. Comparable **between genes and across samples**. Unlike standard TPM, it's not distorted by differences in accessory gene content between isolates.

   See [docs/normalisation.md](docs/normalisation.md) for formulas and rationale.
3. **kallisto_quant** directory containing per-sample kallisto output and the merged count and effective-length matrices.
4. **umap_samples** directory containing a UMAP of samples, computed on log2(`norm_counts` + 1) of core genes.
