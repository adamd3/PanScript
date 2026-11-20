# PanScript: usage

- [Quick start](#quick-start)
- [Inputs](#inputs)
- [Parameters](#parameters)
- [Software profiles](#software-profiles)
- [Running on a cluster](#running-on-a-cluster)

See [output.md](output.md) for what the pipeline produces, and
[normalisation.md](normalisation.md) for how the expression values are calculated.

## Quick start

You need [Nextflow](https://www.nextflow.io/) 24.10.0 or later and one of Docker, Singularity,
Apptainer or conda.

    nextflow run adamd3/PanScript \
        --meta_file metadata.tsv \
        --gpa_file gene_presence_absence.csv \
        --outdir results \
        -profile docker

Check your setup with the small test dataset first:

    nextflow run adamd3/PanScript -profile test,docker --outdir test_results

`nextflow run adamd3/PanScript --help` lists every parameter. Parameter values are checked before
the run starts (types, allowed values and that input files exist), and misspelt parameter names are
reported. Add `-resume` to restart a run from where it stopped, reusing finished steps.

To run a specific release, add `-r <version>` (for example `-r v2.0.0`). This is recommended for
analyses you'll publish, so the exact pipeline version is recorded.

## Inputs

You need three things for each isolate: its genes (a FASTA of gene sequences), its entry in a
Panaroo pangenome, and its RNA-seq reads.

### Preparing inputs from a Panaroo run

The pipeline needs one gene FASTA per isolate, with locus IDs matching Panaroo's
`gene_presence_absence.csv`. `scripts/prepare_inputs_from_panaroo.py` builds these from a standard
Panaroo output directory. It uses `gene_data.csv`, which holds the DNA sequence of every annotated
and refound gene, so refound pseudogenes are included. It also writes a metadata template:

    python scripts/prepare_inputs_from_panaroo.py \
        --panaroo_dir /path/to/panaroo_output --outdir pipeline_inputs

Both Panaroo files can be gzipped (`gene_data.csv.gz`, `gene_presence_absence.csv.gz`).
`gene_data.csv` is large (about 10 GB compressed for 3,000 genomes) but is streamed, so memory use
stays small. Add `--isolates ID1 ID2 ...` to write FASTAs only for the isolates you have RNA-seq
data for; it still reads the whole file (a few minutes for 3,000 genomes).

This creates `pipeline_inputs/fasta/<isolate>.fna` and `pipeline_inputs/metadata_template.tsv`
(with `dna_sample_id` and `fasta` filled in). Complete the template with `rna_sample_id`,
`sample_name`, `fastq1` and `fastq2` (one row per RNA-seq sample), and pass Panaroo's
`gene_presence_absence.csv` as `--gpa_file`.

### Panaroo's locus ID suffixes

In `gene_presence_absence.csv`, Panaroo adds a suffix to some locus IDs to flag possible problems
with that gene in that isolate:

- `_len`: the gene's length differs markedly from the usual length in its gene group (for example
  a truncated or extended gene).
- `_pseudo`: the sequence isn't a valid gene (for example it has internal stop codons), which
  usually means a pseudogene.
- `_stop` (older Panaroo versions): a premature stop codon. It can be combined with `_len`, as in
  `_len_stop`.

These suffixes aren't part of the gene's sequence ID in the isolate FASTA, so the pipeline matches
each locus on its exact ID first, then with the suffixes removed. Flagged genes are therefore
quantified like any other gene. If you'd rather exclude them, filter them out downstream: the
suffixes are visible in `gene_presence_absence.csv` (`gene_presence_absence_roary.csv` and
`gene_presence_absence.Rtab` list the same loci without suffixes).

### Input file formats

- **Metadata file**: tab-delimited (TSV) file, which must contain at least the following named columns:

  - `rna_sample_id`: name of the RNA-Seq data file
  - `dna_sample_id`: gene sequence sample identifier (must match the column name in the Panaroo gene presence/absence file)
  - `sample_name`: unique identifier for the sample
  - `fastq1`: path to fastq file for RNA-Seq data
  - `fastq2`: (optional) path to 2nd fastq file for paired-end RNA-Seq data - leave blank for single-end data. All samples in a run must be the same type: the pipeline stops if some samples have `fastq2` and others don't
  - `fasta`: path to fasta file containing gene sequences
  - Additional columns are optional

  One row per RNA-seq sample. Several samples can come from the same isolate (replicates,
  conditions or time points): give them the same `dna_sample_id` and the same `fasta`. The
  pipeline builds one kallisto index per isolate and shares it between that isolate's samples,
  and each isolate counts once in the median gene length used for length scaling (see
  [normalisation.md](normalisation.md)). The run stops if a `sample_name` is repeated or if one
  `dna_sample_id` is given more than one `fasta`.

  Example:

  ```console
  rna_sample_id	dna_sample_id	sample_name	majority_ST	level7000	acsA	aroE	guaA	mutL	nuoD	ppsA	trpE	patient	collection_date	cohort	infection_type	fastq1	fastq2	fasta
  SRX5123742	SRR8737281	ZG301975	313	2	47	8	7	6	8	11	40	ZG301975		hzi_amr		https://raw.githubusercontent.com/adamd3/PanScript/main/test_data/SRX5123742_T1_sub.fq.gz		https://raw.githubusercontent.com/adamd3/PanScript/main/test_data/SRR8737281.fna
  SRX5123741	SRR8737282	ZG205864	111	3	17	5	5	4	4	4	3	ZG205864		hzi_amr		https://raw.githubusercontent.com/adamd3/PanScript/main/test_data/SRX5123741_T1_sub.fq.gz		https://raw.githubusercontent.com/adamd3/PanScript/main/test_data/SRR8737282.fna
  SRX5123744	SRR8737283	ZG302367	274	4	23	5	11	7	1	12	7	ZG302367		hzi_amr	ear infection	https://raw.githubusercontent.com/adamd3/PanScript/main/test_data/SRX5123744_T1_sub.fq.gz		https://raw.githubusercontent.com/adamd3/PanScript/main/test_data/SRR8737283.fna
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

## Parameters

- `meta_file`: metadata file (see below).
- `gpa_file`: gene presence/absence file from Panaroo (see below).
- `strandedness`: Strandedness of the reads. 0 = unstranded, 1 = forward-stranded, 2 = reverse-stranded, `auto` = inferred separately for each sample by the strandedness check (see `strandedness_check_reads`; the run stops if a sample cannot be classified). Default = `auto`.
- `fragment_len`: Estimated average fragment length for kallisto transcript quantification (only required for single-end reads). Default = 150.
- `fragment_sd`: Estimated standard deviation of fragment length for kallisto transcript quantification (only required for single-end reads). Default = 20.
- `strandedness_check_reads`: the pipeline checks each library's strandedness by pseudo-aligning this many reads (the first N) in both orientations, and **warns** if the result disagrees with `--strandedness`. Results are in `pipeline_info/strandedness_check.tsv`. A wrong setting can discard most reads: with the test data, the reverse setting pseudo-aligns about 5% of reads versus about 50% with the correct forward setting. 0 skips the check. Default = 200000.
- `max_missing_loci`: the pipeline checks that every locus Panaroo lists for an isolate is in that isolate's FASTA. Panaroo's status suffixes on locus IDs (`_len`, `_pseudo`, `_stop`) are ignored when matching, because they aren't part of the sequence ID. Missing loci can't be quantified, so their gene groups are `NA` for that sample. They're always listed in `pipeline_info/missing_loci.tsv` (with a per-sample summary in `pipeline_info/fasta_loci_summary.tsv`), and the run **stops** if more than this fraction of an isolate's loci are missing, which usually means the wrong FASTA or a mismatched `dna_sample_id`. Default = 0.5.
- `group`: group for plots - must be one of the columns in metadata file.
- `skip_trimming`: do not trim reads (skip fastp).
- `fastp_args`: extra options for fastp, which otherwise runs with its **default settings**. These include adapter trimming, discarding reads with more than 40% of bases below Q15 or with too many Ns, a 15 bp minimum length, and automatic polyG trimming for NovaSeq/NextSeq data. To approximate the Trim Galore defaults used before 2026 (3′ quality trimming at Q20, 20 bp minimum length), use `--fastp_args '--cut_tail --cut_tail_mean_quality 20 --length_required 20'`. Default = none.
- `save_trimmed`: save the trimmed reads in `<outdir>/fastp/trimmed/`. Default = false.
- `seed`: random seed for the UMAP, so the plot is reproducible. Default = 42.
- `outdir`: the output directory where the results will be saved (Default: `./results`).
- `publish_dir_mode`: how results are saved in `outdir`. `copy` (the default) makes independent files, so you can delete Nextflow's `work/` directory afterwards. `symlink` saves disk space but the results then break if `work/` is deleted. Other options: `rellink`, `link`, `copyNoFollow`, `move`.
- `skip_multiqc`: don't make the MultiQC report. Default = false.

## Software profiles

Choose how the software is provided with `-profile docker`, `-profile singularity`,
`-profile apptainer` or `-profile conda`. The container profiles all use the same versioned image
(`adamd3/panscript`), whose tool versions are pinned in `environment.yml`. Singularity and Apptainer
cache images in the directory given by `NXF_SINGULARITY_CACHEDIR` / `NXF_APPTAINER_CACHEDIR`.

The `test` profile adds the small test dataset; combine it with a software profile
(`-profile test,docker`).

## Running on a cluster

The pipeline has no built-in site profiles. Put your cluster's settings in a small config file and pass it with `-c`. For example, for SLURM:

    // my_cluster.config
    process {
        executor       = 'slurm'
        queue          = 'my_partition'
        clusterOptions = '--account=my_account'
        resourceLimits = [cpus: 32, memory: 128.GB, time: 24.h]   // your largest node
    }
    executor.queueSize = 100

    nextflow run adamd3/PanScript -profile singularity -c my_cluster.config ...

By default each kallisto task asks for 8 CPUs and 8 GB of memory and the other steps for less. On a smaller machine, cap the requests with `process.resourceLimits` in the same way (for example `[cpus: 4, memory: 8.GB]`).

Your institution may already have a Nextflow config (for example in [nf-core/configs](https://github.com/nf-core/configs)); its `process`, `executor` and container settings work here too.
