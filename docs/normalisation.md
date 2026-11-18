# How expression values are normalised

This document explains how the pipeline turns kallisto output into the expression matrices it
reports, and why. It's the reference for the README's Output section and for the manuscript
Methods.

> **Status (2026-11-18).** Sections marked *[planned]* describe agreed changes that are not
> implemented yet (see PROGRESS.md: P2–P6). Everything else describes current behaviour.

## The problem

Each isolate is quantified against its **own** gene set, and genes are then matched across isolates
through Panaroo orthologous groups. That causes three complications that single-reference RNA-Seq
doesn't have:

1. **Ortholog length varies.** The same gene group can differ in length between isolates, so raw
   counts for a group aren't comparable across isolates even at equal expression.
2. **Gene content varies.** Accessory genes are absent (NA) in some isolates. Any library-size
   measure that sums over "all genes" is therefore summing over a *different set of genes* in each
   sample.
3. **Genes can be split.** In some isolates a gene group is represented by 2+ loci
   (`;`-separated in the Panaroo table). In others it's a single locus.

## Notation

For gene group *g* in sample *s* (from isolate *i(s)*):

- *c<sub>gs</sub>*: kallisto `est_counts` for that isolate's locus or loci in group *g*
- *L<sub>gs</sub>*: kallisto `eff_length` for that locus or loci
- *M<sub>g</sub>* = median over samples of *L<sub>gs</sub>* (ignoring NA): median length
- *sf<sub>s</sub>*: DESeq2 size factor for sample *s*

## Steps

### 1. Split genes
When a gene group maps to several loci in an isolate, their counts and effective lengths are
combined. Currently both are **averaged**. The per-isolate rate *c/L* is correct either way, because
mean/mean = sum/sum. But averaging makes *L<sub>gs</sub>* (and so the median *M<sub>g</sub>*) the
length of a *fragment* rather than of the whole gene.
*[planned, P5]* Both will be **summed**, so *L<sub>gs</sub>* is the total length of the gene in that
isolate.

### 2. Ortholog length scaling → `raw_counts.tsv`
*r<sub>gs</sub>* = *c<sub>gs</sub>* / *L<sub>gs</sub>* × *M<sub>g</sub>*

This is the number of reads the gene would have received if it had the group's median length in
every isolate. It removes length differences **between isolates for the same gene**. It is **not**
normalised for library size, and **not** comparable between different genes.

### 3. Size factors from the core genome
DESeq2 median-of-ratios size factors *sf<sub>s</sub>* are estimated from the matrix of rounded
*r<sub>gs</sub>*, restricted to genes with **no NA in any sample** (the core genome). The size factors are then applied to every non-NA value in that sample, accessory genes included.

The reasoning: median-of-ratios needs the same genes in every sample. The core genome usually has
thousands of genes, which is plenty for a robust estimate. And median-of-ratios isn't thrown off
by a few highly expressed or strongly changing genes.

### 4. Library-size normalisation → `norm_counts.tsv`
*n<sub>gs</sub>* = *r<sub>gs</sub>* / *sf<sub>s</sub>*

These values are corrected for ortholog length variation and library size, so they're comparable
**across samples for the same gene**. They're still not comparable between genes, because longer
genes get more reads. This is the recommended matrix for between-sample analyses of a given gene
(clustering, ordination, plotting one gene across isolates).

### 5. Between-gene comparable expression → nTPM *[planned, P2]*
This replaces the current `rpkm_counts.tsv`.

*x<sub>gs</sub>* = *c<sub>gs</sub>* / *L<sub>gs</sub>* / *sf<sub>s</sub>* (the same as
*n<sub>gs</sub>* / *M<sub>g</sub>*, because the median length cancels)

nTPM<sub>gs</sub> = *K* · *x<sub>gs</sub>*, with *K* = 10<sup>6</sup> / mean<sub>s</sub>( Σ<sub>g</sub> *x<sub>gs</sub>* )

Here the sum is over the genes present (non-NA) in each sample. *K* is one constant for the whole
dataset. It only puts values on a familiar TPM-like scale: per-sample totals average 10<sup>6</sup>.

**Why not ordinary TPM?** TPM divides each sample by its *own* total, Σ<sub>g</sub>
*c<sub>gs</sub>*/*L<sub>gs</sub>*, which is summed over that isolate's genes. Isolates have
different accessory genes, so that total depends on gene content. An isolate carrying many (or
highly expressed) accessory genes would get smaller TPMs for *every* core gene, even with
identical core expression. Using the core-genome size factors avoids this.

**What nTPM can be compared across:**
- *Between genes within a sample*: yes. Within a sample, the value is proportional to reads per
  unit of the isolate's own transcript length, exactly as TPM is.
- *Across samples*: yes. Library size is handled by the core-genome size factors, so it doesn't
  depend on gene content.
- *Difference from strict TPM*: per-sample totals are **approximately**, not exactly, 10<sup>6</sup>.

**Caveats** (the same as for any TPM, plus one specific to pangenomes): reads that map to several
similar genes (paralogs); how accurately effective lengths are estimated (kallisto subtracts about
one fragment length per locus, so for split genes the summed effective length is slightly
shorter than the intact gene); sequence and position bias; and sequence divergence within an
ortholog group.

**Name:** "nTPM" (normalised TPM), following the Human Protein Atlas's use of the term for
TPM normalised across samples. Our size-factor method should be stated explicitly wherever the
term is used.

### 6. Output scale *[planned, P6]*
All output matrices will be **linear** (not log-transformed). Currently `norm_counts.tsv` and
`rpkm_counts.tsv` are log2(*x*+1). Users should apply their own transform (such as log2(*x*+1) or
DESeq2's VST) for analysis.

**UMAP requirement:** `bin/umap.R` currently runs on the log2 `norm_counts` (core genes only, since
NA rows are dropped). Once outputs are linear, `umap.R` **must apply log2(*x*+1) itself** so the
plot is unchanged.

## Output summary

| File | Ortholog length | Library size | Between-gene comparable | Scale |
|------|-----------------|--------------|-------------------------|-------|
| `raw_counts.tsv` | ✓ (scaled to median length) | ✗ | ✗ | counts (non-integer) |
| `norm_counts.tsv` | ✓ | ✓ (core-genome DESeq2) | ✗ | counts; currently log2, *[planned]* linear |
| nTPM *[planned, replaces `rpkm_counts.tsv`]* | ✓ (own effective length) | ✓ (core-genome DESeq2) | ✓ | TPM-like, linear |

Absent genes are NA in all matrices.
