#!/usr/bin/env python3
"""Draw the PanScript pipeline diagram (docs/images/panscript_pipeline.svg).

Run from the repo root:
    python3 docs/images/make_pipeline_diagram.py
    rsvg-convert -z 2 docs/images/panscript_pipeline.svg -o docs/images/panscript_pipeline.png
"""
from pathlib import Path

W, H = 1180, 1390
FS = 1.3                             # font scale
GEN, RNA, MRG, GREY = "#1f7a8c", "#e4572e", "#3d405b", "#8a8f98"
XL, XR, XM = 380, 760, 570           # genome track, RNA track, merged track
FONT = "Helvetica, Arial, sans-serif"
out = []


def add(s):
    out.append(s)


def text(x, y, lines, anchor="start", size=13, weight="normal", colour="#222", mono=False):
    size = size * FS
    fam = "Menlo, Consolas, monospace" if mono else FONT
    for i, (t, *style) in enumerate(lines):
        w = style[0] if style else weight
        c = style[1] if len(style) > 1 else colour
        s = style[2] * FS if len(style) > 2 else size
        add(f'<text x="{x}" y="{y + i * (s + 5)}" font-family="{fam}" font-size="{s}" '
            f'font-weight="{w}" fill="{c}" text-anchor="{anchor}">{t}</text>')


def file_icon(x, y, label, colour):
    """Document icon centred on (x, y) with a label band."""
    w, h, f = 46, 56, 13
    x0, y0 = x - w / 2, y - h / 2
    add(f'<path d="M{x0},{y0} h{w - f} l{f},{f} v{h - f} h{-w} z" fill="#fff" '
        f'stroke="{colour}" stroke-width="2.5" stroke-linejoin="round"/>')
    add(f'<path d="M{x0 + w - f},{y0} v{f} h{f}" fill="none" stroke="{colour}" stroke-width="2"/>')
    add(f'<rect x="{x0 + 1.5}" y="{y0 + 30}" width="{w - 3}" height="15" fill="{colour}"/>')
    add(f'<text x="{x}" y="{y0 + 41.5}" font-family="Menlo, Consolas, monospace" font-size="12" '
        f'fill="#fff" text-anchor="middle">{label}</text>')


def station(x, y, colour, r=11):
    add(f'<circle cx="{x}" cy="{y}" r="{r}" fill="#fff" stroke="{colour}" stroke-width="4"/>')


def line(d, colour, width=7, dash=None):
    da = f' stroke-dasharray="{dash}"' if dash else ""
    add(f'<path d="{d}" fill="none" stroke="{colour}" stroke-width="{width}" '
        f'stroke-linecap="round"{da}/>')


add(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">')
add(f'<rect width="{W}" height="{H}" fill="#fff"/>')

# ---- PanScript box -------------------------------------------------------
BOX_TOP, BOX_BOT = 395, 1215
add(f'<rect x="40" y="{BOX_TOP}" width="{W - 80}" height="{BOX_BOT - BOX_TOP}" rx="18" '
    f'fill="#f5f6f8" stroke="#c9ccd3" stroke-width="1.5"/>')
text(64, BOX_TOP + 34, [("PanScript", "bold", MRG, 22)])
text(64, BOX_TOP + 56, [("Nextflow pipeline", "normal", GREY, 12.5)])

# ---- track lines (drawn first, under the icons) ----------------------------
Y_ANN, Y_PAN, Y_FA, Y_FQ = 70, 190, 310, 70
Y_CHK, Y_IDX, Y_QNT = 475, 585, 705
Y_MRG, Y_LEN, Y_SF, Y_NTPM, Y_END = 815, 905, 995, 1085, 1165
Y_DOWN = 1300

line(f"M{XL},{Y_ANN} V{Y_IDX}", GEN)
line(f"M{XR},{Y_FQ} V{Y_IDX}", RNA)
line(f"M{XL},{Y_IDX} C{XL},{Y_IDX + 75} {XM},{Y_QNT - 75} {XM},{Y_QNT}", GEN)
line(f"M{XR},{Y_IDX} C{XR},{Y_IDX + 75} {XM},{Y_QNT - 75} {XM},{Y_QNT}", RNA)
line(f"M{XM},{Y_QNT} V{Y_END}", MRG)
line(f"M{XM},{Y_END} V{Y_DOWN}", GREY, width=4, dash="2 10")

# ---- upstream: genomes -> Panaroo -> per-isolate gene sets ------------------
text(XL, 22, [("Genomes (per isolate)", "bold", GEN, 14)], anchor="middle")
text(XR, 22, [("Transcriptomes (per sample)", "bold", RNA, 14)], anchor="middle")

file_icon(XL, Y_ANN, "gff", GEN)
text(XL - 40, Y_ANN - 4, [("Annotated genomes", "bold"),
                         ("Bakta or Prokka", "normal", "#555", 12.5)], anchor="end")
file_icon(XL, Y_PAN, "csv", GEN)
text(XL - 40, Y_PAN - 12, [("Pangenome", "bold"),
                          ("Panaroo: gene groups", "normal", "#555", 12.5),
                          ("(gene_presence_absence.csv,", "normal", "#555", 12.5),
                          ("gene_data.csv)", "normal", "#555", 12.5)], anchor="end")
file_icon(XL, Y_FA, "fasta", GEN)
text(XL - 40, Y_FA - 12, [("Isolate gene sets", "bold"),
                         ("each isolate's own gene", "normal", "#555", 12.5),
                         ("sequences, from Panaroo", "normal", "#555", 12.5),
                         ("(prepare_inputs_from_panaroo.py)", "normal", "#555", 11.5)], anchor="end")
file_icon(XR, Y_FQ, "fastq", RNA)
text(XR + 40, Y_FQ - 4, [("RNA-seq reads", "bold"),
                        ("single- or paired-end; one or", "normal", "#555", 12.5),
                        ("more samples per isolate", "normal", "#555", 12.5)])
file_icon(XR, Y_FA, "tsv", RNA)
text(XR + 40, Y_FA - 4, [("Sample metadata", "bold"),
                        ("links each sample to its isolate", "normal", "#555", 12.5)])

# ---- PanScript steps ---------------------------------------------------------
# input checks: an interchange across both tracks
add(f'<rect x="{XL - 26}" y="{Y_CHK - 20}" width="{XR - XL + 52}" height="40" rx="20" '
    f'fill="#fff" stroke="{MRG}" stroke-width="4"/>')
text(XM, Y_CHK + 5, [("Input checks", "bold", MRG, 14)], anchor="middle")
text(XR + 40, Y_CHK - 14, [("metadata (single/paired, names)", "normal", "#555", 12),
                          ("Panaroo loci present in each FASTA", "normal", "#555", 12),
                          ("strandedness from a read subsample", "normal", "#555", 12)])

station(XL, Y_IDX, GEN)
text(XL - 26, Y_IDX - 2, [("Index", "bold"),
                         ("kallisto, one per isolate", "normal", "#555", 12.5)], anchor="end")
station(XR, Y_IDX, RNA)
text(XR + 26, Y_IDX - 2, [("Trimming and read QC", "bold"),
                         ("fastp", "normal", "#555", 12.5)])

def output(y, name, colour=MRG, w=205):
    """File chip left of a station: an output (dark) or an extra input (teal)."""
    add(f'<rect x="{XM - 36 - w}" y="{y - 13}" width="{w}" height="26" rx="5" fill="#fff" '
        f'stroke="{colour}" stroke-width="1.5"/>')
    text(XM - 36 - w / 2, y + 4.5, [(name, "normal", colour, 12)], anchor="middle", mono=True)
    line(f"M{XM - 36},{y} H{XM - 12}", colour, width=2)


STEP_X = XM + 34
station(XM, Y_QNT, MRG, r=13)
text(STEP_X, Y_QNT - 2, [("Pseudo-alignment", "bold"),
                         ("kallisto: each sample against its own", "normal", "#555", 12.5),
                         ("isolate's genes (no single reference)", "normal", "#555", 12.5)])

station(XM, Y_MRG, MRG)
text(STEP_X, Y_MRG - 2, [("Merge into gene groups", "bold"),
                         ("loci → Panaroo groups; split genes summed;", "normal", "#555", 12.5),
                         ("gene absent from an isolate = NA, not 0", "normal", "#555", 12.5)])
output(Y_MRG, "gene_presence_absence", GEN, 250)


station(XM, Y_LEN, MRG)
text(STEP_X, Y_LEN - 2, [("Ortholog length scaling", "bold"),
                         ("to the gene group's median length", "normal", "#555", 12.5)])
output(Y_LEN, "raw_counts.tsv")

station(XM, Y_SF, MRG)
text(STEP_X, Y_SF - 2, [("Library size: core-genome size factors", "bold"),
                        ("DESeq2 median-of-ratios on core genes", "normal", "#555", 12.5)])
output(Y_SF, "norm_counts.tsv")

station(XM, Y_NTPM, MRG)
text(STEP_X, Y_NTPM - 2, [("nTPM", "bold"),
                          ("per own gene length, core size factors;", "normal", "#555", 12.5),
                          ("comparable between genes and samples", "normal", "#555", 12.5)])
output(Y_NTPM, "ntpm.tsv")

station(XM, Y_END, MRG)
text(STEP_X, Y_END - 2, [("Reports", "bold"),
                         ("sample UMAP · MultiQC · software versions", "normal", "#555", 12.5)])

# ---- downstream ---------------------------------------------------------------
add(f'<circle cx="{XM}" cy="{Y_DOWN}" r="13" fill="#fff" stroke="{GREY}" stroke-width="3" '
    f'stroke-dasharray="4 4"/>')
text(STEP_X, Y_DOWN - 10, [("Downstream analysis", "bold", "#555"),
                           ("differential expression, core and accessory", "normal", GREY, 12.5),
                           ("gene expression, genotype–expression links", "normal", GREY, 12.5)])

add("</svg>")
Path(__file__).with_name("panscript_pipeline.svg").write_text("\n".join(out) + "\n")
