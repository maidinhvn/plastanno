#!/usr/bin/env python3
"""Schematic figures for Plastanno (matplotlib box-and-arrow).

Fig1  the pipeline as 3.0.0 actually runs it
Fig2  the candidate-selection layer

Both are drawn from the code, not from the prose: step numbering follows the
`Step N` markers in plastanno/pipeline.py, and Fig2 follows the docstring of
core.reconcile.reconcile_pooled. The previous versions of both drew the
pre-3.0.0 architecture -- a cross-engine reconciliation layer whose claimed
benefit the H3 ablation did not find -- so they described a mechanism that is
neither the default nor a claim we still make.

Writes 300-dpi PNG + vector PDF to docs/figures/."""
import os
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(REPO, "docs", "figures"); os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "figure.dpi": 150})

TEAL, ORANGE, GREY = "#2c7fb8", "#e6550d", "#969696"
ENGA, ENGB = "#3690c0", "#41ab5d"        # engine A (ref) / B (model)
PURPLE = "#7b6bb0"                        # boundary refinement
LIGHT = "#f0f4f8"


def box(ax, x, y, w, h, text, fc=LIGHT, ec=GREY, fs=9, bold=False, tc="black", lw=1.2):
    p = FancyBboxPatch((x-w/2, y-h/2), w, h, boxstyle="round,pad=0.02,rounding_size=0.08",
                       fc=fc, ec=ec, lw=lw, zorder=2)
    ax.add_patch(p)
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, zorder=3,
            color=tc, fontweight="bold" if bold else "normal")


def arrow(ax, x1, y1, x2, y2, color="black", lw=1.4, style="-|>"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style,
                 mutation_scale=14, color=color, lw=lw, zorder=5,
                 shrinkA=2, shrinkB=2))


DPI = 600      # the resolution the committed PNGs are at


def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=DPI, bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.pdf", bbox_inches="tight")
    plt.close(fig); print(f"  wrote {name}.png / .pdf")


# ============================ Fig 1 ============================
def fig1():
    fig, ax = plt.subplots(figsize=(8.8, 9.8))
    ax.set_xlim(0, 12); ax.set_ylim(0, 13.4); ax.axis("off")
    cx = 6.0

    box(ax, cx, 12.7, 7.8, 0.68, "1  Read FASTA  (single plastome record)",
        fc="white", bold=True, fs=9)
    box(ax, cx, 11.55, 7.8, 0.8,
        "2  IR detection  (self-BLASTN)\nlongest minus-strand HSP ≥ 10 kb", fc=LIGHT, fs=8.8)
    ax.text(cx+4.15, 11.55, "{LSC, IRb,\nSSC, IRa}", ha="left", va="center",
            fontsize=7.5, style="italic", color=GREY)
    box(ax, cx, 10.4, 7.8, 0.8,
        "3  Closest relatives  (BLAST genus_reps)\nranked neighbours for the tRNA search", fc=LIGHT, fs=8.8)

    ax_left, ax_right = 3.05, 8.95
    box(ax, ax_left, 8.45, 5.3, 1.55,
        "4A  Engine A — reference\n" + r"$\bf{Exonerate}$ protein2genome" +
        "\nper gene, both IR copies\nrRNA via BLAST  →  $s_{ref}$",
        fc="#deebf7", ec=ENGA, fs=8.5)
    box(ax, ax_right, 8.45, 5.3, 1.55,
        "4B  Engine B — model\n6-frame → " + r"$\bf{hmmsearch}$ (CDS)" +
        "\nARAGORN + tRNAscan-SE + BLAST (tRNA)\nBLAST (rRNA)  →  $s_{model}$",
        fc="#e2f3e6", ec=ENGB, fs=8.5)

    box(ax, cx, 6.25, 8.4, 1.15,
        "5  Candidate selection   (--mode pooled)\n"
        "pool both engines → B-primary / A-rescue → ORF validation → one per locus\n"
        "engine agreement earns no score bonus",
        fc="#fde6d6", ec=ORANGE, fs=8.4, bold=False)

    box(ax, cx, 4.75, 8.4, 0.95,
        "6  Special cases\nCAU naming · rps12 trans-splice · short exons · internal-stop QC",
        fc=LIGHT, fs=8.4)

    box(ax, cx, 3.2, 8.4, 1.2,
        "6b  Boundary refinement   (inventory already fixed)\n"
        "multi-exon CDS splice sites · intron-free tRNA (tRNAscan-SE)\n"
        "intron-bearing tRNA: glocal exon placement (7 genes)",
        fc="#ece7f6", ec=PURPLE, fs=8.4)

    box(ax, cx, 1.5, 8.4, 1.15,
        "7  Output\n.gb  .gff3  .tbl  .faa  .ffn  .frn  .report\n"
        ".provenance.json  .trna_alternatives.tsv  + circular map",
        fc="white", bold=True, fs=8.5)

    arrow(ax, cx, 12.36, cx, 11.95)
    arrow(ax, cx, 11.15, cx, 10.8)
    arrow(ax, cx, 10.0, ax_left, 9.23); arrow(ax, cx, 10.0, ax_right, 9.23)
    arrow(ax, ax_left, 7.67, cx-0.8, 6.83, color=ENGA)
    arrow(ax, ax_right, 7.67, cx+0.8, 6.83, color=ENGB)
    arrow(ax, cx, 5.67, cx, 5.23, color=ORANGE)
    arrow(ax, cx, 4.27, cx, 3.8)
    arrow(ax, cx, 2.6, cx, 2.08, color=PURPLE)

    ax.text(0.2, 13.1, "Fig. 1", fontsize=11, fontweight="bold")
    ax.text(cx, 0.45,
            "Two engines supply candidates; a single selector keeps one feature per locus.\n"
            "Boundaries are refined only after the inventory is fixed, so a boundary change\n"
            "cannot add or remove a gene.",
            ha="center", va="center", fontsize=8, style="italic", color=GREY)
    save(fig, "Fig1_pipeline")


# ============================ Fig 2 ============================
def fig2():
    fig, ax = plt.subplots(figsize=(9.8, 7.0))
    ax.set_xlim(0, 13.5); ax.set_ylim(0, 10.4); ax.axis("off")

    box(ax, 2.35, 9.1, 4.1, 0.95, "Engine A candidates\n$s_{ref}$ (Exonerate)",
        fc="#deebf7", ec=ENGA, fs=8.6)
    box(ax, 2.35, 7.6, 4.1, 0.95,
        "Engine B candidates\n$s_{model}$  (HMM · ARAGORN\ntRNAscan-SE · BLAST)",
        fc="#e2f3e6", ec=ENGB, fs=8.2)

    box(ax, 7.2, 8.35, 4.3, 1.1, "Pool by gene name\n(pooled_candidates)",
        fc="#fde6d6", ec=ORANGE, fs=8.7, bold=True)
    arrow(ax, 4.45, 8.9, 5.0, 8.6, color=ENGA)
    arrow(ax, 4.45, 7.8, 5.0, 8.1, color=ENGB)

    rules = [
        ("both engines",
         ["keep B's coordinates when its ORF check passes, or is",
          "at least as good as A's; substitute A only otherwise"]),
        ("one engine",
         ["keep it — this is the rescue that lifts inventory"]),
        ("same name, far apart",
         ["an IR-duplicated gene → keep both copies;",
          "otherwise ORF validity decides"]),
    ]
    y = 9.95
    for lab, lines in rules:
        ax.text(9.62, y, "• " + lab, fontsize=8.0, va="center",
                fontweight="bold", color="#555555")
        y -= 0.42
        for ln in lines:
            ax.text(9.89, y, ln, fontsize=7.0, va="center", color="#777777")
            y -= 0.40
        y -= 0.12

    steps = [
        ("ORF validation   (validate_orf)",
         "start/stop codon · internal stops · length-vs-expected,\non the spliced coding sequence", LIGHT),
        ("Selection score   (compute_confidence)",
         "single-engine signals only: {ref, orf} or {model, orf}.\n"
         "No cross-engine agreement term — H3 found it bought nothing.", "#fff4ec"),
        ("Locus selection   (_select)",
         "union-find clusters (duplicate fragments / paralog cross-hits);\n"
         "keep best per cluster; drop CDS < 0.6× expected length", LIGHT),
    ]
    y = 6.15
    for title, desc, fc in steps:
        p = FancyBboxPatch((7.2-4.6, y-0.545), 9.2, 1.09,
                           boxstyle="round,pad=0.02,rounding_size=0.08",
                           fc=fc, ec=GREY, lw=1.2, zorder=2)
        ax.add_patch(p)
        ax.text(7.2, y+0.28, title, ha="center", va="center", fontsize=8.6,
                fontweight="bold", zorder=3)
        ax.text(7.2, y-0.19, desc, ha="center", va="center", fontsize=7.6,
                color="#444444", zorder=3)
        y -= 1.52

    arrow(ax, 7.2, 7.8, 7.2, 6.72, color=ORANGE)
    arrow(ax, 7.2, 5.6, 7.2, 5.18)
    arrow(ax, 7.2, 4.08, 7.2, 3.66)

    box(ax, 7.2, 1.35, 9.2, 1.15,
        "One feature per locus\n"
        "flag:  ≥ 0.8 HIGH   ·   ≥ 0.5 MEDIUM   ·   else NEEDS_REVIEW\n"
        "provenance (engine, component scores) → GenBank /note + GFF3",
        fc="white", ec="black", fs=8.3)
    arrow(ax, 7.2, 2.55, 7.2, 1.98)

    ax.text(0.2, 10.0, "Fig. 2", fontsize=11, fontweight="bold")
    ax.text(7.2, 0.25,
            "The score ranks candidates for one locus. It is deliberately not called a confidence.",
            ha="center", va="center", fontsize=7.6, style="italic", color=GREY)
    save(fig, "Fig2_selection")


if __name__ == "__main__":
    print("Fig1 …"); fig1()
    print("Fig2 …"); fig2()
    print("done.")
