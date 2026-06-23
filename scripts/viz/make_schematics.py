#!/usr/bin/env python3
"""Schematic figures for the Plastanno v2 manuscript (matplotlib box-and-arrow).
Fig1  the 7-step hybrid pipeline (dual engines converging on reconciliation)
Fig2  the reconciliation layer in detail
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
LIGHT = "#f0f4f8"

def box(ax, x, y, w, h, text, fc=LIGHT, ec=GREY, fs=9, bold=False, tc="black", lw=1.2):
    p = FancyBboxPatch((x-w/2, y-h/2), w, h, boxstyle="round,pad=0.02,rounding_size=0.08",
                       fc=fc, ec=ec, lw=lw, zorder=2)
    ax.add_patch(p)
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, zorder=3,
            color=tc, fontweight="bold" if bold else "normal")

def arrow(ax, x1, y1, x2, y2, color="black", lw=1.4, style="-|>"):
    # zorder above the boxes (2) and their text (3) so arrowheads are never
    # covered by a box drawn on top of them.
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style,
                 mutation_scale=14, color=color, lw=lw, zorder=5,
                 shrinkA=2, shrinkB=2))

def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.pdf", bbox_inches="tight")
    plt.close(fig); print(f"  wrote {name}.png / .pdf")

# ============================ Fig 1 ============================
def fig1():
    fig, ax = plt.subplots(figsize=(8.4, 8.6)); ax.set_xlim(0, 12); ax.set_ylim(0, 12); ax.axis("off")
    cx = 6.0
    box(ax, cx, 11.3, 7.6, 0.7, "1  Read FASTA  (single plastome record)", fc="white", bold=True)
    box(ax, cx, 10.1, 7.6, 0.8, "2  IR detection  (self-BLASTN)\nlongest minus-strand HSP ≥ 10 kb", fc=LIGHT)
    ax.text(cx+4.05, 10.1, "{LSC, IRb,\nSSC, IRa}", ha="left", va="center", fontsize=7.5, style="italic", color=GREY)
    box(ax, cx, 8.9, 7.6, 0.8, "3  Closest relatives  (BLAST genus_reps)\ntaxonomy → tRNA-DB tier", fc=LIGHT)
    # split to two engines
    ax_left, ax_right = 3.1, 8.9
    box(ax, ax_left, 7.0, 5.2, 1.5,
        "4A  Engine A — reference\n" + r"$\bf{Exonerate}$ protein2genome" + "\nper gene, both IR copies\nrRNA via BLAST  →  $s_{ref}$",
        fc="#deebf7", ec=ENGA, fs=8.5)
    box(ax, ax_right, 7.0, 5.2, 1.5,
        "4B  Engine B — model\n6-frame → " + r"$\bf{hmmsearch}$ (CDS)" + "\nARAGORN + BLAST (tRNA)\nBLAST (rRNA)  →  $s_{model}$",
        fc="#e2f3e6", ec=ENGB, fs=8.5)
    box(ax, cx, 4.9, 8.0, 1.0,
        "5  Reconciliation\nbin → confidence → ORF validation → locus selection",
        fc="#fde6d6", ec=ORANGE, fs=8.8, bold=True)
    box(ax, cx, 3.4, 8.0, 0.95,
        "6  Special cases\nCAU naming · rps12 trans-splice · short exons\nsplice refinement · internal-stop QC", fc=LIGHT, fs=8.5)
    box(ax, cx, 1.9, 8.0, 0.95,
        "7  Output\n.gb  .gff3  .faa  .ffn  .frn  .report\n+ circular plastome map  (.png / .pdf / .svg)", fc="white", bold=True, fs=8.6)
    # arrows
    arrow(ax, cx, 10.95, cx, 10.5)
    arrow(ax, cx, 9.7, cx, 9.3)
    arrow(ax, cx, 8.5, ax_left, 7.75); arrow(ax, cx, 8.5, ax_right, 7.75)
    arrow(ax, ax_left, 6.25, cx-0.7, 5.4, color=ENGA); arrow(ax, ax_right, 6.25, cx+0.7, 5.4, color=ENGB)
    arrow(ax, cx, 4.4, cx, 3.88, color=ORANGE)
    arrow(ax, cx, 2.92, cx, 2.3)
    ax.text(0.2, 11.7, "Fig. 1", fontsize=11, fontweight="bold")
    ax.text(cx, 0.85, "Two independent engines (reference-based + model-based) are\nmerged by a scoring/reconciliation layer — the defining hybrid design.",
            ha="center", va="center", fontsize=8, style="italic", color=GREY)
    save(fig, "Fig1_pipeline")

# ============================ Fig 2 ============================
def fig2():
    fig, ax = plt.subplots(figsize=(9.6, 6.6)); ax.set_xlim(0, 13.5); ax.set_ylim(0, 10); ax.axis("off")
    # inputs
    box(ax, 2.3, 8.8, 4.0, 0.95, "Engine A features\n$s_{ref}$ (Exonerate)", fc="#deebf7", ec=ENGA, fs=8.6)
    box(ax, 2.3, 7.3, 4.0, 0.95, "Engine B features\n$s_{model}$ (HMM/ARAGORN/BLAST)", fc="#e2f3e6", ec=ENGB, fs=8.6)
    # binning
    box(ax, 7.0, 8.05, 4.2, 1.15, "Binning  (match_features)\nJaccard overlap by gene name", fc="#fde6d6", ec=ORANGE, fs=8.7, bold=True)
    arrow(ax, 4.35, 8.6, 4.85, 8.3, color=ENGA); arrow(ax, 4.35, 7.5, 4.85, 7.85, color=ENGB)
    # bins (single line each, to the right of the binning box)
    bins = [("AB", "overlap ≥ 0.5, both engines"), ("A_only", "Engine A only"),
            ("B_only", "Engine B only"), ("conflict", "same name, low overlap → IR copies")]
    for i, (b, d) in enumerate(bins):
        y = 9.25 - i*0.6
        ax.text(9.5, y, f"• {b}", fontsize=8, va="center", fontweight="bold", color="#555555")
        ax.text(10.7, y, d, fontsize=7.2, va="center", color="#777777")
    # pipeline down
    steps = [
        ("Adaptive confidence  (compute_confidence)", "weighted mean over the signals actually available,\nrenormalised to 1  ·  {overlap, ref, model, orf}", "#fff4ec"),
        ("ORF validation  (validate_orf)", "start/stop codon · internal stops · length-vs-expected,\non the spliced coding sequence", LIGHT),
        ("Locus selection  (_select)", "union-find clusters (duplicate fragments / paralog cross-hits);\nkeep best per cluster; drop CDS < 0.6× expected length", LIGHT),
    ]
    y = 6.0
    for title, desc, fc in steps:
        p = FancyBboxPatch((7.0-4.5, y-0.525), 9.0, 1.05, boxstyle="round,pad=0.02,rounding_size=0.08",
                           fc=fc, ec=GREY, lw=1.2, zorder=2)
        ax.add_patch(p)
        ax.text(7.0, y+0.27, title, ha="center", va="center", fontsize=8.6, fontweight="bold", zorder=3)
        ax.text(7.0, y-0.18, desc, ha="center", va="center", fontsize=7.8, color="#444444", zorder=3)
        y -= 1.5
    arrow(ax, 7.0, 7.5, 7.0, 6.55, color=ORANGE)
    arrow(ax, 7.0, 5.48, 7.0, 5.05)
    arrow(ax, 7.0, 3.98, 7.0, 3.55)
    # output with flags
    box(ax, 7.0, 1.4, 9.0, 1.2,
        "Unified annotation set\nper-feature flag:  C ≥ 0.8 HIGH   ·   ≥ 0.5 MEDIUM   ·   else NEEDS_REVIEW\nfull provenance (engine, component scores) → GenBank /note + GFF3",
        fc="white", ec="black", fs=8.4, bold=False)
    arrow(ax, 7.0, 2.48, 7.0, 2.05)
    ax.text(0.2, 9.6, "Fig. 2", fontsize=11, fontweight="bold")
    save(fig, "Fig2_reconciliation")

# ============================ Fig 7 ============================
def fig7():
    fig, ax = plt.subplots(figsize=(9.6, 7.4)); ax.set_xlim(0, 14); ax.set_ylim(0, 11); ax.axis("off")
    # source
    box(ax, 7.0, 10.3, 6.8, 0.7, "Public RefSeq plastomes", fc="white", bold=True, fs=9.6)
    # two branches: development collection (builds DBs) and an independent test set
    box(ax, 3.4, 8.55, 4.8, 1.05,
        "Development collection\n8,807 plastomes\n(2,899 genus-representatives)", fc="#deebf7", ec=ENGA, fs=8.2)
    pT = FancyBboxPatch((10.6 - 2.6, 8.45 - 0.625), 5.2, 1.25, boxstyle="round,pad=0.03,rounding_size=0.08",
                        fc="#fdd0a2", ec=ORANGE, lw=1.4, zorder=2); ax.add_patch(pT)
    ax.text(10.6, 8.72, "Independent test set  (n = 2,151)", ha="center", fontweight="bold", fontsize=8.8, zorder=3)
    ax.text(10.6, 8.18, "not represented in any reference database", ha="center",
            fontsize=7.0, style="italic", color="#a63603", zorder=3)
    ax.plot([7.0, 7.0], [9.95, 9.6], color="black", lw=1.4, zorder=5)
    arrow(ax, 7.0, 9.6, 3.4, 9.03, color=ENGA)
    arrow(ax, 7.0, 9.6, 10.6, 9.08, color=ORANGE)
    # ---- left: build → databases ----
    box(ax, 3.4, 7.1, 4.8, 0.85, "build_all.py · build_protein_db.py\n(genes per gene_catalog)", fc=LIGHT, fs=8.3)
    arrow(ax, 3.4, 8.0, 3.4, 7.55, color=ENGA)
    dbs = ("• genus-representative BLAST (blast_db)\n"
           "• per-gene proteins (protein_db)\n"
           "• profile HMMs (hmm_db)\n"
           "• tiered tRNA (trna_db) + intron exons (exon_db)\n"
           "• full-length rRNA (rrna_db)\n"
           "• exon panel + length templates (boundary_db)\n"
           "• gene_catalog.json")
    p = FancyBboxPatch((3.4-2.7, 4.35-1.45), 5.4, 2.9, boxstyle="round,pad=0.03,rounding_size=0.08",
                       fc="#eef6ff", ec=ENGA, lw=1.4, zorder=2); ax.add_patch(p)
    ax.text(3.4, 5.5, "Reference databases", ha="center", fontweight="bold", fontsize=8.8, zorder=3)
    ax.text(3.4, 4.5, dbs, ha="center", va="center", fontsize=7.6, color="#333333", zorder=3)
    arrow(ax, 3.4, 6.67, 3.4, 5.85, color=ENGA)
    # ---- right: test-set composition (two evaluation strata, one box) ----
    pC = FancyBboxPatch((10.6 - 2.6, 6.15 - 0.72), 5.2, 1.44, boxstyle="round,pad=0.03,rounding_size=0.08",
                        fc=LIGHT, ec=GREY, lw=1.2, zorder=2); ax.add_patch(pC)
    ax.text(10.6, 6.55, "Comprising two strata", ha="center", fontweight="bold", fontsize=8.2, zorder=3)
    ax.text(10.6, 6.10, "• Held-out:  498", ha="center", fontsize=8.0, zorder=3)
    ax.text(10.6, 5.68, "• Prospective (published later):  1,653", ha="center", fontsize=7.8, zorder=3)
    arrow(ax, 10.6, 7.82, 10.6, 6.92, color=ORANGE)
    # ---- converge: evaluation ----
    box(ax, 7.0, 1.7, 9.6, 1.05,
        "Leakage-free evaluation\n"
        "Global F1 = 92.5%   (held-out 90.5 · prospective 93.1)",
        fc="white", ec="black", fs=9.0, bold=False)
    arrow(ax, 3.4, 2.85, 5.6, 2.3, color=ENGA)
    arrow(ax, 10.6, 5.4, 8.4, 2.3, color=ORANGE)
    ax.text(0.2, 10.7, "Fig. 7", fontsize=11, fontweight="bold")
    ax.text(7.0, 0.75, "Reference databases are built from a development collection; the tool is evaluated on an\n"
            "independent set of plastomes not represented in any database — measuring true generalisation.",
            ha="center", va="center", fontsize=8, style="italic", color=GREY)
    save(fig, "Fig7_database_design")

if __name__ == "__main__":
    print("Fig1 …"); fig1()
    print("Fig2 …"); fig2()
    print("Fig7 …"); fig7()
    print("done.")
