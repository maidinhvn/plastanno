#!/usr/bin/env python3
"""Fig5 — noise in public reference annotations (the problem Plastanno targets).

This file used to also build Fig3 (held-out performance) and Fig4 (head-to-head
vs PGA). Both are gone from the public repository, for two separate reasons:

  * they read `bench_runs/`, which is not published here, so nobody outside
    could reproduce them;
  * Fig3's "leakage-free held-out set" was not leakage-free. The reference
    databases were built before the split, so most of that set was represented
    in them and the figure's headline F1 was measured on a contaminated
    denominator. Fig4 compared against one other tool with a scorer that has
    since been replaced.

Replacement evaluation figures will be published with the paper, once the
scorer and protocol have had an independent review. Nothing measured is being
withheld here; what was wrong has been withdrawn rather than left standing.

Fig5 is unaffected: it counts properties of public reference annotations
(committed in fig5_data.json), not Plastanno's performance.

Writes 300-dpi PNG + vector PDF to docs/figures/."""
import json, os
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(REPO, "docs", "figures"); os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "font.size": 9, "font.family": "DejaVu Sans", "axes.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8, "figure.dpi": 150,
})
TEAL, ORANGE, GREY = "#2c7fb8", "#e6550d", "#969696"


DPI = 600      # the resolution the committed PNGs are at


def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=DPI, bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.pdf", bbox_inches="tight")
    plt.close(fig); print(f"  wrote {name}.png / .pdf")


# ============================ Fig 5 ============================
def fig5():
    d = json.load(open(f"{OUT}/fig5_data.json"))
    fig, ax = plt.subplots(1, 2, figsize=(7.6, 3.4))
    # (A) IR annotation gap (donut)
    sizes = [d["ir_annotated"], d["ir_absent"]]
    cols = [TEAL, ORANGE]
    wedges, _ = ax[0].pie(sizes, colors=cols, startangle=90,
                          wedgeprops=dict(width=0.42, edgecolor="white", linewidth=1.2))
    ax[0].text(0, 0.12, f"{d['ir_absent_pct']:.0f}%", ha="center", fontsize=20, fontweight="bold", color=ORANGE)
    ax[0].text(0, -0.18, "no IR\nannotation", ha="center", fontsize=8, color=ORANGE)
    ax[0].legend(wedges, [f"IR annotated ({d['ir_annotated']})",
                          f"IR absent ({d['ir_absent']})"],
                 loc="lower center", bbox_to_anchor=(0.5, -0.18), frameon=False, fontsize=7.5)
    ax[0].set_title(f"(A) IR annotation in references (n={d['n_refs']})", fontsize=9, loc="left")
    # (B) naming heterogeneity: distinct spellings vs canonical gene count
    types = ["rRNA", "tRNA", "CDS"]
    canon = {"rRNA": 4, "tRNA": 30, "CDS": 80}  # approx canonical gene counts in a plastome
    spell = [d["rrna_total_spellings"], d["trna_total_spellings"], d["cds_total_spellings"]]
    x = np.arange(len(types)); w = 0.38
    ax[1].bar(x - w/2, [canon[t] for t in types], w, label="Canonical genes", color=GREY, edgecolor="white")
    ax[1].bar(x + w/2, spell, w, label="Distinct name spellings\nin references", color=ORANGE, edgecolor="white")
    for i, t in enumerate(types):
        ax[1].text(i + w/2, spell[i] + 8, f"{spell[i]}", ha="center", fontsize=8, fontweight="bold")
    ax[1].set_xticks(x); ax[1].set_xticklabels(types); ax[1].set_ylabel("Count")
    ax[1].legend(frameon=False, fontsize=7.5, loc="upper left")
    ax[1].set_title("(B) Gene-name inconsistency across references", fontsize=9, loc="left")
    fig.tight_layout(); save(fig, "Fig5_reference_noise")

if __name__ == "__main__":
    print("Fig5 …"); fig5()
    print("done.")
