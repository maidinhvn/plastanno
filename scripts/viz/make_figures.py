#!/usr/bin/env python3
"""Publication figures for the Plastanno v2 manuscript.

Fig3  held-out performance (by type, per-genome distribution, by stratum)
Fig4  head-to-head vs PGA (paired distribution + scatter)
Reads the committed benchmark artefacts; writes 300-dpi PNG + vector PDF to docs/figures/.
"""
import json, glob, csv, os
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from collections import Counter

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(REPO, "docs", "figures"); os.makedirs(OUT, exist_ok=True)

# ---- style ----
plt.rcParams.update({
    "font.size": 9, "font.family": "DejaVu Sans", "axes.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8, "figure.dpi": 150,
})
TEAL, ORANGE, GREY = "#2c7fb8", "#e6550d", "#969696"
TYPE_C = {"CDS": "#1b9e77", "tRNA": "#7570b3", "rRNA": "#d95f02"}

def f1(tp, ref, pred):
    sn = tp/ref if ref else 0; pr = tp/pred if pred else 0
    return (200*sn*pr/(sn+pr) if sn+pr else 0), 100*sn, 100*pr

def load_heldout():
    strata = dict(l.split() for l in open(f"{REPO}/splits/heldout_v2_strata.tsv"))
    recs = []
    for cp in glob.glob(f"{REPO}/bench_runs/reconfirm/chunk_0[012]/results.json") + \
              glob.glob(f"{REPO}/bench_runs/heldout_v2_final/chunk_*/results.json"):
        for r in json.load(open(cp))["results"]:
            if r.get("status") == "OK":
                r["stratum"] = strata.get(r["acc"], "?"); recs.append(r)
    return recs

def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.pdf", bbox_inches="tight")
    plt.close(fig); print(f"  wrote {name}.png / .pdf")

# ============================ Fig 3 ============================
def fig3(recs):
    fig, ax = plt.subplots(1, 3, figsize=(10, 3.2))
    # (A) Sn/Pr/F1 by type
    types = ["CDS", "tRNA", "rRNA"]; metrics = ["Sn", "Pr", "F1"]
    vals = {}
    for t in types:
        tp = sum(r["tp_by_type"].get(t, 0) for r in recs)
        rf = sum(r["ref_by_type"].get(t, 0) for r in recs)
        pd = sum(r["pred_by_type"].get(t, 0) for r in recs)
        f, sn, pr = f1(tp, rf, pd); vals[t] = [sn, pr, f]
    x = np.arange(len(types)); w = 0.26
    for i, m in enumerate(metrics):
        ax[0].bar(x + (i-1)*w, [vals[t][i] for t in types], w, label=m,
                  color=[TEAL, ORANGE, "#31a354"][i], edgecolor="white", linewidth=0.5)
    for i, t in enumerate(types):
        ax[0].text(x[i]+w, vals[t][2]+0.6, f"{vals[t][2]:.1f}", ha="center", fontsize=7.5, fontweight="bold")
    ax[0].set_xticks(x); ax[0].set_xticklabels(types); ax[0].set_ylim(70, 101)
    ax[0].set_ylabel("Score (%)")
    ax[0].legend(frameon=False, fontsize=7.5, ncol=3, loc="lower center",
                 bbox_to_anchor=(0.5, -0.32))
    ax[0].set_title("(A) Held-out accuracy by feature type", fontsize=9, loc="left")
    # (B) per-genome global F1 distribution
    pg = [f1(r["tp"], r["ref"], r["pred"])[0] for r in recs]
    ax[1].hist(pg, bins=30, color=TEAL, edgecolor="white", linewidth=0.4)
    med = np.median(pg)
    ax[1].axvline(med, color=ORANGE, lw=1.5, ls="--", label=f"median {med:.1f}")
    ax[1].set_xlabel("Per-genome global F1 (%)"); ax[1].set_ylabel("Genomes")
    ax[1].legend(frameon=False, fontsize=7.5)
    ax[1].set_title(f"(B) Per-genome F1 (n={len(recs)})", fontsize=9, loc="left")
    # (C) by stratum
    strata = [("heldout508", "Held-out\n(unseen split)"), ("inscope", "Prospective\n(new RefSeq)")]
    gv = []
    allf = f1(sum(r["tp"] for r in recs), sum(r["ref"] for r in recs), sum(r["pred"] for r in recs))[0]
    labs = ["All\n(n=%d)" % len(recs)]; gv = [allf]; cols = [GREY]
    for key, lab in strata:
        sub = [r for r in recs if r["stratum"] == key]
        gv.append(f1(sum(r["tp"] for r in sub), sum(r["ref"] for r in sub), sum(r["pred"] for r in sub))[0])
        labs.append(f"{lab}\n(n={len(sub)})"); cols.append(TEAL)
    xb = np.arange(len(gv))
    ax[2].bar(xb, gv, 0.55, color=cols, edgecolor="white")
    for i, v in enumerate(gv): ax[2].text(i, v+0.4, f"{v:.1f}", ha="center", fontsize=8, fontweight="bold")
    ax[2].set_xticks(xb); ax[2].set_xticklabels(labs, fontsize=7.5); ax[2].set_ylim(80, 96)
    ax[2].set_ylabel("Global F1 (%)")
    ax[2].set_title("(C) Generalisation by stratum", fontsize=9, loc="left")
    fig.tight_layout(); save(fig, "Fig3_performance")

# ============================ Fig 4 ============================
def fig4():
    rows = list(csv.DictReader(open(f"{REPO}/docs/Plastanno_GD4_h2h_paired.csv")))
    pl = np.array([float(r["plastanno_f1"]) for r in rows])
    pg = np.array([float(r["pga_f1"]) for r in rows])
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.4))
    # (A) paired violin + box
    parts = ax[0].violinplot([pl, pg], positions=[1, 2], showmeans=False, showextrema=False, widths=0.8)
    for b, c in zip(parts["bodies"], [TEAL, ORANGE]):
        b.set_facecolor(c); b.set_alpha(0.35); b.set_edgecolor(c)
    bp = ax[0].boxplot([pl, pg], positions=[1, 2], widths=0.25, patch_artist=True,
                       showfliers=False, medianprops=dict(color="black", lw=1.2))
    for patch, c in zip(bp["boxes"], [TEAL, ORANGE]): patch.set_facecolor(c); patch.set_alpha(0.8)
    ax[0].set_xticks([1, 2]); ax[0].set_xticklabels([f"Plastanno\n(med {np.median(pl):.1f})",
                                                     f"PGA\n(med {np.median(pg):.1f})"])
    ax[0].set_ylabel("Per-genome global F1 (%)")
    ax[0].set_title(f"(A) Paired F1 (n={len(rows)})", fontsize=9, loc="left")
    # (B) scatter
    ax[1].scatter(pg, pl, s=10, color=TEAL, alpha=0.5, edgecolor="none")
    lim = [min(pg.min(), pl.min())-2, 100]
    ax[1].plot(lim, lim, color=GREY, lw=1, ls="--")
    wins = int((pl > pg+0.5).sum()); ties = int((abs(pl-pg) <= 0.5).sum())
    ax[1].text(0.05, 0.95, f"Plastanno better: {wins}\nties: {ties}\nPGA better: {len(rows)-wins-ties}",
               transform=ax[1].transAxes, va="top", fontsize=7.5,
               bbox=dict(boxstyle="round", fc="white", ec=GREY, alpha=0.9))
    ax[1].set_xlim(lim); ax[1].set_ylim(lim); ax[1].set_aspect("equal")
    ax[1].set_xlabel("PGA F1 (%)"); ax[1].set_ylabel("Plastanno F1 (%)")
    ax[1].set_title("(B) Per-genome comparison", fontsize=9, loc="left")
    fig.tight_layout(); save(fig, "Fig4_headtohead")

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
    print("Fig3 …"); fig3(load_heldout())
    print("Fig4 …"); fig4()
    print("Fig5 …"); fig5()
    print("done.")
