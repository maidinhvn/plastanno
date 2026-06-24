#!/usr/bin/env python3
"""
plastome_circular_map.py  —  Standalone circular plastome (chloroplast genome) map.

Draws an OGDRAW-style circular map from a single GenBank file:
  * genes as block arrows showing transcription direction
      - (+) strand genes on the OUTER ring, (-) strand on the INNER ring
  * genes coloured by functional class
  * a thin MIDDLE ring showing the quadripartite structure (LSC / IRb / SSC / IRa)
    at the REAL inverted-repeat boundary coordinates, with curved region labels
  * radial gene labels with leader lines, spread to avoid overlap (labels ALL genes)

Inverted-repeat detection (for the quadripartite ring), in order of preference:
  1. --ir  jlb,jsb,jsa,jla            (you supply the four junction coords, 1-based)
  2. self-BLAST  (needs `blastn` on PATH)  -> longest minus-strand self-hit = IR
  3. annotated `repeat_region` features in the GenBank file
  4. none found -> genes are drawn without the region ring

Dependencies: biopython, numpy, matplotlib   (+ NCBI blastn for option 2)

Usage:
  python plastome_circular_map.py genome.gb
  python plastome_circular_map.py genome.gb --out mymap --rotate 180
  python plastome_circular_map.py genome.gb --ir 84244,111291,128863,155910

License: free to use/modify (CC0 / public domain).
"""
import argparse, math, os, re, subprocess, sys, tempfile
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Polygon
from Bio import SeqIO

# ---- Journal-quality figure defaults --------------------------------------
# Embed TrueType fonts in PDF/PS (journals reject Type-3); keep SVG text as text
# (editable); use a sans-serif face. Vector PDF/SVG are resolution-independent;
# the PNG is rasterised at --dpi (default 600, line-art standard).
matplotlib.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "axes.linewidth": 0.6,
})


def _slug(s):
    """Filesystem-safe token from a free-text name (organism)."""
    return re.sub(r"[^A-Za-z0-9]+", "_", (s or "")).strip("_")


# ----------------------------------------------------------------------------
# Gene functional classification -> (label, colour).  Edit freely.
# ----------------------------------------------------------------------------
def categorize(gene_name):
    n = gene_name.lower()
    if n.startswith(("psa", "psb", "pet", "atp")) or n == "rbcl":
        return ("Photosynthesis", "#41ab5d")
    if n.startswith("ndh"):
        return ("NADH dehydrogenase", "#a1d99b")
    if n.startswith(("rpl", "rps")):
        return ("Ribosomal proteins", "#fb6a4a")
    if n.startswith("rpo"):
        return ("RNA polymerase", "#807dba")
    if n.startswith("trn"):
        return ("Transfer RNA (tRNA)", "#969696")
    if n.startswith("rrn"):
        return ("Ribosomal RNA (rRNA)", "#3182bd")
    return ("Other genes (ycf, matK, ...)", "#fdae61")


# ----------------------------------------------------------------------------
# Inverted-repeat detection
# ----------------------------------------------------------------------------
def ir_from_blast(seq):
    """Self-BLAST; return (jlb, jsb, jsa, jla, lsc, ir, ssc) 1-based, or None."""
    L = len(seq)
    try:
        with tempfile.TemporaryDirectory() as td:
            fa = os.path.join(td, "s.fa")
            open(fa, "w").write(f">q\n{seq}\n")
            out = subprocess.run(
                ["blastn", "-query", fa, "-subject", fa, "-evalue", "1e-30",
                 "-word_size", "20", "-dust", "no", "-perc_identity", "90",
                 "-outfmt", "6 length qstart qend sstart send sstrand"],
                capture_output=True, text=True, timeout=600).stdout
    except Exception:
        return None
    best = None
    for ln in out.strip().splitlines():
        f = ln.split("\t")
        if len(f) < 6 or f[5] != "minus":
            continue
        length = int(f[0])
        if best is None or length > best[0]:
            best = (length, (int(f[1]), int(f[2])), (int(f[3]), int(f[4])))
    if not best or best[0] < 4000:
        return None
    _, q, s = best
    a = (min(q), max(q)); b = (min(s), max(s))
    first, second = sorted([a, b])
    between = second[0] - first[1] - 1
    outside = L - second[1] + first[0] - 1
    ir = best[0]
    if between <= outside:                 # SSC between the two IR copies
        irb, ira = first, second
        jlb, jsb, jsa, jla = irb[0], irb[1], ira[0], ira[1]
        ssc, lsc = between, outside
    else:                                  # SSC wraps the origin
        ira, irb = first, second
        jsa, jla, jlb, jsb = ira[1], ira[0], irb[1], irb[0]
        ssc, lsc = outside, between
    return (jlb, jsb, jsa, jla, lsc, ir, ssc)


def ir_from_annotation(rec):
    """Use annotated repeat_region features (>=2) as IR, or None."""
    irs = [(int(f.location.start) + 1, int(f.location.end))
           for f in rec.features if f.type == "repeat_region"]
    irs = [x for x in irs if x[1] - x[0] > 4000]
    if len(irs) < 2:
        return None
    irs.sort(); first, second = irs[0], irs[1]; L = len(rec.seq)
    between = second[0] - first[1] - 1
    return (first[0], first[1], second[0], second[1], L - 2*(second[1]-second[0]) - between,
            second[1] - second[0], between)


# ----------------------------------------------------------------------------
# Drawing
# ----------------------------------------------------------------------------
# Radii (outer -> inner): forward genes | structural ring (MIDDLE) | reverse genes
FWD_IN, FWD_OUT = 0.965, 1.03         # (+) strand arrows  (outer), close to ring
REG_R, REG_LW   = 0.93, 15            # quadripartite ring (MIDDLE), thinner than the two gene rings
REV_IN, REV_OUT = 0.83, 0.895         # (-) strand arrows  (inner), close to ring
LAB_OUT = 1.10                         # forward gene labels (outside)
LAB_IN  = 0.66                         # reverse gene labels (inside)


def _spread(angles, dmin):
    """Push angularly-close labels apart (relaxation) so radial text doesn't collide."""
    a = list(angles); n = len(a)
    for _ in range(300):
        moved = False
        for i in range(n - 1):
            gap = a[i + 1] - a[i]
            if gap < dmin:
                sh = (dmin - gap) / 2.0
                a[i] -= sh; a[i + 1] += sh; moved = True
        if not moved:
            break
    return a


def draw_map(rec, ir, organism, rotate_deg, out, dpi=600):
    L = len(rec.seq)
    ROT = math.radians(rotate_deg)
    ang = lambda p: 2 * math.pi * (p % L) / L
    xy = lambda r, a: (r * math.sin(a + ROT), r * math.cos(a + ROT))

    def arc(ax, r, p0, p1, color, lw, z=1):
        if p1 < p0:
            p1 += L
        th = np.linspace(ang(p0), ang(p0) + 2*math.pi*(p1-p0)/L, max(2, int((p1-p0)/40))) + ROT
        ax.plot(r*np.sin(th), r*np.cos(th), color=color, lw=lw, solid_capstyle="butt", zorder=z)

    def gene_arrow(ax, p0, p1, strand, r_in, r_out, color):
        if p1 < p0:
            p1 += L
        span = p1 - p0
        ah = min(span * 0.55, L * 0.0035)
        r_mid = (r_in + r_out) / 2
        ns = max(2, int(span / 25))
        pts = []
        if strand >= 0:
            be = p1 - ah
            pts += [xy(r_out, ang(p)) for p in np.linspace(p0, be, ns)]
            pts.append(xy(r_mid, ang(p1)))
            pts += [xy(r_in, ang(p)) for p in np.linspace(be, p0, ns)]
        else:
            bs = p0 + ah
            pts.append(xy(r_mid, ang(p0)))
            pts += [xy(r_out, ang(p)) for p in np.linspace(bs, p1, ns)]
            pts += [xy(r_in, ang(p)) for p in np.linspace(p1, bs, ns)]
        ax.add_patch(Polygon(pts, closed=True, facecolor=color, edgecolor="black", lw=0.3, zorder=3))

    def radial_labels(ax, items, r_edge, r_text, fontsize=7.0):
        """items: list of (true_angle, name). Radial labels + leader lines, spread to avoid overlap.
        r_edge = gene-band edge the leader starts from; r_text = radius where the label text sits."""
        if not items:
            return
        items = sorted(items, key=lambda t: t[0])
        true_a = [t[0] for t in items]
        span = 2 * XLIM
        unit_pt = (FIGSIZE * 72.0) / span          # points per data unit
        dmin = (fontsize * 1.35 / unit_pt) / abs(r_text)
        adj = _spread(true_a, dmin)
        outward = r_text > r_edge
        for (ta, nm), aa in zip(items, adj):
            r_lead = r_text - 0.02 if outward else r_text + 0.02
            x0, y0 = xy(r_edge, ta)
            x1, y1 = xy(r_lead, aa)
            ax.plot([x0, x1], [y0, y1], color="#999", lw=0.4, zorder=2)
            thp = aa + ROT
            rot = math.degrees(math.atan2(math.cos(thp), math.sin(thp)))
            ha = "left" if outward else "right"
            if rot > 90 or rot < -90:
                rot += 180
                ha = "right" if outward else "left"
            tx, ty = xy(r_text, aa)
            ax.text(tx, ty, nm, fontsize=fontsize, rotation=rot, rotation_mode="anchor",
                    ha=ha, va="center", style="italic", zorder=4)

    def curved_text(ax, s, a_center, radius, fontsize=11):
        """Draw text bent along the ring (each glyph rotated tangentially)."""
        unit_pt = (FIGSIZE * 72.0) / (2 * XLIM)
        cw = (fontsize * 0.80 / unit_pt) / radius          # radians per character
        flip = radius * math.cos(a_center + ROT) < 0       # bottom half -> flip to stay upright
        chars = s[::-1] if flip else s
        a0 = a_center - cw * (len(chars) - 1) / 2
        for i, ch in enumerate(chars):
            a = a0 + i * cw
            thp = a + ROT
            rot = math.degrees(math.atan2(math.cos(thp), math.sin(thp))) - 90 + (180 if flip else 0)
            x, y = xy(radius, a)
            ax.text(x, y, ch, fontsize=fontsize, rotation=rot, ha="center", va="center",
                    rotation_mode="anchor", fontweight="bold", color="#222", zorder=6)

    FIGSIZE = 15; XLIM = 1.58
    fig, ax = plt.subplots(figsize=(FIGSIZE, FIGSIZE))
    cats = {}

    # structural ring in the MIDDLE, with region names curved ON the ring
    if ir:
        jlb, jsb, jsa, jla, lsc, irlen, ssc = ir
        for lab, a, b, col in [("IRb", jlb, jsb, "#fdae6b"), ("SSC", jsb, jsa, "#9ecae1"),
                               ("IRa", jsa, jla, "#fdae6b"), ("LSC", jla, jlb, "#c7e9c0")]:
            arc(ax, REG_R, a, b, col, REG_LW, z=1)
            mid = a + ((b - a) % L) / 2
            curved_text(ax, lab, ang(mid % L), REG_R, fontsize=11)

    # gene arrows + collect labels per ring (ALL genes; multi-exon genes drawn exon-by-exon)
    fwd, rev = [], []
    for f in rec.features:
        if f.type != "gene":
            continue
        nm = f.qualifiers.get("gene", f.qualifiers.get("locus_tag", ["?"]))[0]
        cat, col = categorize(nm); cats[cat] = col
        strand = 1 if f.location.strand == 1 else -1
        parts = list(f.location.parts)            # exons (1 part for simple genes)
        for p in parts:                            # draw each exon as an arrow
            s, e = int(p.start), int(p.end)
            if e - s <= 0 or e - s > L * 0.4:
                continue
            if strand == 1:
                gene_arrow(ax, s, e, 1, FWD_IN, FWD_OUT, col)
            else:
                gene_arrow(ax, s, e, -1, REV_IN, REV_OUT, col)
        big = max(parts, key=lambda p: int(p.end) - int(p.start))   # one label, at largest exon
        a_mid = ang((int(big.start) + int(big.end)) / 2)
        (fwd if strand == 1 else rev).append((a_mid, nm))

    radial_labels(ax, fwd, FWD_OUT, LAB_OUT)      # forward -> labels outside
    radial_labels(ax, rev, REV_IN, LAB_IN)        # reverse -> labels inside

    # centre: scientific name in italics (space preserved), accession, total size
    org_it = "$\\mathit{" + organism.replace(" ", r"\ ") + "}$"
    ax.text(0, 0, f"{org_it}\n{rec.id}\n{L:,} bp", ha="center", va="center", fontsize=13)

    order = ["Photosynthesis", "NADH dehydrogenase", "Ribosomal proteins", "RNA polymerase",
             "Transfer RNA (tRNA)", "Ribosomal RNA (rRNA)", "Other genes (ycf, matK, ...)"]
    leg = [Patch(facecolor=cats[c], label=c) for c in order if c in cats]
    if ir:
        leg += [Patch(facecolor="#c7e9c0", label=f"Large single-copy, LSC ({lsc:,} bp)"),
                Patch(facecolor="#fdae6b", label=f"Inverted repeat, IRa/IRb ({irlen:,} bp each)"),
                Patch(facecolor="#9ecae1", label=f"Small single-copy, SSC ({ssc:,} bp)")]
    ax.legend(handles=leg, loc="upper left", fontsize=11.5, frameon=False,
              bbox_to_anchor=(-0.02, 0.04), ncol=2,
              title="Gene functional class / region", title_fontsize=12.5)
    ax.text(0, XLIM*0.97, "Genes outside the circle: (+) strand   |   inside: (-) strand",
            ha="center", fontsize=12, color="#444")
    ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-XLIM, XLIM); ax.set_ylim(-XLIM, XLIM)
    ax.set_title(f"Plastome map of {org_it} ({rec.id})", fontsize=17, pad=8)
    plt.savefig(f"{out}.png", dpi=dpi, bbox_inches="tight", pad_inches=0.1)   # raster (high-dpi)
    plt.savefig(f"{out}.pdf", bbox_inches="tight", pad_inches=0.1)            # vector
    plt.savefig(f"{out}.svg", bbox_inches="tight", pad_inches=0.1)            # vector
    plt.close()
    print(f"Saved {out}.{{png,pdf,svg}}  ({L:,} bp; PNG @ {dpi} dpi; IR ring: {'yes' if ir else 'no'})")


def main():
    ap = argparse.ArgumentParser(description="Standalone circular plastome map from a GenBank file.")
    ap.add_argument("genbank", help="input GenBank (.gb) file with sequence + gene features")
    ap.add_argument("--out", help="output basename (default: <accession>_map)")
    ap.add_argument("--rotate", type=float, default=180,
                    help="rotate whole map by degrees (default 180 -> LSC on the left)")
    ap.add_argument("--ir", help="supply IR junctions instead of detecting: jlb,jsb,jsa,jla (1-based)")
    ap.add_argument("--organism", help="override organism name in title")
    ap.add_argument("--dpi", type=int, default=600,
                    help="raster (PNG) resolution; default 600 (line-art journal standard). "
                         "PDF/SVG are vector and resolution-independent.")
    a = ap.parse_args()

    rec = next(SeqIO.parse(a.genbank, "genbank"))
    seq = str(rec.seq).upper()
    gb_org = (rec.annotations.get("organism") or "").strip()
    organism = a.organism or (gb_org if _slug(gb_org) else rec.id)

    if a.ir:
        jlb, jsb, jsa, jla = (int(x) for x in a.ir.split(","))
        ir = (jlb, jsb, jsa, jla, 0, jsb - jlb, jsa - jsb)   # lengths approximate
    else:
        ir = ir_from_blast(seq) or ir_from_annotation(rec)
        if ir is None:
            print("[warn] IR not detected (no blastn / no repeat_region) -> drawing without region ring",
                  file=sys.stderr)

    # Descriptive output basename: <Organism>_<accession>_plastome_map (organism
    # from the GenBank ORGANISM header when present), e.g.
    # Nicotiana_tabacum_NC_001879.2_plastome_map — never a generic name.
    if a.out:
        out = a.out
    else:
        org_slug = _slug(gb_org)
        out = (f"{org_slug}_{rec.id}_plastome_map" if org_slug
               else f"{rec.id}_plastome_map")
    draw_map(rec, ir, organism, a.rotate, out, dpi=a.dpi)


if __name__ == "__main__":
    main()
