#!/usr/bin/env python3
"""The circular map italicises the organism name with mathtext, so the name has
to be escaped first.

Without escaping, every character mathtext treats as markup was interpreted as
markup. The visible symptom: when `--organism` is not given the accession is
used in its place, and `NC_053537.1` was drawn as "NC" with a subscript zero
followed by "53537.1" — in the map title and in the centre label, on every
output format. Underscores are common in accessions, so this was not a corner
case.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from plastanno.viz.plastome_circular_map import _mathit          # noqa: E402

RUN = [0]
FAIL = []


def check(label, got, want):
    RUN[0] += 1
    ok = got == want
    if not ok:
        FAIL.append("%s (got %r, want %r)" % (label, got, want))
    print("  %s %s" % ("ok  " if ok else "FAIL", label))


# The case that was actually broken.
check("underscore is escaped, not left as a subscript operator",
      _mathit("NC_053537.1"), r"$\mathit{NC\_053537.1}$")

# A real binomial still italicises with its space preserved.
check("a space becomes a mathtext space",
      _mathit("Panax ginseng"), r"$\mathit{Panax\ ginseng}$")

# Every other character mathtext would swallow.
for ch in "${}^_#&%~":
    out = _mathit("a" + ch + "b")
    check("%r is escaped" % ch, "\\" + ch in out, True)

check("a backslash becomes the mathtext symbol",
      _mathit("a\\b"), r"$\mathit{a\backslash b}$")

# Escaping is worthless if the result no longer renders. Parse each one the way
# matplotlib will when it draws the figure.
for name in ["NC_053537.1", "Panax ginseng", "a_b^c{d}$e",
             "Daucus carota subsp. sativus", "a\\b", "100%_x"]:
    try:
        TextPath((0, 0), _mathit(name), size=10, prop=FontProperties())
        ok = True
    except Exception as exc:                       # noqa: BLE001
        ok = "%s: %s" % (type(exc).__name__, exc)
    check("mathtext renders %r" % name, ok, True)

# Adversarial: the unescaped form must actually render DIFFERENTLY, otherwise
# this test would pass even if _mathit did nothing useful.
raw = "$\\mathit{" + "NC_053537.1" + "}$"
p_raw = TextPath((0, 0), raw, size=10, prop=FontProperties())
p_esc = TextPath((0, 0), _mathit("NC_053537.1"), size=10, prop=FontProperties())
check("escaped and unescaped really differ (the bug was real)",
      p_raw.get_extents().bounds != p_esc.get_extents().bounds, True)

# Guard against someone rebuilding the string inline again, which is how the
# bug got in: the mathtext wrapper must exist in exactly one place.
src = open(os.path.join(ROOT, "plastanno", "viz",
                        "plastome_circular_map.py")).read()
check("the mathtext wrapper is built in one place only",
      src.count('"$\\\\mathit{'), 1)

print()
print("%d checks, %d failed" % (RUN[0], len(FAIL)))
for f in FAIL:
    print("  FAILED: %s" % f)
sys.exit(1 if FAIL else 0)
