#!/usr/bin/env python3
"""'tRNA-???' must cost one candidate, not the whole genome.

ARAGORN writes `tRNA-???` when it finds a tRNA structure it cannot assign an
amino acid to. That is ordinary output. The fail-loud parser added after Gate A
treated it as an unreadable coordinate form and raised, aborting the annotation:
three of sixty genomes died this way in
`benchmark_v3/intron_glocal/baseline_exits.txt`.

These checks pin the narrow behaviour: the unassignable line is skipped, every
other line on either side of it still parses, and the genuinely malformed lines
the raise exists for still raise.
"""
import io
import os
import sys
from contextlib import redirect_stdout

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from plastanno.identify.engine_b import parse_aragorn_output   # noqa: E402

RUN = [0]
FAIL = []


def check(label, got, want):
    RUN[0] += 1
    ok = got == want
    print("  %s %s" % ("ok  " if ok else "FAIL", label))
    if not ok:
        FAIL.append("%s: got %r want %r" % (label, got, want))


def raises(label, fn, *a):
    RUN[0] += 1
    try:
        fn(*a)
    except ValueError as exc:
        print("  ok   %s" % label)
        return str(exc)
    print("  FAIL %s (no ValueError)" % label)
    FAIL.append(label)
    return ""


L = 160000

# The two real lines that killed NC_034893.1 and NC_077615.1, verbatim.
UNK_FWD = "24  tRNA-???                 [99767,99838]\t33  \t(tatg)\n"
UNK_REV = "23  tRNA-???                c[65107,65182]\t35  \t(catt)\n"
GOOD = ("1   tRNA-His                  c[156510,71]\t35  \t(gtg)\n"
        "2   tRNA-Lys                  c[1737,4320]\t33  \t(ttt)i(38,2511)\n")

print("tRNA-??? is skipped, not fatal")

buf = io.StringIO()
with redirect_stdout(buf):
    feats = parse_aragorn_output(GOOD + UNK_FWD, L)
check("a trailing 'tRNA-???' does not abort the parse", len(feats), 2)
check("  ... and it is announced, not silent", "???" in buf.getvalue(), True)

# Sandwiched: the lines AFTER the bad one must survive too. A `break` instead of
# a `continue` would pass the check above and fail this one.
buf = io.StringIO()
with redirect_stdout(buf):
    feats = parse_aragorn_output(UNK_FWD + GOOD + UNK_REV, L)
check("lines after an unassignable one still parse", len(feats), 2)
check("both orientations of '???' are skipped",
      buf.getvalue().count("???"), 2)

# Only '???' is forgiven. The malformed coordinate the raise was added for in
# the first place must still raise, or this fix has widened into the hole it was
# meant to leave open.
BAD = GOOD + "6   tRNA-His                      c[-3,71]\t35  \t(gtg)\n"
m = raises("the linear-mode negative start still raises",
           parse_aragorn_output, BAD, L)
check("  ... still quoting the offending line", "-3" in m, True)

# A '???' line whose coordinates are themselves impossible is still only skipped:
# the amino acid is what makes it unusable, so there is nothing to validate.
buf = io.StringIO()
with redirect_stdout(buf):
    feats = parse_aragorn_output(GOOD + "9  tRNA-???  [0,12]\t33  \t(ttt)\n", L)
check("'???' short-circuits before coordinate checks", len(feats), 2)

check("header and 'N genes found' still ignored",
      len(parse_aragorn_output(">x\n12 genes found\n\n\n", L)), 0)

print()
print("%d checks, %d failed" % (RUN[0], len(FAIL)))
for f in FAIL:
    print("  FAILED: %s" % f)
sys.exit(1 if FAIL else 0)
