#!/usr/bin/env python3
"""An exon panel that cannot be searched is reported once, not skipped in silence.

Run DIRECTLY and check $?. Exit 0 = all pass.

refine_splice places the junctions of the panel genes by BLAST against the packaged exon panel.
When that database was missing or unreadable, _blast_exons returned nothing without a word, and
junction refinement was skipped for every panel gene while the run looked normal. That happened
during development, when a version-5 rebuild lost its index files. Now the first failure prints
one WARNING saying why, and later calls stay quiet.
"""
import contextlib
import io
import os
import shutil
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from plastanno.annotate import refine_splice as RS           # noqa: E402

FAIL, RUN = [], [0]


def check(label, got, want):
    RUN[0] += 1
    ok = got == want
    print("  %-70s %s" % (label, "ok" if ok else "FAIL got=%r want=%r" % (got, want)))
    if not ok:
        FAIL.append(label)


def calls(n, region="ACGT" * 30, gene="ndhA"):
    """Run _blast_exons n times from a fresh state; return (results, printed warnings)."""
    RS._PANEL_PROBLEM.clear()
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        res = [RS._blast_exons(region, gene) for _ in range(n)]
    return res, [l for l in buf.getvalue().splitlines() if "WARNING" in l]


REAL = RS._BDB
tmp = tempfile.mkdtemp(prefix="panelwarn_")
try:
    print("--- a missing panel ---")
    RS._BDB = os.path.join(tmp, "nothing", "exons")
    res, warn = calls(3)
    check("nothing is returned", res, [{}, {}, {}])
    check("one warning for three calls", len(warn), 1)
    check("  ... saying the panel database is missing", "database is missing" in (warn[0] if warn else ""), True)
    check("  ... and that refinement is skipped", "refinement is skipped" in (warn[0] if warn else ""), True)

    if shutil.which("blastn"):
        print("--- an unreadable panel ---")
        bad = os.path.join(tmp, "bad"); os.makedirs(bad)
        for e in ("nhr", "nin", "nsq"):
            open(os.path.join(bad, "exons." + e), "wb").write(b"not a blast database")
        RS._BDB = os.path.join(bad, "exons")
        res, warn = calls(2)
        check("nothing is returned", res, [{}, {}])
        check("one warning for two calls", len(warn), 1)
        check("  ... naming BLAST's own error", "BLAST could not search" in (warn[0] if warn else ""), True)

        print("--- the packaged panel ---")
        RS._BDB = REAL
        exon = None
        for line in open(REAL + ".fasta"):
            if exon is None and line.startswith(">ndhA|e0|0"):
                exon = ""
            elif exon == "" and not line.startswith(">"):
                exon = line.strip()
        region = "GCC" * 30 + exon + "GTGCGTACGT" + "GCC" * 30
        res, warn = calls(1, region=region, gene="ndhA")
        check("a real ndhA exon is found", 0 in res[0], True)
        check("  ... without any warning", warn, [])
    else:
        print("  (blastn not on PATH: the unreadable-panel and packaged-panel cases are not run)")
finally:
    RS._BDB = REAL
    RS._PANEL_PROBLEM.clear()
    shutil.rmtree(tmp, ignore_errors=True)

print()
print("%d checks, %d failed" % (RUN[0], len(FAIL)))
for f in FAIL:
    print("  FAILED: %s" % f)
sys.exit(1 if FAIL else 0)
