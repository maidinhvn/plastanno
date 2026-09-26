#!/usr/bin/env python3
"""Two hardening assertions added after Gate A closed.

Gate A showed that geometry survives a change of representation. It did not
constrain what happens when ARAGORN is invoked wrongly, or when a coordinate
arrives that cannot be a position on the circle. Both were silent failures:
dropping `-c` made every origin-crossing tRNA vanish without a message, and a
negative start was skipped along with the header lines.

These tests check that the two now fail loudly, and — the part that matters for
a tool already measured — that valid input is unaffected.
"""
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from plastanno.identify.engine_b import (          # noqa: E402
    ARAGORN_ARGS, check_aragorn_args, aragorn_coords, parse_aragorn_output,
)

RUN = [0]
FAIL = []


def check(label, got, want):
    RUN[0] += 1
    ok = got == want
    print("  %-66s %s" % (label, "ok" if ok else "FAIL got=%r want=%r" % (got, want)))
    if not ok:
        FAIL.append(label)


def raises(label, fn, *a, **kw):
    """The call must raise, and the message must name the reason."""
    RUN[0] += 1
    try:
        fn(*a, **kw)
    except (ValueError, RuntimeError) as exc:
        msg = str(exc)
        informative = len(msg) > 20
        print("  %-66s %s" % (label, "ok" if informative else "FAIL silent: %r" % msg))
        if not informative:
            FAIL.append(label)
        return msg
    except Exception as exc:                        # noqa: BLE001
        print("  %-66s FAIL wrong type %s" % (label, type(exc).__name__))
        FAIL.append(label)
        return ""
    print("  %-66s FAIL did not raise" % label)
    FAIL.append(label)
    return ""


L = 156_513                                         # a real plastome length


print("--- 1. the production invocation is locked to circular mode ---")
check("ARAGORN_ARGS contains -c", "-c" in ARAGORN_ARGS, True)
check("ARAGORN_ARGS does not contain -l", "-l" in ARAGORN_ARGS, False)
check("ARAGORN_ARGS is exactly the reviewed set",
      tuple(ARAGORN_ARGS), ("-i", "-t", "-c", "-w"))
check("check_aragorn_args accepts the production set",
      tuple(check_aragorn_args()), tuple(ARAGORN_ARGS))

m = raises("dropping -c is refused", check_aragorn_args, ("-i", "-t", "-w"))
check("  ... and the message says which flag", "-c" in m, True)
m = raises("adding -l is refused", check_aragorn_args, ("-i", "-t", "-c", "-l", "-w"))
check("  ... and the message says linear mode is unsupported", "-l" in m, True)

# The lock has to be on the path the pipeline actually takes, not only on a
# constant a caller may ignore. Read the call site back.
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "..", "plastanno", "identify", "engine_b.py")).read()
# The argv expression: from the "aragorn" literal to the end of that call's
# argument, which ends at `capture_output`. Matching to the first ']' would stop
# inside `["aragorn"] + list(...)`.
call = re.search(r'"aragorn".{0,200}?capture_output', src, re.S)
argv = call.group(0) if call else ""
# Count INVOCATIONS, not occurrences of the word. `--exon-mode aragorn` puts the
# string "aragorn" in the source as a mode name too, and a test that cannot tell a
# mode name from a subprocess call is testing the wrong thing.
check("the pipeline invokes the aragorn binary exactly once",
      src.count('["aragorn"]'), 1)
check("run_aragorn builds its argv from check_aragorn_args()",
      "check_aragorn_args()" in argv, True)
check("no hard-coded aragorn flag list survives at the call site",
      '"-c"' not in argv and '"-i"' not in argv, True)


print("--- 2. the parser refuses a coordinate that cannot be on the circle ---")
m = raises("negative start (the -l form 'c[-3,71]')", aragorn_coords, "-3", "71", L)
check("  ... message names the allowed range", str(L) in m, True)
raises("zero start (1-based positions start at 1)", aragorn_coords, "0", "71", L)
raises("start past the end of the genome", aragorn_coords, str(L + 1), "71", L)
raises("end past the end of the genome", aragorn_coords, "100", str(L + 1), L)
raises("both out of range", aragorn_coords, "-5", "-1", L)

print("--- 3. valid coordinates are translated exactly as before ---")
# These are the conversions Gate A locked; the range check must not perturb them.
check("ordinary feature [101,172] -> one arc",
      aragorn_coords("101", "172", 1000), ([(100, 172)], False))
check("origin-crossing [971,40] -> two arcs, wrapped",
      aragorn_coords("971", "40", 1000), ([(970, 1000), (0, 40)], True))
check("first base of the genome is accepted",
      aragorn_coords("1", "72", 1000), ([(0, 72)], False))
check("last base of the genome is accepted",
      aragorn_coords("929", "1000", 1000), ([(928, 1000)], False))
check("a feature ending exactly at the origin seam",
      aragorn_coords("1000", "40", 1000), ([(999, 1000), (0, 40)], True))


print("--- 4. a hit-shaped line the parser cannot read is refused, not skipped ---")
GOOD = (">NC_073016.1 Trevesia palmata chloroplast, complete genome\n"
        "37 genes found\n"
        "1   tRNA-His                  c[156510,71]\t35  \t(gtg)\n"
        "2   tRNA-Lys                  c[1737,4320]\t33  \t(ttt)i(38,2511)\n"
        "5   tRNA-Ser                 [10171,10940]\t712 \t(cga)i(33,678)\n")
feats = parse_aragorn_output(GOOD, L)
check("the real ARAGORN -c output parses", len(feats), 3)
check("the wrapped tRNA keeps start > end", feats[0].start > feats[0].end, True)

BAD = GOOD + "6   tRNA-His                      c[-3,71]\t35  \t(gtg)\n"
m = raises("a line with the linear-mode negative start", parse_aragorn_output, BAD, L)
check("  ... the message quotes the offending line", "-3" in m, True)

# The header, the count line and blanks must still be ignored: a fail-loud
# parser that chokes on ARAGORN's own preamble would be worse than the hole.
check("header and 'N genes found' still ignored",
      len(parse_aragorn_output(">x\n12 genes found\n\n\n", L)), 0)
check("an unrelated line is still ignored",
      len(parse_aragorn_output("Total 37 tRNA genes\n", L)), 0)


print()
print("%d checks, %d failed" % (RUN[0], len(FAIL)))
for f in FAIL:
    print("  FAILED: %s" % f)
sys.exit(1 if FAIL else 0)
