#!/usr/bin/env python3
"""Every CLI flag is documented, and every documented flag exists.

The README described the pre-3.0.0 tool for three months. It called tRNAscan-SE
optional after it became required, omitted numpy and scipy, and said "four
external tools" when there were five. None of that was caught by anything,
because nothing compared the documentation to the code.

This compares them in both directions:

  * a flag the CLI accepts and the README never mentions is undocumented — users
    cannot find it, and a flag that changes annotation output is worth finding;
  * a flag the README shows and the CLI does not accept is a broken instruction.

The second direction needs an allowlist, because the README legitimately shows
commands for other programs. Each entry names the program it belongs to, so the
allowlist stays honest instead of becoming a place to hide failures.
"""
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

RUN = [0]
FAIL = []


def check(label, got, want):
    RUN[0] += 1
    ok = got == want
    if not ok:
        FAIL.append("%s (got %r, want %r)" % (label, got, want))
    print("  %s %s" % ("ok  " if ok else "FAIL", label))


# Flags the CLI defines, read from the source rather than by running argparse,
# so this does not depend on an importable environment.
cli_src = open(os.path.join(ROOT, "plastanno", "cli.py")).read()
cli_flags = set(re.findall(r"add_argument\(\s*[\"'](--[a-z0-9-]+)", cli_src))
check("the CLI defines flags at all", len(cli_flags) > 5, True)

readme = open(os.path.join(ROOT, "README.md")).read()
readme_flags = set(re.findall(r"--[a-z0-9][a-z0-9-]+", readme))

# Flags belonging to other programs, each with the program named.
FOREIGN = {
    "--help": "universal",
    "--override-channels": "conda",
    "--set": "scripts/benchmark/multi_genome_bench.py",
    "--sim": "scripts/benchmark/benchmark_gene_by_gene.py",
    "--tol": "scripts/benchmark/benchmark_gene_by_gene.py",
    "--workers": "scripts/benchmark/multi_genome_bench.py",
    "--final": "scripts/benchmark/multi_genome_bench.py",
    "--keep": "scripts/benchmark/multi_genome_bench.py",
    "--no-deps": "pip",
    "--no-build-isolation": "pip",
}

print("every CLI flag is documented")
for f in sorted(cli_flags):
    check("README mentions %s" % f, f in readme_flags, True)

print()
print("every documented flag is real")
for f in sorted(readme_flags - cli_flags):
    check("%s belongs to a named program%s"
          % (f, "" if f in FOREIGN else " — UNKNOWN"),
          f in FOREIGN, True)

print()
print("the allowlist does not hide real flags")
for f in sorted(FOREIGN):
    check("%s is genuinely not a plastanno CLI flag" % f, f in cli_flags, False)

print()
print("%d checks, %d failed" % (RUN[0], len(FAIL)))
for f in FAIL:
    print("  FAILED: %s" % f)
sys.exit(1 if FAIL else 0)
