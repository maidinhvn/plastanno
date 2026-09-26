#!/usr/bin/env python3
"""`--trna-mode hybrid` is the default, and its dependency fails early and loudly.

The evidence for the flip is in benchmark_v3/holdout/RESULT_HYBRID.md: on 99
held-out genomes from 89 families, the intron-free acceptor stem improves by
+1.124 [+1.065, +1.182], and by +1.163 restricted to the 67 genomes from families
the development set never contained, with intron-bearing loci, CDS, rRNA and the
whole tRNA inventory unmoved.

The flip has a cost that must not be paid silently: tRNAscan-SE becomes a hard
dependency. Two things are pinned here.

  1. The default is hybrid in every place that declares one, so there is no second
     "default" waiting to diverge. `pipeline.run()` had kept `legacy` in its own
     signature while the CLI said otherwise, which is the same shape of defect as
     the missing train.csv: two sources of truth, one of them silent.
  2. The dependency is checked BEFORE any work begins. Learning that the binary is
     absent after Exonerate, the HMM search and reconciliation have all run is
     several wasted minutes per genome, and across a batch it is hours.
"""
import os
import shutil
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

RUN = [0]
FAIL = []


def check(label, got, want):
    RUN[0] += 1
    ok = got == want
    print("  %-68s %s" % (label, "ok" if ok else "FAIL got=%r want=%r" % (got, want)))
    if not ok:
        FAIL.append(label)


print("--- 1. every declared default says hybrid ---")
import plastanno.cli as cli                              # noqa: E402
import inspect                                           # noqa: E402
from plastanno.pipeline import run, _write_provenance    # noqa: E402

parser = cli.build_parser() if hasattr(cli, "build_parser") else None
if parser is None:
    # The parser is built inside main(); read the declarations from the source
    # instead of restructuring production code to suit a test.
    src = inspect.getsource(cli)
    check("no parser still declares legacy as the tRNA default",
          'choices=["legacy", "hybrid"], default="legacy"' in src, False)
    check("both parsers declare hybrid",
          src.count('choices=["legacy", "hybrid"], default="hybrid"'), 2)

check("pipeline.run() signature default",
      inspect.signature(run).parameters["trna_mode"].default, "hybrid")
check("_write_provenance would not stamp the wrong mode",
      inspect.signature(_write_provenance).parameters["trna_mode"].default, "hybrid")

print("--- 2. the missing dependency is reported usefully ---")
from plastanno.identify.trna_hybrid import (             # noqa: E402
    require_trnascan, TrnascanUnavailable,
)

real_which = shutil.which
try:
    shutil.which = lambda *a, **k: None                  # pretend it is absent
    msg = ""
    try:
        require_trnascan()
    except TrnascanUnavailable as exc:
        msg = str(exc)
    RUN[0] += 1
    ok = bool(msg)
    print("  %-68s %s" % ("require_trnascan raises when the binary is absent",
                          "ok" if ok else "FAIL did not raise"))
    if not ok:
        FAIL.append("require_trnascan raises when the binary is absent")
    check("the message names the tool", "tRNAscan-SE" in msg, True)
    check("it says how to install it", "conda install" in msg, True)
    check("it says how to opt out", "--trna-mode legacy" in msg, True)
    check("it says hybrid is the default, so the user knows why it fired",
          "DEFAULT" in msg, True)

    print("--- 3. the check happens before any work ---")
    # A fasta that does not exist. If the dependency were checked at step 6b we
    # would get a file error; getting TrnascanUnavailable proves the ordering.
    with tempfile.TemporaryDirectory() as d:
        kind = ""
        try:
            run(os.path.join(d, "no-such-genome.fasta"), os.path.join(d, "out"),
                trna_mode="hybrid")
        except TrnascanUnavailable:
            kind = "dependency"
        except Exception as exc:                          # noqa: BLE001
            kind = "%s: %s" % (type(exc).__name__, str(exc)[:60])
        check("run() refuses on the dependency, not on the missing input",
              kind, "dependency")

    print("--- 4. legacy still runs without the dependency ---")
    with tempfile.TemporaryDirectory() as d:
        kind = ""
        try:
            run(os.path.join(d, "no-such-genome.fasta"), os.path.join(d, "out"),
                trna_mode="legacy")
        except TrnascanUnavailable:
            kind = "dependency"
        except Exception:                                 # noqa: BLE001
            kind = "got past the dependency check"
        check("legacy is not gated on tRNAscan-SE",
              kind, "got past the dependency check")
finally:
    shutil.which = real_which

print("--- 5. an unknown mode is still refused ---")
with tempfile.TemporaryDirectory() as d:
    kind = ""
    try:
        run(os.path.join(d, "x.fasta"), os.path.join(d, "out"), trna_mode="turbo")
    except ValueError as exc:
        kind = "rejected" if "turbo" in str(exc) else "rejected without naming it"
    except Exception:                                     # noqa: BLE001
        kind = "wrong exception"
    check("a typo does not fall through to a default", kind, "rejected")


print()
print("%d checks, %d failed" % (RUN[0], len(FAIL)))
for f in FAIL:
    print("  FAILED: %s" % f)
sys.exit(1 if FAIL else 0)
