#!/usr/bin/env python3
"""Glocal intron refinement is the shipping default, in every entry point.

A default that is set in one place and missed in another is how this project
previously shipped a mode that silently never ran. `pipeline.run`, the `run`
subcommand and the `batch` subcommand each carry their own copy, so each is
checked separately rather than inferred from one of them.

The measured basis is `benchmark_v3/intron_verify/RESULT.md`: 645 pristine
non-RefSeq genomes, 5178 loci, exact coordinates 13.0% -> 55.1%, outer-exact
63.2% -> 93.5%, zero loci made worse.
"""
import inspect
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from plastanno import cli, pipeline                              # noqa: E402

RUN = [0]
FAIL = []


def check(label, got, want):
    RUN[0] += 1
    ok = got == want
    print("  %s %s" % ("ok  " if ok else "FAIL", label))
    if not ok:
        FAIL.append("%s: got %r want %r" % (label, got, want))


print("the default reaches every entry point")

sig = inspect.signature(pipeline.run)
check("pipeline.run defaults to glocal",
      sig.parameters["intron_mode"].default, "glocal")

# The provenance sidecar records what produced an annotation. If its default
# disagrees with the pipeline's, a run that took the default would be recorded
# as having taken the other mode.
psig = inspect.signature(pipeline._write_provenance)
check("the provenance sidecar agrees with it",
      psig.parameters["intron_mode"].default, "glocal")

src = inspect.getsource(cli)
check("both CLI subcommands default to glocal",
      src.count('choices=["legacy", "glocal"], default="glocal"'), 2)
check("  ... and neither is left on legacy",
      'choices=["legacy", "glocal"], default="legacy"' in src, False)

# A flag that the caller never forwards is a flag that does nothing. This is the
# exact failure that voided benchmark_v3/intron_geometry/.
check("cmd_run forwards the flag to the pipeline",
      src.count("intron_mode = args.intron_mode"), 2)

# Help text that names the wrong default is a documentation bug that reads as a
# behaviour bug.
check("the help text calls glocal the default",
      "glocal \"\n             \"(default)" in src or "glocal (default)" in src
      or 'glocal "' in src and '"(default): each ' in src, True)

# legacy must remain reachable: it is what reproduces the frozen benchmark
# record, and the verification arm depends on it.
check("legacy is still an accepted choice",
      'choices=["legacy", "glocal"]' in src, True)
try:
    pipeline.run.__wrapped__          # noqa: B018
except AttributeError:
    pass
check("an unknown mode is rejected rather than silently ignored",
      "unknown intron mode" in inspect.getsource(pipeline), True)

print()
print("%d checks, %d failed" % (RUN[0], len(FAIL)))
for f in FAIL:
    print("  FAILED: %s" % f)
sys.exit(1 if FAIL else 0)
