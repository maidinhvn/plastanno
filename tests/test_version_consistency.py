#!/usr/bin/env python3
"""The version is written twice; it must be written the same twice.

A dynamic `attr: plastanno.__version__` in pyproject.toml would remove the
duplication, and the module docstring used to claim that was how it worked. It
is not: this repository ships both a `plastanno/` package and a `plastanno.py`
entry script, so setuptools resolves the name to the script and the build fails
with "'plastanno' is not a package". The two literals therefore stay, and this
test is what keeps them from drifting — a stale version once made `--version`
report an old release while the code was newer.
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
    print("  %s %s" % ("ok  " if ok else "FAIL", label))
    if not ok:
        FAIL.append("%s: got %r want %r" % (label, got, want))


pyproject = open(os.path.join(ROOT, "pyproject.toml")).read()
m = re.search(r'^version\s*=\s*"([^"]+)"', pyproject, re.M)
RUN[0] += 1
if m:
    print("  ok   pyproject.toml declares a version")
else:
    print("  FAIL pyproject.toml declares a version")
    FAIL.append("no version in pyproject.toml")
    print("\n%d checks, %d failed" % (RUN[0], len(FAIL)))
    sys.exit(1)
declared = m.group(1)

from plastanno import __version__                              # noqa: E402

check("pyproject.toml and plastanno.__version__ agree", declared, __version__)
check("  ... and it is a three-part version",
      bool(re.fullmatch(r"\d+\.\d+\.\d+", __version__)), True)

# The CHANGELOG must carry an entry for whatever version ships, or a release
# goes out with no record of what changed in it.
changelog = open(os.path.join(ROOT, "CHANGELOG.md")).read()
check("CHANGELOG.md has an entry for this version",
      ("[%s]" % __version__) in changelog, True)

print()
print("%d checks, %d failed" % (RUN[0], len(FAIL)))
for f in FAIL:
    print("  FAILED: %s" % f)
sys.exit(1 if FAIL else 0)
