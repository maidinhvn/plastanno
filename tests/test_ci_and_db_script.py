#!/usr/bin/env python3
"""The download location, and the environment CI builds, must not drift.

Two failures this guards against, both of which already happened:

  * `scripts/get_database.sh` carried its own copy of the Zenodo URL and MD5.
    When the database was republished, `plastanno/fetch_db.py` was updated and
    the script was not. The stale record still existed and its stale checksum
    still matched, so the script downloaded the *previous* database and exited
    0 — CI installed a bundle missing eight genes and said nothing.

  * The CI environment listed only the dependencies 2.x needed. 3.0.0 added
    scipy/numpy (reconcile) and made `--trna-mode hybrid` the default, which
    requires tRNAscan-SE before step 1. CI went red on the release commit.

So: the script must *derive* the constants rather than repeat them, and the
workflow must install everything pyproject declares plus every external tool.
"""
import os
import re
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(ROOT, "scripts", "get_database.sh")
WORKFLOW = os.path.join(ROOT, ".github", "workflows", "ci.yml")

RUN = [0]
FAIL = []


def check(label, got, want):
    RUN[0] += 1
    ok = got == want
    if not ok:
        FAIL.append("%s (got %r, want %r)" % (label, got, want))
    print("  %s %s" % ("ok  " if ok else "FAIL", label))


script_src = open(SCRIPT).read()
workflow_src = open(WORKFLOW).read()

# ---------------------------------------------------------------- the script
print("scripts/get_database.sh")

# No second copy of the constants. A bare "zenodo.org/records/<digits>" or a
# 32-hex literal in the script means someone reintroduced the duplication.
check("no hard-coded Zenodo record in the script",
      bool(re.search(r"zenodo\.org/records/\d+", script_src)), False)
check("no hard-coded MD5 literal in the script",
      bool(re.search(r'"[0-9a-f]{32}"', script_src)), False)

# Run the real extraction lines — not a reimplementation of them — against a
# synthetic project root, so the test exercises the shell the script ships.
head = script_src.split("TAR=")[0] + '\necho "URL:$URL"\necho "MD5:$MD5"\n'


def run_head(fetch_db_body):
    """Run the script's header against a fake fetch_db.py; return (rc, out)."""
    tmp = tempfile.mkdtemp()
    os.makedirs(os.path.join(tmp, "scripts"))
    os.makedirs(os.path.join(tmp, "plastanno"))
    with open(os.path.join(tmp, "scripts", "get_database.sh"), "w") as fh:
        fh.write(head)
    with open(os.path.join(tmp, "plastanno", "fetch_db.py"), "w") as fh:
        fh.write(fetch_db_body)
    p = subprocess.run(["bash", os.path.join(tmp, "scripts", "get_database.sh")],
                       capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


# Adversarial: a URL the script has never seen. If it were hard-coded, or if
# the sed pattern were anchored to the real record, this would not come back.
rc, out = run_head('URL = "https://example.invalid/records/999/x.tar.gz"\n'
                   'MD5 = "ffffffffffffffffffffffffffffffff"\n')
check("reads an arbitrary URL out of fetch_db.py", rc, 0)
check("  ... and it is the one fetch_db.py declares",
      "URL:https://example.invalid/records/999/x.tar.gz" in out, True)
check("  ... and the MD5 likewise",
      "MD5:ffffffffffffffffffffffffffffffff" in out, True)

# Adversarial: constants absent or renamed. Silently continuing with an empty
# URL would curl nothing and then "verify" an empty file, so it must abort.
rc, out = run_head("# no constants here\n")
check("aborts when the constants cannot be read", rc != 0, True)
check("  ... and says which file it looked in", "fetch_db.py" in out, True)

# And against the real file: what it extracts must equal what fetch-db uses.
sys.path.insert(0, ROOT)
from plastanno import fetch_db                                  # noqa: E402

real = subprocess.run(
    ["bash", "-c",
     'cd %s && sed -n \'s/^URL[[:space:]]*=[[:space:]]*"\\(.*\\)".*/\\1/p\' '
     'plastanno/fetch_db.py' % ROOT],
    capture_output=True, text=True).stdout.strip()
check("extraction from the real file matches fetch_db.URL", real, fetch_db.URL)

# No absolute path from a developer's machine may ship in the public tree. The
# repository rule has existed since 2.x ("remove personal paths -> $PLASTANNO_DATA")
# and was broken again by copying a benchmark test across from the dev tree.
print()
print("no personal paths")
# This rule is about the PUBLISHED tree only. The development tree legitimately
# contains absolute paths to a private dataset, and it is identifiable because
# it is the one that CONTAINS the public checkout as a subdirectory.
import glob
if os.path.isdir(os.path.join(ROOT, "forgit")):
    print("  skip  development tree — personal paths are allowed here")
    leaked = None
else:
    leaked = []
    for pat in ("plastanno/**/*.py", "scripts/**/*.py", "scripts/**/*.sh",
                "tests/*.py", "*.py", "*.sh"):
        for f in glob.glob(os.path.join(ROOT, pat), recursive=True):
            body = open(f, encoding="utf-8", errors="replace").read()
            # assembled at runtime so this file does not contain the literals it
            # searches for, which would make the check flag itself
            for marker in ("/data06/" + "users/", "/data06/" + "biotools/",
                           "/home/" + "vutrinh"):
                if marker in body:
                    leaked.append("%s -> %s" % (os.path.relpath(f, ROOT), marker))
    check("no developer-machine path in the published tree", leaked, [])

# --------------------------------------------------------------- the workflow
print()
print(".github/workflows/ci.yml")

m = re.search(r"conda create[^\n]*\n((?:\s+[^\n]*\\\n)*\s+[^\n]*)", workflow_src)
env_block = (m.group(0) if m else "")
env_pkgs = set(re.findall(r"[A-Za-z][A-Za-z0-9._-]+", env_block.replace("\\", " ")))

# Every runtime dependency pyproject declares must be installed by CI.
pyproject = open(os.path.join(ROOT, "pyproject.toml")).read()
deps_block = re.search(r"dependencies = \[(.*?)\]", pyproject, re.S).group(1)
declared = re.findall(r'"([A-Za-z][A-Za-z0-9._-]*)"', deps_block)
check("pyproject declares a non-empty dependency list", len(declared) > 0, True)
for d in sorted(declared):
    check("CI installs declared dependency %s" % d, d.lower() in
          {p.lower() for p in env_pkgs}, True)

# Every external tool the pipeline shells out to.
for tool, pkg in [("blastn", "blast"), ("exonerate", "exonerate"),
                  ("hmmsearch", "hmmer"), ("aragorn", "aragorn"),
                  ("tRNAscan-SE", "trnascan-se")]:
    check("CI installs %s (for %s)" % (pkg, tool), pkg in env_pkgs, True)

# The tool check must name tRNAscan-SE with exact case: macOS filesystems are
# case-insensitive, so a lowercase probe would pass there and mask a real gap.
check("tool check probes tRNAscan-SE with exact case",
      "tRNAscan-SE" in workflow_src, True)

# A regression test only protects anything if something runs it. CI used to
# run no tests at all — only the single-genome smoke check.
check("CI runs the test suite", "tests/test_*.py" in workflow_src, True)

# The default tRNA mode is far slower than the one the last green run tested;
# an unbounded job would burn the 6 h default before failing.
check("the job is time-bounded", "timeout-minutes:" in workflow_src, True)

print()
print("%d checks, %d failed" % (RUN[0], len(FAIL)))
for f in FAIL:
    print("  FAILED: %s" % f)
sys.exit(1 if FAIL else 0)
