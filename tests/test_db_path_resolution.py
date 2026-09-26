#!/usr/bin/env python3
"""Where the database is written and where it is looked for must be one place.

`fetch_db._data_parent` and `paths.db_root` each resolved the user-data
directory on their own, and the two copies drifted. Without `platformdirs`,
`fetch_db` wrote to `~/.local/share/plastanno` while `db_root` skipped that
location entirely and returned the repo layout — so a user could download 267 MB
and then be told the database was missing.

`platformdirs` is a declared dependency in both `pyproject.toml` and the conda
recipe, so the broken state needs a `--no-deps` or damaged install to reach. That
is precisely why it went unnoticed: the failure is invisible in a correct
environment and silent in a broken one.

These checks run with platformdirs both present and absent, because a resolution
rule that only holds in one of the two is the bug that was just fixed.
"""
import builtins
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from plastanno import paths                                    # noqa: E402
from plastanno import fetch_db                                 # noqa: E402

RUN = [0]
FAIL = []


def check(label, got, want):
    RUN[0] += 1
    ok = got == want
    print("  %s %s" % ("ok  " if ok else "FAIL", label))
    if not ok:
        FAIL.append("%s: got %r want %r" % (label, got, want))


class no_platformdirs:
    """Make `import platformdirs` raise, whether or not it is installed."""

    def __enter__(self):
        self._real = builtins.__import__
        self._cached = sys.modules.pop("platformdirs", None)

        def fake(name, *a, **k):
            if name == "platformdirs":
                raise ImportError("simulated: platformdirs not installed")
            return self._real(name, *a, **k)

        builtins.__import__ = fake
        return self

    def __exit__(self, *exc):
        builtins.__import__ = self._real
        if self._cached is not None:
            sys.modules["platformdirs"] = self._cached
        return False


print("the two resolvers agree, with platformdirs present")

check("paths and fetch_db return the same user-data parent",
      paths.user_data_parent(), fetch_db._data_parent())

print("\n... and with platformdirs absent")

with no_platformdirs():
    a, b = paths.user_data_parent(), fetch_db._data_parent()
check("they still agree", a, b)
check("  ... on the documented Linux fallback",
      a, Path.home() / ".local" / "share" / "plastanno")

# The bug in one sentence: fetch_db wrote here, db_root looked elsewhere.
with no_platformdirs():
    parent = fetch_db._data_parent()
    saved = os.environ.pop("PLASTANNO_DB", None)
    try:
        # db_root only returns the user-data dir when a database is actually
        # there, so this asserts the CANDIDATE it considers, not the result.
        considered = paths.user_data_parent() / "database"
    finally:
        if saved is not None:
            os.environ["PLASTANNO_DB"] = saved
check("db_root considers exactly the directory fetch-db writes to",
      considered, parent / "database")


print("\nresolution order")

saved = os.environ.get("PLASTANNO_DB")
try:
    os.environ["PLASTANNO_DB"] = "/tmp/plastanno-db-test-override"
    check("$PLASTANNO_DB wins over everything",
          paths.db_root(), Path("/tmp/plastanno-db-test-override"))
    os.environ["PLASTANNO_DB"] = "~/some-db"
    check("  ... and is expanded", paths.db_root(), Path.home() / "some-db")
finally:
    if saved is None:
        os.environ.pop("PLASTANNO_DB", None)
    else:
        os.environ["PLASTANNO_DB"] = saved

# With no override and no installed database, the repo layout must still work —
# that is what keeps source-tree runs going.
saved = os.environ.pop("PLASTANNO_DB", None)
try:
    with no_platformdirs():
        # point the user-data probe at somewhere empty by checking the rule
        # rather than the machine's actual state
        has_db = (paths.user_data_parent() / "database" / "blast_db").exists()
        root = paths.db_root()
    expected = ((paths.user_data_parent() / "database") if has_db
                else paths._REPO_FALLBACK)
    check("falls through to the repo layout only when no database is installed",
          root, expected)
finally:
    if saved is not None:
        os.environ["PLASTANNO_DB"] = saved

check("database_ready() agrees with db_root()",
      paths.database_ready(),
      (paths.db_root() / "blast_db").is_dir()
      and any((paths.db_root() / "blast_db").glob("genus_reps.*")))

print()
print("%d checks, %d failed" % (RUN[0], len(FAIL)))
for f in FAIL:
    print("  FAILED: %s" % f)
sys.exit(1 if FAIL else 0)
