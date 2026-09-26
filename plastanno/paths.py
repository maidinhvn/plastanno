"""Central resolution of the reference-database location.

Resolution order (first hit wins):
  1. ``$PLASTANNO_DB``                         — explicit override
  2. platform user-data dir (if it has a DB)   — populated by ``plastanno fetch-db``
     e.g. ``~/.local/share/plastanno/database`` (Linux),
          ``~/Library/Application Support/plastanno/database`` (macOS)
  3. repo-layout fallback ``<repo>/database``   — keeps source-tree / dev runs working

The user-data directory is resolved by :func:`user_data_parent`, which
``fetch_db`` imports as well so that the place the database is written to and the
place it is looked for cannot drift apart.

Behaviour is unchanged when running from the source tree with ``$PLASTANNO_DB``
unset and no installed data dir: it falls through to (3), the historical location.
"""
import os
from pathlib import Path

_REPO_FALLBACK = Path(__file__).resolve().parent.parent / "database"
_PKG_DATA      = Path(__file__).resolve().parent / "data"


def user_data_parent() -> Path:
    """Directory whose ``database/`` subdir ``plastanno fetch-db`` populates.

    Single source of truth, imported by ``fetch_db`` as well. It used to be
    duplicated there, and the two copies drifted: without ``platformdirs``,
    ``fetch_db`` wrote to ``~/.local/share/plastanno`` while ``db_root`` skipped
    that location entirely and fell through to the repo layout. A user in that
    state downloaded 267 MB and was then told the database was missing.
    """
    try:
        import platformdirs
        return Path(platformdirs.user_data_dir("plastanno"))
    except Exception:
        # platformdirs is a declared dependency, so this branch is for broken or
        # --no-deps installs. It must match what fetch_db would use, or the two
        # disagree again.
        return Path.home() / ".local" / "share" / "plastanno"


def config_dir() -> Path:
    """Small runtime configs (gene_catalog.json, exon_templates.json, boundary_db/).

    These ship inside the wheel (``plastanno/data``); fall back to the big-DB root
    for source-tree runs where they may live under ``database/`` instead.
    """
    if (_PKG_DATA / "gene_catalog.json").exists():
        return _PKG_DATA
    return db_root()


def db_root() -> Path:
    """Return the reference-database root directory."""
    env = os.environ.get("PLASTANNO_DB")
    if env:
        return Path(env).expanduser()
    cand = user_data_parent() / "database"
    if (cand / "blast_db").exists():
        return cand
    return _REPO_FALLBACK


def database_ready() -> bool:
    """True if the reference database has been fetched at the resolved location.

    Checks the closest-relative BLAST db (`blast_db/genus_reps.*`), which is the
    first thing the pipeline needs — so a friendly preflight message can replace
    the cryptic BLAST error a fresh install would otherwise hit."""
    bd = db_root() / "blast_db"
    return bd.is_dir() and any(bd.glob("genus_reps.*"))
