"""Central resolution of the reference-database location.

Resolution order (first hit wins):
  1. ``$PLASTANNO_DB``                         — explicit override
  2. platform user-data dir (if it has a DB)   — populated by ``plastanno fetch-db``
     e.g. ``~/.local/share/plastanno/database`` (Linux),
          ``~/Library/Application Support/plastanno/database`` (macOS)
  3. repo-layout fallback ``<repo>/database``   — keeps source-tree / dev runs working

Behaviour is unchanged when running from the source tree with ``$PLASTANNO_DB``
unset and no installed data dir: it falls through to (3), the historical location.
"""
import os
from pathlib import Path

_REPO_FALLBACK = Path(__file__).resolve().parent.parent / "database"
_PKG_DATA      = Path(__file__).resolve().parent / "data"


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
    try:
        import platformdirs
        cand = Path(platformdirs.user_data_dir("plastanno")) / "database"
        if (cand / "blast_db").exists():
            return cand
    except Exception:
        pass  # platformdirs absent or unusable -> fall back to repo layout
    return _REPO_FALLBACK
