"""Download and install the Plastanno reference database from Zenodo.

The database (~266 MB compressed, ~700 MB extracted) is not shipped with the
package. ``plastanno fetch-db`` downloads it into the platform user-data dir, where
:func:`plastanno.paths.db_root` then finds it automatically. Override the location
with ``$PLASTANNO_DB``.
"""
import hashlib
import sys
import tarfile
import urllib.request
from pathlib import Path

# Zenodo concept DOI 10.5281/zenodo.20807994 -> version record 20807995
URL = "https://zenodo.org/records/20807995/files/plastanno-database.tar.gz"
MD5 = "d366f7ea58a78a5cb43cb42a8639d2e2"


def _data_parent() -> Path:
    """Directory whose ``database/`` subdir db_root() resolves to."""
    try:
        import platformdirs
        return Path(platformdirs.user_data_dir("plastanno"))
    except Exception:
        return Path.home() / ".local" / "share" / "plastanno"


def _progress(blocks, bsize, total):
    if total > 0:
        pct = min(100, blocks * bsize * 100 // total)
        sys.stdout.write(f"\r      {pct}%")
        sys.stdout.flush()


def fetch_db(force: bool = False) -> Path:
    parent = _data_parent()
    dbdir = parent / "database"
    if (dbdir / "blast_db").exists() and not force:
        print(f"Database already present at {dbdir} (use --force to re-download).")
        return dbdir

    parent.mkdir(parents=True, exist_ok=True)
    tar = parent / "plastanno-database.tar.gz"

    print(f"[1/3] Downloading reference database (~266 MB) from Zenodo ...")
    urllib.request.urlretrieve(URL, tar, _progress)
    print()

    print("[2/3] Verifying checksum ...")
    h = hashlib.md5()
    with open(tar, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    if h.hexdigest() != MD5:
        tar.unlink()
        sys.exit(f"ERROR: checksum mismatch ({h.hexdigest()} != {MD5})")

    print(f"[3/3] Extracting into {parent} ...")
    with tarfile.open(tar) as t:
        try:
            t.extractall(parent, filter="data")   # Python >= 3.12
        except TypeError:
            t.extractall(parent)
    tar.unlink()
    print(f"Done. Database installed at {dbdir}")
    return dbdir
