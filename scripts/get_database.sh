#!/usr/bin/env bash
# Download and install the Plastanno reference databases (~266 MB) from Zenodo
# into ./database/. Run once after cloning:  bash scripts/get_database.sh
set -e
cd "$(dirname "$0")/.."                 # project root

# The download location is defined in exactly one place — plastanno/fetch_db.py,
# which `plastanno fetch-db` also uses. Read it from there rather than repeating
# it here: the two copies drifted once already (this script kept pointing at a
# superseded Zenodo record, whose checksum still matched, so CI silently
# installed an older database instead of failing).
SRC="plastanno/fetch_db.py"
URL=$(sed -n 's/^URL[[:space:]]*=[[:space:]]*"\(.*\)".*/\1/p' "$SRC")
MD5=$(sed -n 's/^MD5[[:space:]]*=[[:space:]]*"\(.*\)".*/\1/p' "$SRC")
if [ -z "$URL" ] || [ -z "$MD5" ]; then
    echo "ERROR: could not read URL/MD5 from $SRC." >&2
    echo "       Expected lines of the form:  URL = \"https://...\"" >&2
    exit 1
fi
TAR="plastanno-database.tar.gz"

echo "[1/3] Downloading reference databases (~266 MB) from Zenodo ..."
echo "      $URL"
if command -v curl >/dev/null 2>&1; then
    curl -L -o "$TAR" "$URL"
elif command -v wget >/dev/null 2>&1; then
    wget -O "$TAR" "$URL"
else
    echo "ERROR: need 'curl' or 'wget' on PATH." >&2; exit 1
fi

echo "[2/3] Verifying checksum ..."
if command -v md5sum >/dev/null 2>&1; then            # Linux
    echo "$MD5  $TAR" | md5sum -c - || { echo "ERROR: checksum mismatch." >&2; exit 1; }
elif command -v md5 >/dev/null 2>&1; then             # macOS
    [ "$(md5 -q "$TAR")" = "$MD5" ] || { echo "ERROR: checksum mismatch." >&2; exit 1; }
else
    echo "  (no md5 tool found — skipping checksum)"
fi

echo "[3/3] Unpacking into database/ ..."
tar xzf "$TAR"
rm -f "$TAR"
echo "Done. Reference databases are installed in database/  (verify with: bash diagnose.sh)"
