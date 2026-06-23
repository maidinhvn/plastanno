#!/usr/bin/env bash
# Download and install the Plastanno reference databases (~266 MB) from Zenodo
# into ./database/. Run once after cloning:  bash scripts/get_database.sh
set -e
cd "$(dirname "$0")/.."                 # project root

URL="https://zenodo.org/records/20807995/files/plastanno-database.tar.gz"
MD5="d366f7ea58a78a5cb43cb42a8639d2e2"
TAR="plastanno-database.tar.gz"

echo "[1/3] Downloading reference databases (~266 MB) from Zenodo ..."
if command -v curl >/dev/null 2>&1; then
    curl -L -o "$TAR" "$URL"
elif command -v wget >/dev/null 2>&1; then
    wget -O "$TAR" "$URL"
else
    echo "ERROR: need 'curl' or 'wget' on PATH." >&2; exit 1
fi

echo "[2/3] Verifying checksum ..."
if command -v md5sum >/dev/null 2>&1; then
    echo "$MD5  $TAR" | md5sum -c - || { echo "ERROR: checksum mismatch." >&2; exit 1; }
else
    echo "  (md5sum not found — skipping checksum)"
fi

echo "[3/3] Unpacking into database/ ..."
tar xzf "$TAR"
rm -f "$TAR"
echo "Done. Reference databases are installed in database/  (verify with: bash diagnose.sh)"
