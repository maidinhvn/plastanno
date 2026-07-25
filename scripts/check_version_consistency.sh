#!/usr/bin/env bash
# Release guard: pyproject.toml and plastanno/__init__.py both declare the version
# (pyproject must stay static because the root plastanno.py shim shadows the
# package, which breaks setuptools' dynamic `attr:` lookup at build time).
# Run this before tagging a release.
set -u
cd "$(dirname "$0")/.." || exit 1

PY=$(grep -m1 '^version *= *"' pyproject.toml        | sed 's/.*"\(.*\)".*/\1/')
IN=$(grep -m1 '^__version__ *= *"' plastanno/__init__.py | sed 's/.*"\(.*\)".*/\1/')

echo "  pyproject.toml        : ${PY:-<missing>}"
echo "  plastanno/__init__.py : ${IN:-<missing>}"

if [ -n "$PY" ] && [ "$PY" = "$IN" ]; then
    echo "  OK — versions agree ($PY)"
    exit 0
fi
echo "  MISMATCH — bump both before tagging."
exit 1
