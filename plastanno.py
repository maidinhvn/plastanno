#!/usr/bin/env python3
"""Backward-compatible entry point.

The CLI now lives in ``plastanno/cli.py`` (installed as the ``plastanno`` command).
This thin shim keeps ``python3 plastanno.py run ...`` working from the source tree.
"""
import os
import sys

# Ensure the source-tree package is importable when run as a loose script.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from plastanno.cli import main

if __name__ == "__main__":
    main()
