"""Plastanno — dual-evidence, confidence-aware plastome annotation.

The version appears twice: here and in ``pyproject.toml``. A dynamic
``attr: plastanno.__version__`` would remove the duplication, but this repository
ships both a ``plastanno/`` package and a ``plastanno.py`` entry script, and
setuptools resolves the name to the script — the build then fails with
"'plastanno' is not a package". So the two literals stay, and
``tests/test_version_consistency.py`` asserts they agree; a stale egg-info once
made ``--version`` report an old release while the code was newer.
"""

__version__ = "3.0.0"
