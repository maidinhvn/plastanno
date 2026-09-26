# Tests shipped with Plastanno

These cover the annotation tool itself and need no external dataset:

    for t in tests/test_*.py; do python3 "$t" || echo "FAILED: $t"; done

Each prints `N checks, M failed` and exits non-zero on failure.

The benchmark scoring machinery — and its own tests — is deliberately **not** in
this repository. It scores predictions against reference GenBank files from a
dataset that is not distributed, so those tests would only fail with ImportError
here. They live with the evaluation record, not with the tool.
