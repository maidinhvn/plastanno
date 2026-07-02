#!/usr/bin/env python3
"""
Plastanno v2 — Plastome Annotation Tool (command-line interface)
Hybrid: Reference-based (Engine A) + Model-based (Engine B)

Usage:
    plastanno run <input.fasta> [options]
    plastanno batch <input_dir> [options]
    plastanno fetch-db [--force]

Examples:
    plastanno run genome.fasta
    plastanno run genome.fasta --output results/ --threads 8
    plastanno batch genomes/ --output results/ --threads 8
"""
import argparse
import sys
import time
from pathlib import Path


def _require_database():
    """Preflight: make sure the reference database is present before annotating.

    A fresh `conda install plastanno` / `pip install` ships only the code, not the
    ~266 MB database, so without this guard the first `plastanno run` dies deep in
    the pipeline with an opaque BLAST error. Fail early with the exact fix instead.
    """
    from plastanno.paths import database_ready, db_root
    if not database_ready():
        print("ERROR: reference database not found (looked in: %s)" % db_root())
        print()
        print("Plastanno needs its ~266 MB reference database before it can annotate.")
        print("Download it once with:")
        print()
        print("    plastanno fetch-db")
        print()
        print("(or set $PLASTANNO_DB to a directory that already holds it).")
        sys.exit(1)


def cmd_run(args):
    """Annotate a single plastome."""
    from plastanno.pipeline import run

    if not Path(args.input).exists():
        print(f"ERROR: Input file not found: {args.input}")
        sys.exit(1)
    _require_database()

    prefix = args.prefix or Path(args.input).stem
    output = args.output or "plastanno_output"

    print("=" * 60)
    print("  Plastanno v2.0 — Plastome Annotation Tool")
    print("=" * 60)
    print(f"  Input  : {args.input}")
    print(f"  Output : {output}")
    print(f"  Threads: {args.threads}")
    print("=" * 60)

    result = run(
        input_fasta = args.input,
        output_dir  = output,
        prefix      = prefix,
        threads     = args.threads,
        no_plot     = args.no_plot,
        reference   = args.reference,
        use_trnascan= args.trnascan,
    )
    return result


def cmd_batch(args):
    """Annotate multiple plastomes in a directory."""
    input_dir = Path(args.input_dir)
    fastas    = sorted(input_dir.glob("*.fasta")) + \
                sorted(input_dir.glob("*.fa"))

    if not fastas:
        print(f"ERROR: No FASTA files in {input_dir}")
        sys.exit(1)
    _require_database()

    print(f"Found {len(fastas)} FASTA files")
    output = Path(args.output or "plastanno_batch_output")
    output.mkdir(parents=True, exist_ok=True)

    results = []
    for i, fa in enumerate(fastas, 1):
        print(f"\n[{i}/{len(fastas)}] {fa.name}")
        try:
            from plastanno.pipeline import run
            result = run(
                input_fasta = str(fa),
                output_dir  = str(output / fa.stem),
                threads     = args.threads,
                no_plot     = args.no_plot,
                reference   = args.reference,
                use_trnascan= args.trnascan,
            )
            results.append({
                "file"   : fa.name,
                "genes"  : len(result["annotations"]),
                "elapsed": result["elapsed"],
                "status" : "OK",
            })
        except Exception as e:
            print(f"  ERROR: {e}")
            results.append({
                "file"  : fa.name,
                "status": f"FAILED: {e}",
            })

    # Summary
    print(f"\n{'='*60}")
    print(f"Batch complete: {len(fastas)} genomes")
    ok   = sum(1 for r in results if r["status"]=="OK")
    fail = len(results) - ok
    print(f"  Success: {ok}  Failed: {fail}")
    if ok > 0:
        avg_genes = sum(
            r["genes"] for r in results if r["status"]=="OK"
        ) / ok
        avg_time  = sum(
            r["elapsed"] for r in results if r["status"]=="OK"
        ) / ok
        print(f"  Avg genes  : {avg_genes:.0f}")
        print(f"  Avg runtime: {avg_time:.1f}s")
    print(f"{'='*60}")


def cmd_fetch_db(args):
    """Download and install the reference database into the platform data dir."""
    from plastanno.fetch_db import fetch_db
    fetch_db(force=args.force)


def main():
    parser = argparse.ArgumentParser(
        description="Plastanno v2 — Plastome Annotation Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    subparsers = parser.add_subparsers(dest="command")

    # run command
    run_parser = subparsers.add_parser(
        "run", help="Annotate a single plastome"
    )
    run_parser.add_argument("input",
        help="Input FASTA file")
    run_parser.add_argument("--output", "-o",
        help="Output directory (default: plastanno_output)")
    run_parser.add_argument("--prefix", "-p",
        help="Output file prefix (default: input filename)")
    run_parser.add_argument("--threads", "-t",
        type=int, default=4,
        help="Number of threads (default: 4)")
    run_parser.add_argument("--no-plot",
        action="store_true",
        help="Skip circular plot generation")
    run_parser.add_argument("--reference",
        help="Annotated reference plastome (GenBank) whose per-gene proteins are "
             "overlaid on the built-in DB for Engine A (automatic fallback to DB)")
    run_parser.add_argument("--trnascan",
        action="store_true",
        help="Also run tRNAscan-SE as an extra tRNA source (must be on PATH)")

    # batch command
    batch_parser = subparsers.add_parser(
        "batch", help="Annotate multiple plastomes"
    )
    batch_parser.add_argument("input_dir",
        help="Directory containing FASTA files")
    batch_parser.add_argument("--output", "-o",
        help="Output directory")
    batch_parser.add_argument("--threads", "-t",
        type=int, default=4)
    batch_parser.add_argument("--no-plot",
        action="store_true")
    batch_parser.add_argument("--reference",
        help="Annotated reference plastome (GenBank) overlaid on Engine A's DB")
    batch_parser.add_argument("--trnascan",
        action="store_true",
        help="Also run tRNAscan-SE as an extra tRNA source")

    # fetch-db command
    fetch_parser = subparsers.add_parser(
        "fetch-db", help="Download the reference database (~266 MB) from Zenodo"
    )
    fetch_parser.add_argument("--force",
        action="store_true",
        help="Re-download even if a database is already present")

    args = parser.parse_args()

    if args.command == "run":
        cmd_run(args)
    elif args.command == "batch":
        cmd_batch(args)
    elif args.command == "fetch-db":
        cmd_fetch_db(args)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
