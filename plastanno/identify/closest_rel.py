"""
Find closest relatives using BLAST vs genus representatives.
Returns top-N relatives with taxonomy info.

Improvements over v1:
- Taxonomy info included in return value
- Cleaner interface
"""
import os
import subprocess
import tempfile
from pathlib import Path
from typing import List, Dict

import pandas as pd


TRAIN_CSV = str(
    Path(__file__).parent.parent.parent /
    "database" / "splits" / "train.csv"
)

# Load taxonomy once at import
_TAX_DICT = None

def _load_taxonomy():
    global _TAX_DICT
    if _TAX_DICT is not None:
        return _TAX_DICT
    train_csv = Path(TRAIN_CSV)
    if not train_csv.exists():
        _TAX_DICT = {}
        return _TAX_DICT
    df = pd.read_csv(train_csv)[
        ["accession","genus","family","order"]
    ]
    _TAX_DICT = {
        row["accession"].split(".")[0]: {
            "genus" : row["genus"],
            "family": row["family"],
            "order" : row["order"],
        }
        for _, row in df.iterrows()
    }
    return _TAX_DICT


def find_closest_relatives(
    query_fasta : str,
    blast_db    : str,
    gene_catalog: dict,
    top_n       : int = 10,
    threads     : int = 4,
) -> Dict:
    """
    BLAST query genome against genus representatives.

    Returns:
        {
          "relatives": [
            {
              "accession": "NC_XXXXXX",
              "pident"   : 99.8,
              "qcovs"    : 99.0,
              "bitscore" : 66000.0,
              "genus"    : "Gynostemma",
              "family"   : "Cucurbitaceae",
              "order"    : "Cucurbitales",
            },
            ...
          ]
        }
    """
    result = subprocess.run([
        "blastn",
        "-query",        query_fasta,
        "-db",           blast_db,
        "-outfmt",
        "6 sseqid pident qcovs bitscore",
        "-max_target_seqs", str(top_n * 2),
        "-num_threads",  str(threads),
    ], capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(
            f"BLASTN failed: {result.stderr}"
        )

    # Parse and deduplicate
    seen      = set()
    relatives = []
    for line in result.stdout.strip().split("\n"):
        if not line: continue
        p = line.split("\t")
        if len(p) < 4: continue

        acc      = p[0].split(".")[0]
        pident   = float(p[1])
        qcovs    = float(p[2])
        bitscore = float(p[3])

        if acc in seen: continue
        seen.add(acc)

        rel = {
            "accession": acc,
            "pident"   : pident,
            "qcovs"    : qcovs,
            "bitscore" : bitscore,
        }

        # Add taxonomy
        tax = _load_taxonomy().get(acc, {})
        rel["genus"]  = tax.get("genus",  "")
        rel["family"] = tax.get("family", "")
        rel["order"]  = tax.get("order",  "")

        relatives.append(rel)
        if len(relatives) >= top_n:
            break

    return {"relatives": relatives}
