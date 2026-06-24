#!/usr/bin/env python3
"""
Partition the frozen held-out set into N stratified-by-FAMILY chunks for batched
final evaluation.

Creating this partition does NOT touch annotation results — it only groups
accessions by taxonomic family read from the GenBank ORGANISM lineage — so it is
safe to run before the tool is frozen. The chunks must only be EVALUATED once,
with the tool frozen:

    python3 scripts/benchmark/multi_genome_bench.py --set heldout --final \
        --acc-file splits/heldout_chunks/chunk_00.txt --workers 16 --keep
    ... repeat for chunk_01 .. chunk_09, then pool the 10 results.

Stratification: within each family, accessions are distributed round-robin across
the N chunks after a seeded shuffle (with a per-family offset so even singleton
families spread out). Every chunk is therefore a family-balanced ~1/N sample of
the whole held-out set — NOT k-fold cross-validation (the tool is not retrained;
the chunks are disjoint test subsamples used for compute batching and a
batch-to-batch stability estimate).

Outputs (with --write):
    splits/heldout_chunks/chunk_00.txt .. chunk_<N-1>.txt   one accession/line
    splits/heldout_chunks/manifest.json   seed, per-chunk size + sha256 + families
    splits/heldout_families.tsv           cached accession -> family
"""
import os, re, json, hashlib, random, argparse
from collections import defaultdict, Counter

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("PLASTANNO_DATA","benchmark_data")+"/rawdata"
HELDOUT = os.path.join(REPO, "splits", "heldout_set.txt")
OUTDIR = os.path.join(REPO, "splits", "heldout_chunks")
FAMCACHE = os.path.join(REPO, "splits", "heldout_families.tsv")
SEED = 20260621
NCHUNKS = 10

strip = lambda a: re.sub(r"\.\d+$", "", a)


def _family_from_gb(path):
    from Bio import SeqIO
    try:
        rec = next(SeqIO.parse(path, "genbank"))
    except Exception:
        return "UNKNOWN"
    tax = rec.annotations.get("taxonomy", [])
    cand = [t for t in tax if t.endswith("aceae")]   # family rank
    return cand[-1] if cand else (tax[-1] if tax else "UNKNOWN")


def load_families(accs, rawmap):
    if os.path.exists(FAMCACHE):
        fam = {}
        for line in open(FAMCACHE):
            a, f = line.rstrip("\n").split("\t")
            fam[a] = f
        if all(a in fam for a in accs):
            return fam
    fam = {}
    for i, acc in enumerate(accs):
        fn = rawmap.get(strip(acc))
        fam[acc] = _family_from_gb(os.path.join(RAW, fn + ".gb")) if fn else "UNKNOWN"
        if (i + 1) % 500 == 0:
            print("  ... families read %d/%d" % (i + 1, len(accs)))
    return fam


def main(write):
    held = [l.strip() for l in open(HELDOUT) if l.strip()]
    rawmap = {}
    for fn in os.listdir(RAW):
        if fn.endswith(".gb"):
            rawmap[strip(fn[:-3])] = fn[:-3]

    fam = load_families(held, rawmap)
    if write and not os.path.exists(FAMCACHE):
        with open(FAMCACHE, "w") as f:
            for a in held:
                f.write("%s\t%s\n" % (a, fam[a]))

    by_fam = defaultdict(list)
    for a in held:
        by_fam[fam.get(a, "UNKNOWN")].append(a)

    rng = random.Random(SEED)
    chunks = [[] for _ in range(NCHUNKS)]
    for i, famname in enumerate(sorted(by_fam)):
        members = by_fam[famname][:]
        rng.shuffle(members)
        for j, acc in enumerate(members):
            chunks[(i + j) % NCHUNKS].append(acc)   # offset i spreads small families

    # integrity
    assert sum(len(c) for c in chunks) == len(held), "lost accessions"
    flat = [a for c in chunks for a in c]
    assert len(set(flat)) == len(held), "overlap/duplication between chunks"

    print("Held-out: %d accessions across %d families -> %d chunks (seed=%d)"
          % (len(held), len(by_fam), NCHUNKS, SEED))
    for c, ch in enumerate(chunks):
        print("  chunk_%02d: %4d accessions, %3d families"
              % (c, len(ch), len(set(fam[a] for a in ch))))

    if not write:
        print("\n(dry-run) pass --write to create chunk files + manifest")
        return

    os.makedirs(OUTDIR, exist_ok=True)
    manifest = {"seed": SEED, "nchunks": NCHUNKS, "total": len(held),
                "n_families": len(by_fam), "stratified_by": "family",
                "chunks": []}
    for c, ch in enumerate(chunks):
        ch_sorted = sorted(ch)
        body = "\n".join(ch_sorted) + "\n"
        open(os.path.join(OUTDIR, "chunk_%02d.txt" % c), "w").write(body)
        manifest["chunks"].append({
            "chunk": c, "n": len(ch_sorted),
            "sha256": hashlib.sha256(body.encode()).hexdigest(),
            "top_families": dict(Counter(fam[a] for a in ch).most_common(5)),
        })
    json.dump(manifest, open(os.path.join(OUTDIR, "manifest.json"), "w"), indent=2)
    print("\nWrote %d chunks + manifest.json to %s" % (NCHUNKS, OUTDIR))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true",
                    help="create chunk files + manifest (otherwise dry-run)")
    main(ap.parse_args().write)
