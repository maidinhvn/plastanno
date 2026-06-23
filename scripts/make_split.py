#!/usr/bin/env python3
"""
Stratified 70/30 DEV/HELD-OUT split of the plastome catalogue, by structure
`mode`. The held-out set is written ONCE and must never be touched during
development. NC_053537 (used throughout development) is forced into DEV.

Run with --write only after inspecting the printed class distribution.
"""
import json, os, re, random, hashlib, argparse, datetime
from collections import Counter, defaultdict

CLS = os.environ.get("PLASTANNO_DATA","/data06/users/vutrinh/Apiales_Plastomes_20260514")+"/classify_results.json"
RAW = os.environ.get("PLASTANNO_DATA","/data06/users/vutrinh/Apiales_Plastomes_20260514")+"/rawdata"
OUT = "splits"
SEED = 20260619
HELDOUT_FRAC = 0.30
FORCE_DEV = {"NC_053537"}

strip = lambda a: re.sub(r"\.\d+$", "", a)

def main(write):
    cls = json.load(open(CLS))
    mode_of = {strip(k): v.get("mode", "UNKNOWN") for k, v in cls.items()}

    # accession (versioned) present as a raw .gb pair, keyed by version-less id
    raw_versioned = {}
    for fn in os.listdir(RAW):
        if fn.endswith(".gb"):
            acc = fn[:-3]                       # e.g. NC_000927.1
            raw_versioned.setdefault(strip(acc), acc)

    ids = sorted(raw_versioned)                 # version-less ids that have a .gb
    by_mode = defaultdict(list)
    for vid in ids:
        by_mode[mode_of.get(vid, "UNKNOWN")].append(vid)

    print("Tổng genome (raw, có .gb): %d" % len(ids))
    print("Phân bố theo mode:")
    for m in sorted(by_mode, key=lambda x: -len(by_mode[x])):
        n = len(by_mode[m]); print("  %-12s %5d  -> heldout~%d dev~%d"
                                   % (m, n, round(n*HELDOUT_FRAC), n-round(n*HELDOUT_FRAC)))

    rng = random.Random(SEED)
    dev, held = [], []
    for m in sorted(by_mode):
        grp = by_mode[m][:]
        rng.shuffle(grp)
        k = round(len(grp) * HELDOUT_FRAC)
        h = [x for x in grp if x not in FORCE_DEV][:k]   # never put forced ids in held-out
        hs = set(h)
        held += h
        dev  += [x for x in grp if x not in hs]

    assert FORCE_DEV & set(dev) == FORCE_DEV, "forced id not in DEV!"
    assert not (set(dev) & set(held)), "overlap between DEV and HELDOUT!"
    assert len(dev) + len(held) == len(ids)

    print("\nDEV=%d  HELDOUT=%d  (held-out %.1f%%)"
          % (len(dev), len(held), 100*len(held)/len(ids)))
    print("HELDOUT theo mode:", dict(Counter(mode_of.get(x,"UNKNOWN") for x in held)))

    if not write:
        print("\n[DRY-RUN] chưa ghi file. Chạy lại với --write để cố định.")
        return

    # write versioned accessions (the actual filename stems) for direct use
    def dump(path, vids):
        with open(path, "w") as fh:
            for vid in sorted(vids):
                fh.write(raw_versioned[vid] + "\n")
    dump(os.path.join(OUT, "dev_set.txt"), dev)
    dump(os.path.join(OUT, "heldout_set.txt"), held)

    man = {
        "created": datetime.datetime.now().isoformat(timespec="seconds"),
        "seed": SEED, "heldout_frac": HELDOUT_FRAC,
        "n_total": len(ids), "n_dev": len(dev), "n_heldout": len(held),
        "forced_dev": sorted(FORCE_DEV),
        "heldout_by_mode": dict(Counter(mode_of.get(x,"UNKNOWN") for x in held)),
        "dev_by_mode": dict(Counter(mode_of.get(x,"UNKNOWN") for x in dev)),
        "dev_sha256": hashlib.sha256("\n".join(sorted(dev)).encode()).hexdigest(),
        "heldout_sha256": hashlib.sha256("\n".join(sorted(held)).encode()).hexdigest(),
    }
    json.dump(man, open(os.path.join(OUT, "split_manifest.json"), "w"), indent=2)
    print("\nĐã ghi splits/dev_set.txt, splits/heldout_set.txt, splits/split_manifest.json")
    print("HELD-OUT cố định vĩnh viễn — KHÔNG chạm trong lúc phát triển.")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true")
    main(ap.parse_args().write)
