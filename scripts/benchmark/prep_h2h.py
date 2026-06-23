#!/usr/bin/env python3
"""Prepare PGA head-to-head: stratified subset of held-out v2 + fair DEV references."""
import os, re, random
import pandas as pd
norm=lambda a:re.sub(r'\.\d+$','',a.strip())
REPO="/data06/users/vutrinh/Plastanno_v2"
TSV="/data06/users/vutrinh/PlastAnnot/database/all_genes_full.tsv"
N=250; random.seed(20260622)

# DEV taxonomy: genus/family/order -> [dev versioned accs]
dev=set(l.strip() for l in open(f"{REPO}/splits/dev_set.txt") if l.strip())
devn={norm(a) for a in dev}
gen2dev={}; fam2dev={}; ord2dev={}
seen=set()
with open(TSV) as f:
    h=f.readline().rstrip("\n").split("\t"); ix={c:h.index(c) for c in ("accession","genus","family","order","genome_len")}
    for line in f:
        p=line.rstrip("\n").split("\t")
        if len(p)<=ix["genome_len"]: continue
        a=norm(p[ix["accession"]])
        if a not in devn or a in seen: continue
        seen.add(a)
        ver=p[ix["accession"]]
        gen2dev.setdefault(p[ix["genus"]],[]).append(ver)
        fam2dev.setdefault(p[ix["family"]],[]).append(ver)
        ord2dev.setdefault(p[ix["order"]],[]).append(ver)
print(f"DEV refs indexed: {len(seen)} genomes")

# test taxonomy from inventory xlsx
df=pd.read_excel(f"{REPO}/docs/Plastanno_dataset_inventory.xlsx", sheet_name="Test_set (heldout v2 2151)")
df["acc"]=df["accession"]

# stratified subset: proportional by family, both strata
by_fam={}
for _,r in df.iterrows(): by_fam.setdefault(str(r["family"]),[]).append(r)
subset=[]
fams=sorted(by_fam, key=lambda k:-len(by_fam[k]))
# proportional sampling
total=len(df)
for fam in fams:
    rows=by_fam[fam]; k=max(1, round(len(rows)/total*N))
    random.shuffle(rows); subset.extend(rows[:k])
random.shuffle(subset); subset=subset[:N]
print(f"subset: {len(subset)} genomes, {len(set(str(r['family']) for r in subset))} families")
from collections import Counter
print("  per stratum:", dict(Counter(r["set"] for r in subset)))

# fair reference per target: same genus -> family -> order (exclude self/heldout: dev only)
refmap={}; refpool=set()
miss=0
for r in subset:
    tacc=r["acc"]; g=str(r["genus"]); fa=str(r["family"]); od=str(r["order"])
    cands=[a for a in gen2dev.get(g,[]) if norm(a)!=norm(tacc)][:2]
    if not cands: cands=[a for a in fam2dev.get(fa,[]) if norm(a)!=norm(tacc)][:2]
    if not cands: cands=[a for a in ord2dev.get(od,[]) if norm(a)!=norm(tacc)][:2]
    if not cands: miss+=1; continue
    refmap[tacc]=cands; refpool.update(cands)
print(f"targets with a fair DEV reference: {len(refmap)}/{len(subset)} (missing {miss})")
print(f"reference pool (distinct DEV genomes): {len(refpool)}")

os.makedirs(f"{REPO}/bench_runs/h2h", exist_ok=True)
open(f"{REPO}/splits/h2h_subset.txt","w").write("\n".join(refmap.keys())+"\n")
with open(f"{REPO}/splits/h2h_refmap.tsv","w") as o:
    for t,cs in refmap.items(): o.write(f"{t}\t{','.join(cs)}\n")
open(f"{REPO}/splits/h2h_refpool.txt","w").write("\n".join(sorted(refpool))+"\n")
print("WROTE splits/h2h_subset.txt, h2h_refmap.tsv, h2h_refpool.txt")
