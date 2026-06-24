#!/usr/bin/env python3
"""Plastanno vs PGA on the head-to-head subset. Both tools are scored here, in
parallel, with the identical gene-by-gene scorer (+-60 bp ends + 0.6 similarity).
Predictions come from the self-contained head-to-head run produced by run_h2h.sh:
Plastanno -> bench_runs/h2h/plast_out/<acc>/<acc>.gb, PGA -> bench_runs/h2h/pga_out/<acc>.gb."""
import os, re, glob, json, importlib.util, warnings
import statistics as st
from multiprocessing import Pool
warnings.filterwarnings("ignore")
REPO=os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
EVAL=os.environ.get("PLASTANNO_DATA","benchmark_data")+"/eval_v2_raw"
H=f"{REPO}/bench_runs/h2h"
spec=importlib.util.spec_from_file_location("bm", f"{REPO}/scripts/benchmark/benchmark_gene_by_gene.py")
bm=importlib.util.module_from_spec(spec); spec.loader.exec_module(bm)
from collections import Counter

def f1(tp,ref,pred):
    sn=tp/ref if ref else 0; pr=tp/pred if pred else 0
    return 200*sn*pr/(sn+pr) if sn+pr else 0

def score_pga(acc):
    ref_gb=f"{EVAL}/{acc}.gb"; pred_gb=f"{H}/pga_out/{acc}.gb"
    if not (os.path.exists(ref_gb) and os.path.exists(pred_gb)): return None
    try:
        warnings.filterwarnings("ignore")
        ref=bm.load_units(ref_gb); pred=bm.load_units(pred_gb)
        tp,fp,fn,_,_=bm.run(ref,pred,60,0.6)
        return dict(acc=acc, tp=len(tp), ref=len(ref), pred=len(pred),
                    tp_by_type=dict(Counter(rr["type"] for _p,rr,_s in tp)),
                    ref_by_type=dict(Counter(u["type"] for u in ref)),
                    pred_by_type=dict(Counter(u["type"] for u in pred)))
    except Exception:
        return None

if __name__=="__main__":
    subset=[l.strip() for l in open(f"{REPO}/splits/h2h_subset.txt") if l.strip()]
    # Plastanno predictions: dedicated self-contained head-to-head run of the
    # frozen tool on the same subset (bench_runs/h2h/plast_out), scored here with
    # the identical scorer used for PGA.
    def plast_pred(acc):
        p = f"{H}/plast_out/{acc}/{acc}.gb"
        return p if os.path.exists(p) else None
    def score_plast(acc):
        ref=f"{EVAL}/{acc}.gb"; pred=plast_pred(acc)
        if not (pred and os.path.exists(ref)): return None
        try:
            warnings.filterwarnings("ignore")
            r=bm.load_units(ref); pp=bm.load_units(pred); tp,fp,fn,_,_=bm.run(r,pp,60,0.6)
            return dict(acc=acc, tp=len(tp), ref=len(r), pred=len(pp),
                        tp_by_type=dict(Counter(x["type"] for _a,x,_b in tp)),
                        ref_by_type=dict(Counter(u["type"] for u in r)),
                        pred_by_type=dict(Counter(u["type"] for u in pp)))
        except Exception: return None
    with Pool(24) as p:
        plast_l=[x for x in p.map(score_plast, subset) if x]
        pga_list=[x for x in p.map(score_pga, subset) if x]
    plast={x["acc"]:x for x in plast_l}
    pga={x["acc"]:x for x in pga_list}
    common=[a for a in subset if a in plast and a in pga]
    print(f"subset={len(subset)} | Plastanno scored={len(plast)} | PGA scored={len(pga)} | paired={len(common)}")
    def micro(src, keys):
        out={}
        for t in ("CDS","tRNA","rRNA"):
            tp=sum(src[a]["tp_by_type"].get(t,0) for a in keys)
            ref=sum(src[a]["ref_by_type"].get(t,0) for a in keys)
            pred=sum(src[a]["pred_by_type"].get(t,0) for a in keys)
            out[t]=f1(tp,ref,pred)
        tp=sum(src[a]["tp"] for a in keys); ref=sum(src[a]["ref"] for a in keys); pred=sum(src[a]["pred"] for a in keys)
        out["GLOBAL"]=f1(tp,ref,pred); return out
    for tool,src in (("Plastanno",plast),("PGA",pga)):
        m=micro(src,common)
        print(f"\n=== {tool} (micro, n={len(common)}) ===  CDS {m['CDS']:.1f}  tRNA {m['tRNA']:.1f}  rRNA {m['rRNA']:.1f}  GLOBAL {m['GLOBAL']:.1f}")
    dp=[f1(plast[a]['tp'],plast[a]['ref'],plast[a]['pred']) for a in common]
    dg=[f1(pga[a]['tp'],pga[a]['ref'],pga[a]['pred']) for a in common]
    w=sum(1 for a,b in zip(dp,dg) if a>b+0.5); t=sum(1 for a,b in zip(dp,dg) if abs(a-b)<=0.5)
    print(f"\nper-genome macro F1: Plastanno {st.mean(dp):.1f} (med {st.median(dp):.1f}) | PGA {st.mean(dg):.1f} (med {st.median(dg):.1f})")
    print(f"Plastanno wins {w} | ties {t} | PGA wins {len(common)-w-t}  (of {len(common)})")
    with open(f"{H}/h2h_paired.csv","w") as o:
        o.write("accession,plastanno_f1,pga_f1\n")
        for a,pp,gg in zip(common,dp,dg): o.write(f"{a},{pp:.2f},{gg:.2f}\n")
    print(f"WROTE {H}/h2h_paired.csv")
