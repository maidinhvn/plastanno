"""
Engine A: Reference-based annotation using Exonerate.

Key improvements over v1:
- Pre-built protein DB (not on-the-fly from relatives)
- Search BOTH IRs for IR genes
- Returns Feature objects with s_ref score
"""
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import List, Dict

from ..core.feature import Feature


# Exonerate executable, resolved in order: $PLASTANNO_EXONERATE override, then the
# `exonerate` found on PATH, then a last-resort default. Keeps the tool portable
# while still working out-of-the-box where exonerate is already on PATH.
EXONERATE = (os.environ.get("PLASTANNO_EXONERATE")
             or shutil.which("exonerate")
             or "/data06/biotools/anaconda3/bin/exonerate")
BUFFER    = 2000   # bp buffer around search region

# IR genes need to be searched in both IRb and IRa
IR_GENES = {
    "rpl2","rpl23","ndhB","rps7","ycf2",
    "trnI-CAU","trnI-GAU","trnA-UGC",
    "trnR-ACG","trnN-GUU","trnL-CAA","trnV-GAC",
    "rrn16","rrn23","rrn4.5","rrn5",
}

RRNA_PRODUCTS = {
    "rrn16" : "16S ribosomal RNA",
    "rrn23" : "23S ribosomal RNA",
    "rrn4.5": "4.5S ribosomal RNA",
    "rrn5"  : "5S ribosomal RNA",
}

# Minimum gap (bp) between two physical copies of an rRNA gene when separating the
# IR copies by spatial clustering in detect_rrna. Larger than the longest rRNA gene
# (rrn23 ≈ 2.9 kb) so a single copy is never split, yet far below the inverted-repeat
# copy separation (≥10 kb) so the two copies are never merged.
_COPY_GAP = 4000


def get_search_regions(gene_name, ir_boundaries, genome_len):
    """
    Return list of (start, end) regions to search.
    IR genes: search both IRb and IRa.
    Others  : search their expected region only.
    """
    def region_coords(region):
        coords = ir_boundaries.get(region)
        if not coords: return None
        s, e = coords
        return max(0, s-BUFFER), min(genome_len, e+BUFFER)

    if gene_name in IR_GENES:
        regions = []
        for r in ["IRb", "IRa"]:
            c = region_coords(r)
            if c: regions.append(c)
        return regions if regions else [(0, genome_len)]

    # Determine region from gene catalog (passed via genome_len hack)
    # Default: search full genome (safe fallback)
    return [(0, genome_len)]


def run_exonerate_region(genome_seq, protein_seq,
                          gene_name, prot_acc,
                          region_start, region_end):
    """Run Exonerate on a specific genomic region."""
    buf_s = max(0, region_start - BUFFER)
    buf_e = min(len(genome_seq), region_end + BUFFER)
    region = genome_seq[buf_s:buf_e]

    with tempfile.TemporaryDirectory() as tmpdir:
        genome_fa  = os.path.join(tmpdir, "region.fasta")
        protein_fa = os.path.join(tmpdir, "protein.fasta")

        with open(genome_fa,  "w") as f:
            f.write(f">region\n{region}\n")
        with open(protein_fa, "w") as f:
            f.write(f">{prot_acc}\n{protein_seq}\n")

        result = subprocess.run([
            EXONERATE,
            "--model",          "protein2genome",
            "--query",          protein_fa,
            "--target",         genome_fa,
            "--score",          "30",
            "--percent",        "20",
            "--bestn",          "2",
            "--showalignment",  "no",
            "--showvulgar",     "no",
            "--showtargetgff",  "yes",
            "--softmasktarget", "no",
        ], capture_output=True, text=True, timeout=60)

        if result.returncode != 0:
            return []

        hits = _parse_exonerate_gff(
            result.stdout, gene_name, prot_acc, buf_s
        )
    return hits


def _parse_exonerate_gff(gff_text, gene_name, prot_acc, offset):
    """
    Parse an Exonerate protein2genome target-GFF block into Feature objects.

    Exonerate emits, per alignment: one 'gene' line (outer bounds + identity/
    similarity) followed by one 'cds' (and 'exon') line per coding segment, with
    'intron' lines between them. The previous version read only the 'gene' line,
    so multi-exon genes were flattened into a single span that included introns.
    This version collects every coding segment into Feature.exons (0-based,
    half-open, genome coordinates, ascending) so spliced genes are correct.
    """
    hits = []
    cur = None

    def finalize(h):
        if h is None:
            return
        exons = sorted(h["cds"] or h["exon"] or [(h["gstart"], h["gend"])])
        feat = Feature(
            gene_name=gene_name, gene_type="CDS",
            start=exons[0][0], end=exons[-1][1], strand=h["strand"],
            engine="A", s_ref=h["sim"],
        )
        feat.exons = exons
        feat.has_intron = len(exons) > 1
        hits.append(feat)

    for line in gff_text.split("\n"):
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 9:
            continue
        ftype = parts[2].lower()
        if ftype not in ("gene", "cds", "exon"):
            continue
        start = int(parts[3]) - 1 + offset
        end   = int(parts[4]) + offset
        strand = 1 if parts[6] == "+" else -1

        if ftype == "gene":
            finalize(cur)
            sim = 0.0
            for attr in parts[8].split(";"):
                if "similarity" in attr.lower():
                    try:
                        sim = float(attr.split()[-1]) / 100
                    except Exception:
                        pass
            cur = {"gstart": start, "gend": end, "strand": strand,
                   "sim": sim, "cds": [], "exon": []}
        elif cur is not None:
            cur[ftype].append((start, end))

    finalize(cur)
    return hits

def run_exonerate_gene(genome_seq, gene_name,
                        protein_db, ir_boundaries,
                        gene_catalog, threads=4, ref_proteins=None):
    """
    Run Exonerate for one gene using pre-built protein DB.
    If ref_proteins supplies a protein for this gene (from a user-provided
    --reference), it is tried first and the built-in DB proteins act as a
    fallback — so a closer reference is used when available without ever losing
    the default coverage.
    Returns list of Feature objects.
    """
    from Bio import SeqIO
    # Built-in per-gene protein DB (fallback / default)
    prot_fasta = Path(protein_db) / f"{gene_name}.fasta"
    proteins = []
    if prot_fasta.exists():
        try:
            proteins = list(SeqIO.parse(str(prot_fasta), "fasta"))
        except Exception:
            proteins = []

    # Overlay user-reference protein(s) first (auto fallback to DB below)
    if ref_proteins and gene_name in ref_proteins:
        proteins = list(ref_proteins[gene_name]) + proteins

    if not proteins:
        return []

    # Get search regions
    genome_len = len(genome_seq)
    cat        = gene_catalog.get(gene_name, {})
    region_key = cat.get("region", "LSC")

    if gene_name in IR_GENES or region_key in ("IRb","IR"):
        search_regs = []
        for r in ["IRb","IRa"]:
            c = ir_boundaries.get(r)
            if c: search_regs.append(c)
        if not search_regs:                       # IR-lacking: search whole genome
            search_regs = [(0, genome_len)]
    else:
        c = ir_boundaries.get(region_key)
        search_regs = [c] if c else [(0, genome_len)]

    # Run Exonerate on each region
    all_hits = []
    for reg_s, reg_e in search_regs:
        for prot in proteins[:5]:  # max 5 proteins per region
            hits = run_exonerate_region(
                genome_seq  = genome_seq,
                protein_seq = str(prot.seq),
                gene_name   = gene_name,
                prot_acc    = prot.id,
                region_start= reg_s,
                region_end  = reg_e,
            )
            all_hits.extend(hits)

    # Deduplicate by position
    deduped = {}
    for hit in all_hits:
        key = hit.start // 500
        if key not in deduped or \
           hit.s_ref > deduped[key].s_ref:
            deduped[key] = hit

    return list(deduped.values())


def detect_rrna(genome_seq, rrna_dbs, threads=4):
    """Detect rRNA genes using BLAST full-length databases.

    BLAST returns each rRNA copy as one or several HSPs (a long gene is split
    wherever local mismatches break the alignment). The previous version kept a
    single best-pident HSP per 25 kb bin, which had two failure modes: (a) when
    the best HSP covered only part of the gene the feature was truncated (e.g.
    rrn16 reported as 67 bp), and (b) a short spurious opposite-strand HSP could
    out-score the real alignment and flip the strand. Both share one root cause —
    picking one HSP instead of the full-length alignment.

    This version, per 25 kb bin (≈ one IR copy), chooses the dominant strand by
    total bitscore (the real full-length alignment outweighs a short spurious
    hit), then clusters the co-linear HSPs on that strand (rRNA has no introns,
    so the gene is one contiguous block) and reports the highest-scoring cluster
    spanning min start → max end. Boundaries and strand therefore both reflect
    the real alignment.
    """
    import tempfile, os
    from collections import defaultdict

    features = []
    with tempfile.NamedTemporaryFile(
        suffix=".fasta", mode="w", delete=False
    ) as f:
        f.write(f">query\n{genome_seq}\n")
        qfa = f.name

    try:
        for gene, db in rrna_dbs.items():
            if not Path(db + ".nhr").exists(): continue

            r = subprocess.run([
                "blastn", "-query", qfa, "-db", db,
                "-outfmt",
                "6 sseqid qstart qend sstrand pident length bitscore",
                "-max_target_seqs", "50",
                "-perc_identity",   "80",
                "-num_threads",     str(threads),
            ], capture_output=True, text=True)

            # Separate the (usually two IR) copies of this gene by SPATIAL
            # clustering rather than a fixed 25 kb bin. A fixed bin merged the two
            # copies whenever the SSC was small — e.g. the operon-terminal rrn5 /
            # rrn4.5 copies can sit only ~15 kb apart, so both fell in one 25 kb
            # bin and the pipeline emitted a single copy (the other became a false
            # negative). Single-linkage with a gap of _COPY_GAP separates physical
            # copies (IR copies are ≥10 kb apart) without splitting one gene (all
            # of a gene's HSPs lie within its ≤3 kb span, far below _COPY_GAP).
            raw = []
            for line in r.stdout.strip().split("\n"):
                p = line.split("\t")
                if len(p) < 7: continue
                length = int(p[5])
                if length < 50: continue
                raw.append({
                    "sid":p[0], "qs":int(p[1]) - 1, "qe":int(p[2]),
                    "strand":1 if p[3]=="plus" else -1,
                    "pident":float(p[4]), "length":length, "bitscore":float(p[6]),
                })

            raw.sort(key=lambda h: h["qs"])
            copies = []
            for h in raw:
                if copies and h["qs"] - copies[-1]["_end"] <= _COPY_GAP:
                    copies[-1]["hsps"].append(h)
                    copies[-1]["_end"] = max(copies[-1]["_end"], h["qe"])
                else:
                    copies.append({"hsps":[h], "_end":h["qe"]})

            for copy in copies:
                hsps = copy["hsps"]
                # Reconstruct each DB subject's view of this copy: a single good
                # reference covers the whole gene, but its alignment may be split
                # into several HSPs by local mismatches. Merge each subject's
                # co-linear HSPs (same strand, ≤200 bp gap) into one span; do NOT
                # merge across subjects (their differing lengths would inflate the
                # boundaries). Then keep the highest-scoring subject — full-length
                # (merged within subject) yet tight (one reference's boundaries).
                by_sid = defaultdict(list)
                for h in hsps:
                    by_sid[h["sid"]].append(h)

                candidates = []
                for sub_hsps in by_sid.values():
                    # dominant strand for this subject, then cluster co-linear HSPs
                    by_strand = defaultdict(float)
                    for h in sub_hsps:
                        by_strand[h["strand"]] += h["bitscore"]
                    sstrand = max(by_strand, key=by_strand.get)
                    same = sorted((h for h in sub_hsps if h["strand"] == sstrand),
                                  key=lambda h: h["qs"])
                    clusters = []
                    for h in same:
                        if clusters and h["qs"] <= clusters[-1]["qe"] + 200:
                            c = clusters[-1]
                            c["qe"]          = max(c["qe"], h["qe"])
                            c["bitscore"]   += h["bitscore"]
                            c["pident_sum"] += h["pident"] * h["length"]
                            c["len_sum"]    += h["length"]
                        else:
                            clusters.append({
                                "qs":h["qs"], "qe":h["qe"], "bitscore":h["bitscore"],
                                "pident_sum":h["pident"]*h["length"], "len_sum":h["length"],
                                "strand":sstrand,
                            })
                    candidates.append(max(clusters, key=lambda c: c["bitscore"]))

                best = max(candidates, key=lambda c: c["bitscore"])
                pident = best["pident_sum"] / best["len_sum"]
                features.append(Feature(
                    gene_name = gene,
                    gene_type = "rRNA",
                    product   = RRNA_PRODUCTS.get(gene, gene),
                    start     = best["qs"],
                    end       = best["qe"],
                    strand    = best["strand"],
                    engine    = "A",
                    s_ref     = pident / 100,
                ))
    finally:
        os.unlink(qfa)

    return features


def load_reference_proteins(ref_path):
    """Extract per-gene protein sequences from an annotated reference GenBank, for
    overlaying onto the built-in protein DB (the --reference option). Returns
    {gene_name: [SeqRecord, ...]}; empty dict if the file cannot be read as
    GenBank (e.g. a plain FASTA carries no gene annotations to extract)."""
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    out = {}
    try:
        recs = list(SeqIO.parse(ref_path, "genbank"))
    except Exception:
        recs = []
    for rec in recs:
        for f in rec.features:
            if f.type != "CDS":
                continue
            g = (f.qualifiers.get("gene", [None])[0]
                 or f.qualifiers.get("product", [None])[0])
            if not g:
                continue
            tr = (f.qualifiers.get("translation", [None])[0])
            if not tr:
                try:
                    s = str(f.extract(rec.seq))
                    tr = str(Seq(s[: len(s)//3*3]).translate(table=11)).rstrip("*")
                except Exception:
                    continue
            if tr:
                out.setdefault(g, []).append(SeqRecord(Seq(tr), id=f"ref_{g}", description=""))
    return out


def run_engine_a(
    genome_seq,
    ir_boundaries,
    relatives,
    protein_db,
    gene_catalog,
    threads = 4,
    ref_proteins = None,
) -> List[Feature]:
    """
    Run Engine A: reference-based annotation.
    1. CDS via Exonerate (pre-built protein DB, optionally overlaid by ref_proteins)
    2. rRNA via BLAST (full-length DB)
    """
    features = []

    # 1. CDS annotation
    from pathlib import Path
    prot_dir = Path(protein_db)
    for gene_fa in sorted(prot_dir.glob("*.fasta")):
        gene_name = gene_fa.stem
        hits = run_exonerate_gene(
            genome_seq   = genome_seq,
            gene_name    = gene_name,
            protein_db   = protein_db,
            ir_boundaries= ir_boundaries,
            gene_catalog = gene_catalog,
            threads      = threads,
            ref_proteins = ref_proteins,
        )
        features.extend(hits)

    return features
