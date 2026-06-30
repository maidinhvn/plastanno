"""
Engine B: Model-based annotation.

Components:
1. HMM scan for CDS (HMMER nhmmer)
2. tRNA detection (ARAGORN + BLAST hierarchical + exon DB)
3. rRNA detection (BLAST full-length)

Key improvements over v1:
- HMM used for DETECTION (not just naming)
- Hierarchical tRNA DB with taxonomy awareness
- Actual exon sequences for intron tRNA
- Returns Feature objects with s_model score
"""
import os
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import List, Dict

from Bio.Seq import Seq

from ..core.feature import Feature


# ── HMM-based CDS detection ───────────────────────────────────────────────────

def run_hmm_scan(genome_seq, hmm_dir, threads=4):
    """
    HMM-based CDS detection (clean rewrite, Option B).

    The profile DB is amino-acid (ALPH amino), so the previous nhmmer
    (nucleotide-only) call always failed with 'Invalid alphabet'. Here the genome
    is translated in all six reading frames (plastid table 11), split into
    stop-free ORF segments, and searched with hmmsearch (protein). Domain
    coordinates are mapped back to genome nucleotide coordinates.

    Returns list[Feature] (gene_type='CDS', engine='B', s_model set).
    """
    try:
        from Bio.Seq import Seq
    except Exception:
        return []

    hmm_db = Path(hmm_dir).parent / "all_profiles.hmm"
    if not hmm_db.exists():
        profiles = sorted(Path(hmm_dir).glob("*.hmm"))
        if not profiles:
            return []
        with open(hmm_db, "w") as out:
            for prof in profiles:
                out.write(open(prof).read())

    L   = len(genome_seq)
    fwd = genome_seq.upper()
    rev = str(Seq(fwd).reverse_complement())
    MIN_AA = 20

    targets = []  # (name, protein_segment)
    for strand_char, seqstr in (("f", fwd), ("r", rev)):
        for k in range(3):
            sub = seqstr[k:]
            sub = sub[: len(sub) // 3 * 3]
            if not sub:
                continue
            prot = str(Seq(sub).translate(table=11))
            off = 0
            for piece in prot.split("*"):
                if len(piece) >= MIN_AA:
                    targets.append(("%s%d_%d" % (strand_char, k, off), piece))
                off += len(piece) + 1
    if not targets:
        return []

    qfa = tempfile.NamedTemporaryFile("w", suffix=".faa", delete=False)
    dom = tempfile.NamedTemporaryFile("w", suffix=".domtbl", delete=False)
    dom.close()
    raw = []
    try:
        for name, seq in targets:
            qfa.write(">%s\n%s\n" % (name, seq))
        qfa.close()

        subprocess.run(
            ["hmmsearch", "--cpu", str(threads), "-E", "1e-5",
             "--domtblout", dom.name, "--noali", str(hmm_db), qfa.name],
            capture_output=True, text=True, timeout=600,
        )

        for line in open(dom.name):
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            if len(p) < 23:
                continue
            gene     = p[3]
            qlen     = int(p[5])
            i_eval   = float(p[12])
            dscore   = float(p[13])
            env_from = int(p[19])
            env_to   = int(p[20])
            if i_eval > 1e-5:
                continue
            if (env_to - env_from + 1) < 0.5 * qlen:
                continue
            frame_id, seg_off = p[0].rsplit("_", 1)
            gs, ge, strand = _map_to_genome(frame_id, int(seg_off), env_from, env_to, L)
            if gs < 0 or ge > L or ge <= gs:
                continue
            raw.append((gene, gs, ge, strand, min(1.0, dscore / 200.0)))
    except subprocess.TimeoutExpired:
        return []
    finally:
        for fp in (qfa.name, dom.name):
            try:
                os.unlink(fp)
            except OSError:
                pass

    # de-duplicate per gene: keep best-scoring, non-overlapping copies (handles IR)
    features, by_gene = [], {}
    for rec in raw:
        by_gene.setdefault(rec[0], []).append(rec)
    for gene, recs in by_gene.items():
        recs.sort(key=lambda r: r[4], reverse=True)
        kept = []
        for _, gs, ge, strand, sm in recs:
            if any(not (ge <= kgs or gs >= kge) for kgs, kge in kept):
                continue
            kept.append((gs, ge))
            features.append(Feature(gene_name=gene, gene_type="CDS",
                                    start=gs, end=ge, strand=strand,
                                    engine="B", s_model=sm))
    return features


def _map_to_genome(frame_id, seg_off, env_from, env_to, L):
    """Map AA domain env coords (1-based within an ORF segment) to genome nt coords."""
    strand_char, k = frame_id[0], int(frame_id[1:])
    a0 = seg_off + (env_from - 1)
    a1 = seg_off + (env_to - 1)
    if strand_char == "f":
        return k + a0 * 3, k + (a1 + 1) * 3, 1
    rc_s = k + a0 * 3
    rc_e = k + (a1 + 1) * 3
    return L - rc_e, L - rc_s, -1

def _aragorn_intron_exons(start, end, strand, off, ilen):
    """Given ARAGORN's 1-based-inclusive [start,end], strand, intron offset (from
    the 5' end of the mature tRNA) and intron length, return 0-based half-open
    exon list (ascending genome coord) or None if the geometry is degenerate.

    ARAGORN reports the intron as i(off,ilen); the mature length is the span
    minus the intron, split into a 5' exon of length `off` and a 3' exon of the
    rest. The 5' exon sits at the low-coord end on +strand, high-coord end on -.
    """
    mature = (end - start + 1) - ilen
    e2len = mature - off
    if off <= 0 or e2len <= 0 or ilen <= 0:
        return None
    if strand == 1:
        ex5 = (start, start + off - 1)          # 5' exon, low coord
        ex3 = (end - e2len + 1, end)            # 3' exon, high coord
    else:
        ex5 = (end - off + 1, end)             # 5' exon, high coord
        ex3 = (start, start + e2len - 1)       # 3' exon, low coord
    lo, hi = sorted([ex5, ex3])
    return [(lo[0] - 1, lo[1]), (hi[0] - 1, hi[1])]  # 0-based half-open, ascending


def run_aragorn(genome_seq, genome_len):
    """
    Run ARAGORN for tRNA detection, including intron-containing tRNAs (-i).

    Clean-rewrite parser for ARAGORN v1.2.x batch (-w) output. Lines look like:
        1   tRNA-His                     c[31,105]      35      (gtg)
        4   tRNA-Arg                  [9981,10052]      34      (tct)
        2   tRNA-Lys                  c[1696,4271]      33      (ttt)i(38,2504)
    Columns: index, tRNA-<AA3>, [c]?[start,end] (1-based), length, (anticodon),
    and — with -i — an optional i(intron_offset,intron_length). A leading 'c' on
    the coordinate means the minus strand. The anticodon is converted to RNA for
    the plastome gene name (tRNA-His + gtg -> trnH-GUG).

    Returns list[Feature] (gene_type='tRNA', engine='B'), 0-based [start, end).
    For intron tRNAs the structural call (coords + exons) is from ARAGORN; the
    anticodon ARAGORN reports can be wrong for a few intron tRNAs (trnI-GAU,
    trnG-UCC), so the name is corrected by BLAST downstream
    (`_rename_intron_trnas_by_blast`).
    """
    AA3TO1 = {
        "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
        "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
        "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
        "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
        "Sec": "U",
    }

    features = []
    tmp = tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False)
    try:
        tmp.write(">query\n%s\n" % genome_seq)
        tmp.close()

        result = subprocess.run(
            ["aragorn", "-i", "-t", "-l", "-w", tmp.name],
            capture_output=True, text=True, timeout=120,
        )

        line_re = re.compile(
            r"^\s*\d+\s+"                  # index
            r"(?:tRNA|tmRNA)-(\w+?)\d*\s+" # tRNA-<AA3> (drop trailing digit)
            r"(c?)\[(\d+),(\d+)\]\s+"      # [c][start,end]
            r"\d+\s+"                       # length
            r"\(([a-zA-Z]+)\)"              # (anticodon)
            r"(?:i\((\d+),(\d+)\))?"        # optional intron i(offset,length)
        )

        for line in result.stdout.splitlines():
            m = line_re.match(line)
            if not m:
                continue
            aa3, comp, s_str, e_str, anti, i_off, i_len = m.groups()
            start, end = int(s_str), int(e_str)
            if end < start:
                start, end = end, start
            strand = -1 if comp == "c" else 1
            anti_rna = anti.upper().replace("T", "U")

            if aa3 in ("fMet", "fMeth"):
                code = "fM"
            else:
                code = AA3TO1.get(aa3.capitalize())
            gene_name = ("trn%s-%s" % (code, anti_rna)) if code else ("trn?-%s" % anti_rna)

            exons = None
            if i_off is not None:
                exons = _aragorn_intron_exons(
                    start, end, strand, int(i_off), int(i_len)
                )

            f = Feature(
                gene_name=gene_name,
                gene_type="tRNA",
                product="tRNA-%s" % aa3,
                start=(exons[0][0] if exons else start),
                end=(exons[-1][1] if exons else end),
                strand=strand,
            )
            f.engine = "B"
            if exons:
                f.exons = exons
                f.has_intron = True
                f.notes.append("ARAGORN intron i(%s,%s)" % (i_off, i_len))
            else:
                f.notes.append("ARAGORN")
            features.append(f)
    finally:
        try:
            os.unlink(tmp.name)
        except OSError:
            pass

    return features

def _normalize_trna_name(raw):
    """Convert ARAGORN tRNA names to standard format."""
    ANTICODON_MAP = {
        "Ala":"A","Arg":"R","Asn":"N","Asp":"D",
        "Cys":"C","Gln":"Q","Glu":"E","Gly":"G",
        "His":"H","Ile":"I","Leu":"L","Lys":"K",
        "Met":"M","Phe":"F","Pro":"P","Ser":"S",
        "Thr":"T","Trp":"W","Tyr":"Y","Val":"V",
        "fMet":"fM","Sec":"Sec",
    }
    # Extract anticodon from parentheses
    m = re.search(r'\(([A-Za-z]{3})\)', raw)
    if not m: return None
    anticodon = m.group(1).upper().replace("T","U")

    # Extract amino acid
    m2 = re.search(r'tRNA-([A-Za-z]+)', raw)
    if not m2: return None
    aa = m2.group(1)

    # Map 3-letter to 1-letter
    aa1 = ANTICODON_MAP.get(aa, aa[0].upper())
    return f"trn{aa1}-{anticodon}"


def blast_trna(genome_seq, trna_db, min_length=60,
               threads=4):
    """
    BLAST genome against tRNA DB.
    Returns list of Feature objects.
    """
    features = []

    with tempfile.NamedTemporaryFile(
        suffix=".fasta", mode="w", delete=False
    ) as f:
        f.write(f">genome\n{genome_seq}\n")
        qfa = f.name

    try:
        result = subprocess.run([
            "blastn", "-query", qfa, "-db", trna_db,
            "-outfmt",
            "6 qstart qend sseqid pident length sstrand",
            "-perc_identity", "75",
            "-word_size",     "11",
            "-dust",          "no",
            "-num_threads",   str(threads),
        ], capture_output=True, text=True, timeout=60)

        # Deduplicate by position bin
        best = {}
        for line in result.stdout.strip().split("\n"):
            if not line: continue
            p = line.split("\t")
            if len(p) < 6: continue
            qs     = int(p[0]) - 1
            qe     = int(p[1])
            gene   = _extract_gene_from_sid(p[2])
            pident = float(p[3])
            length = int(p[4])
            strand = 1 if p[5]=="plus" else -1

            if length < min_length or not gene: continue

            key = (gene, qs // 100)
            if key not in best or pident > best[key]["pident"]:
                best[key] = {
                    "gene":gene, "qs":qs, "qe":qe,
                    "strand":strand, "pident":pident
                }

        for h in best.values():
            features.append(Feature(
                gene_name = h["gene"],
                gene_type = "tRNA",
                start     = h["qs"],
                end       = h["qe"],
                strand    = h["strand"],
                engine    = "B",
                s_model   = h["pident"] / 100,
            ))

    except subprocess.TimeoutExpired:
        pass
    finally:
        os.unlink(qfa)

    return features


# Largest plausible tRNA exon (bp). Some exon-DB entries are malformed — a few
# "exon1" subjects bundle part of the intron and align as 400+ bp blocks; the real
# exon is the short portion adjacent to the intron, so an over-long hit is trimmed
# to its intron-proximal end rather than trusted whole.
_MAX_TRNA_EXON = 55

# Per-gene cap on the intron length (gap between the two exon hits). All plastid
# intron tRNAs have introns < ~1 kb (largest: trnI-GAU ≈ 950 bp) EXCEPT trnK-UUU,
# whose group-I intron carries matK and runs ~2.5 kb. A single permissive global
# cap (3.5 kb, set for trnK) also let a 5' exon mis-pair with a wrong 3' exon
# 3+ kb away, producing a giant spurious tRNA (e.g. a 3.5 kb "trnL-UAA") that then
# swallowed an adjacent real tRNA in _merge_trna_sources (trnF-GAA was lost this
# way). Capping each gene to its biologically plausible intron kills the mis-pair
# without touching trnK.
_INTRON_MAX_BY_GENE = {"trnK-UUU": 3000}
_INTRON_MAX_DEFAULT = 1200


def blast_intron_trna(genome_seq, exon_db,
                       min_intron=200, max_intron=3500,
                       threads=4):
    """
    Detect intron-containing tRNAs using actual exon sequences.

    Pairs a 5' and 3' exon hit of the same gene separated by an intron. Three
    robustness fixes over the naive version: (a) the minimum hit length is 15 bp
    so genuinely short first exons (e.g. trnG-UCC exon1 ≈ 24 bp) are not dropped;
    (b) an over-long hit (malformed DB entry that swallows intron sequence) is
    trimmed to its intron-proximal end, recovering the true short exon; (c) for
    each gene/strand the *tightest* valid intron is chosen instead of the first
    pair encountered, which stops a spurious upstream hit from mis-pairing with
    the real 3' exon (the trnI-GAU / trnA-UGC operon confusion).
    """
    if not Path(exon_db + ".nhr").exists():
        return []

    with tempfile.NamedTemporaryFile(
        suffix=".fasta", mode="w", delete=False
    ) as f:
        f.write(f">genome\n{genome_seq}\n")
        qfa = f.name

    raw_hits = []
    try:
        result = subprocess.run([
            "blastn", "-query", qfa, "-db", exon_db,
            "-outfmt",
            "6 qstart qend sseqid pident length sstrand",
            "-perc_identity", "80",
            "-word_size",     "11",
            "-dust",          "no",
            "-num_threads",   str(threads),
        ], capture_output=True, text=True, timeout=120)

        for line in result.stdout.strip().split("\n"):
            if not line: continue
            p = line.split("\t")
            if len(p) < 6: continue
            qs     = int(p[0]) - 1
            qe     = int(p[1])
            gene   = _extract_gene_from_sid(p[2])
            exon   = _exon_num_from_sid(p[2])     # 1 (5') / 2 (3') / None
            pident = float(p[3])
            length = int(p[4])
            strand = 1 if p[5]=="plus" else -1

            if length < 15 or not gene: continue
            raw_hits.append({
                "gene":gene, "exon":exon, "qs":qs, "qe":qe,
                "strand":strand, "pident":pident
            })

    except subprocess.TimeoutExpired:
        pass
    finally:
        os.unlink(qfa)

    # Group by gene + strand, then pick the single best (tightest-intron) exon
    # pair per locus. The intron sits between the left hit's right end and the
    # right hit's left end, so each exon is trimmed at that shared (proximal) side.
    groups = defaultdict(list)
    for h in raw_hits:
        groups[(h["gene"], h["strand"])].append(h)

    features = []
    for (gene, strand), hits in groups.items():
        # Biologically plausible intron length for THIS gene (trnK is the only
        # large one); never exceed the caller-supplied absolute ceiling.
        cap = min(max_intron, _INTRON_MAX_BY_GENE.get(gene, _INTRON_MAX_DEFAULT))
        hits = sorted(hits, key=lambda x: x["qs"])
        best = None  # (gap, -combined_pident, left_exon, right_exon)
        for i, a in enumerate(hits):
            for b in hits[i + 1:]:
                gap = b["qs"] - a["qe"]
                if gap < min_intron:
                    continue
                if gap > cap:
                    break                      # hits sorted by qs → no closer b later
                # A real intron tRNA spans exactly one 5' exon and one 3' exon.
                # The exon-DB subjects are tagged _e1/_e2, so reject same-exon
                # pairs (e1+e1, e2+e2): those are spurious cross-hits that the
                # "tightest intron" rule would otherwise prefer over the true
                # wider pair (the trnA-UGC / trnI-GAU operon mis-pairing).
                if (a["exon"] and b["exon"] and a["exon"] == b["exon"]):
                    continue
                # trim each exon to its intron-proximal end (left hit's right side,
                # right hit's left side)
                la = (max(a["qs"], a["qe"] - _MAX_TRNA_EXON), a["qe"])
                rb = (b["qs"], min(b["qe"], b["qs"] + _MAX_TRNA_EXON))
                key = (gap, -(a["pident"] + b["pident"]), la, rb)
                if best is None or key < best:
                    best = key
        if best is None:
            continue
        _, _, la, rb = best
        features.append(Feature(
            gene_name=gene, gene_type="tRNA",
            start=la[0], end=rb[1], strand=strand,
            exons=[la, rb], has_intron=True, engine="B",
            s_model=1.0,
        ))

    # Deduplicate per locus (two IR copies stay separate; they don't overlap)
    seen = {}
    deduped = []
    for f in sorted(features, key=lambda x: x.s_model, reverse=True):
        key = (f.gene_name, f.start // 200)
        if key not in seen:
            seen[key] = True
            deduped.append(f)

    return deduped


def _exon_num_from_sid(sid):
    """Extract exon number (1 or 2) from an exon-DB subject ID like
    'NC_003386_trnG-UCC_e1'. Returns None if untagged."""
    m = re.search(r'_e([12])\b', sid)
    return int(m.group(1)) if m else None


def _extract_gene_from_sid(sid):
    """Extract gene name from BLAST subject ID."""
    # Try trnX-NNN pattern
    m = re.search(r'(trn[A-Z]-[A-Z]{3})', sid)
    if m: return m.group(1)
    # Try last underscore-delimited field
    parts = sid.split("_")
    if len(parts) >= 2:
        candidate = parts[-1]
        if candidate.startswith("trn"):
            return candidate
    return None


def select_trna_db(relatives, trna_db_dir):
    """
    Select best tRNA DB based on taxonomy.
    Priority: genus → family → global
    """
    import json
    idx_path = Path(trna_db_dir) / "index.json"
    if not idx_path.exists():
        return str(Path(trna_db_dir) / "global")

    idx = json.load(open(idx_path))

    if relatives:
        genus  = relatives[0].get("genus", "")
        family = relatives[0].get("family","")

        genus_db = Path(trna_db_dir) / "genus" / genus
        if genus and genus in idx.get("genus_dbs",[]) and \
           genus_db.with_suffix(".nhr").exists():
            return str(genus_db)

        family_db = Path(trna_db_dir) / "family" / family
        if family and family in idx.get("family_dbs",[]) and \
           family_db.with_suffix(".nhr").exists():
            return str(family_db)

    return str(Path(trna_db_dir) / "global")


def detect_rrna_b(genome_seq, rrna_dbs, threads=4):
    """Detect rRNA using BLAST (same as Engine A)."""
    from .engine_a import detect_rrna
    feats = detect_rrna(genome_seq, rrna_dbs, threads)
    for f in feats:
        f.engine  = "B"
        f.s_model = f.s_ref
    return feats


# ── Main Engine B ─────────────────────────────────────────────────────────────

def run_trnascan(genome_seq):
    """Optional extra tRNA source: tRNAscan-SE in organellar mode (-O), widely
    regarded as the strongest tRNA predictor. Contributes intronless tRNAs;
    intron-containing hits are left to ARAGORN/exon-DB (tRNAscan does not model
    plastid group-II introns). Returns list[Feature]; silently returns [] if
    tRNAscan-SE is not on PATH or fails."""
    if not shutil.which("tRNAscan-SE"):
        print("        tRNAscan-SE not found on PATH — skipping")
        return []
    AA3TO1 = {
        "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
        "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
        "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
        "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
        "Sec": "U", "SeC": "U", "fMet": "fM",
    }
    features = []
    fa  = tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False)
    out = fa.name + ".trnascan.out"
    try:
        fa.write(">query\n%s\n" % genome_seq); fa.close()
        subprocess.run(["tRNAscan-SE", "-O", "-q", "-o", out, fa.name],
                       capture_output=True, text=True, timeout=600)
        if not os.path.exists(out):
            return []
        for line in open(out):
            parts = line.split()
            if len(parts) < 9 or not parts[1].isdigit():
                continue          # header / malformed
            try:
                begin, end = int(parts[2]), int(parts[3])
                aa, codon  = parts[4], parts[5]
                i_b, i_e   = int(parts[6]), int(parts[7])
            except ValueError:
                continue
            if i_b or i_e:        # intron present -> leave to ARAGORN/exon-DB
                continue
            if aa in ("Undet", "Sup"):
                continue
            strand = -1 if begin > end else 1
            s0, e0 = min(begin, end) - 1, max(begin, end)
            code = AA3TO1.get(aa, AA3TO1.get(aa.capitalize()))
            anti = codon.upper().replace("T", "U")
            gene = ("trn%s-%s" % (code, anti)) if code else ("trn?-%s" % anti)
            f = Feature(gene_name=gene, gene_type="tRNA",
                        product="tRNA-%s" % aa, start=s0, end=e0, strand=strand)
            f.engine = "B"
            f.notes.append("tRNAscan-SE")
            features.append(f)
    except Exception:
        return features
    finally:
        for p in (fa.name, out):
            try:
                os.unlink(p)
            except OSError:
                pass
    return features


def run_engine_b(
    genome_seq,
    ir_boundaries,
    relatives,
    hmm_dir,
    trna_db_dir,
    exon_db,
    rrna_dbs,
    gene_catalog,
    threads = 4,
    use_trnascan = False,
) -> List[Feature]:
    """
    Run Engine B: model-based annotation.
    1. HMM scan for CDS
    2. ARAGORN + BLAST for tRNA (with intron support)
    3. BLAST for rRNA
    """
    features = []
    genome_len = len(genome_seq)

    # 1. HMM scan for CDS
    hmm_feats = run_hmm_scan(genome_seq, hmm_dir, threads)
    features.extend(hmm_feats)
    print(f"        HMM: {len(hmm_feats)} CDS candidates")

    # 2. tRNA detection
    # 2a. ARAGORN
    aragorn_feats = run_aragorn(genome_seq, genome_len)
    print(f"        ARAGORN: {len(aragorn_feats)} tRNAs")

    # 2b. BLAST hierarchical tRNA DB
    trna_db = select_trna_db(relatives, trna_db_dir)
    db_name = Path(trna_db).name
    print(f"        tRNA DB: {db_name}")

    if Path(trna_db + ".nhr").exists():
        blast_feats = blast_trna(
            genome_seq, trna_db, threads=threads
        )
        print(f"        BLAST tRNA: {len(blast_feats)} candidates")
    else:
        blast_feats = []

    # 2c. Intron tRNA detection
    intron_feats = blast_intron_trna(
        genome_seq, exon_db, threads=threads
    )
    print(f"        Intron tRNA: {len(intron_feats)} detected")

    # 2d. tRNAscan-SE (optional extra source)
    trnascan_feats = run_trnascan(genome_seq) if use_trnascan else []
    if use_trnascan:
        print(f"        tRNAscan-SE: {len(trnascan_feats)} tRNAs")

    # Merge tRNA: intron > tRNAscan-SE > ARAGORN > BLAST
    trna_feats = _merge_trna_sources(
        aragorn_feats, blast_feats, intron_feats, trnascan_feats
    )
    # Correct intron-tRNA names (ARAGORN mis-reads a few anticodons) by BLAST.
    if Path(trna_db + ".nhr").exists():
        _rename_intron_trnas_by_blast(trna_feats, genome_seq, trna_db)
    print(f"        tRNA total: {len(trna_feats)}")
    features.extend(trna_feats)

    # 3. rRNA detection
    rrna_feats = detect_rrna_b(genome_seq, rrna_dbs, threads)
    print(f"        rRNA: {len(rrna_feats)}")
    features.extend(rrna_feats)

    return features


def _merge_trna_sources(aragorn, blast, intron, trnascan=None):
    """
    Merge tRNA candidates from the three sources, removing duplicates by LOCUS
    rather than by name. Two candidates overlapping by more than half of the
    shorter one are treated as the same physical tRNA (e.g. ARAGORN and BLAST
    disagreeing on the anticodon, as with trnM-CAU vs trnT-UGU at one locus), and
    only the most trustworthy one is kept. Spatially separate copies (the two IR
    copies of a tRNA) do not overlap and are both retained.

    Source trust: intron-aware (exon DB) > ARAGORN (structural) > BLAST (similarity).
    """
    SRC_INTRON, SRC_TRNASCAN, SRC_ARAGORN, SRC_BLAST = 3, 2.5, 2, 1
    cand = ([(f, SRC_INTRON)   for f in intron]
            + [(f, SRC_TRNASCAN) for f in (trnascan or [])]
            + [(f, SRC_ARAGORN)  for f in aragorn]
            + [(f, SRC_BLAST)    for f in blast])

    def same_locus(a, b):
        inter   = max(0, min(a.end, b.end) - max(a.start, b.start))
        shorter = min(a.end - a.start, b.end - b.start)
        return shorter > 0 and inter / shorter > 0.5

    def rank(item):
        f, pr = item
        return (pr, round(getattr(f, "s_model", 0.0) or 0.0, 3), f.end - f.start)

    cand.sort(key=rank, reverse=True)   # best source/score/length first
    kept = []
    for f, pr in cand:
        if any(same_locus(f, k) for k in kept):
            continue
        kept.append(f)
    return kept


def _spliced_seq(genome_seq, exons, strand):
    """Concatenate exon sequences (ascending coord) and reverse-complement for
    the minus strand, giving the mature 5'→3' tRNA sequence."""
    seq = "".join(genome_seq[s:e] for s, e in sorted(exons))
    if strand == -1:
        seq = str(Seq(seq).reverse_complement())
    return seq


def _rename_intron_trnas_by_blast(feats, genome_seq, trna_db):
    """Correct the gene name of intron-containing tRNAs by BLASTing their spliced
    (mature) sequence against the tRNA DB and taking the best-hit name. ARAGORN's
    structural intron calls give the right coordinates but can mis-read the
    anticodon (e.g. trnI-GAU → tRNA-Glu, trnG-UCC → tRNA-Ser); the DB carries the
    true names, so similarity restores them. Mutates feats in place."""
    targets = [f for f in feats
               if f.gene_type == "tRNA" and getattr(f, "has_intron", False)
               and getattr(f, "exons", None)]
    if not targets or not Path(trna_db + ".nhr").exists():
        return
    qfa = None
    try:
        with tempfile.NamedTemporaryFile(
            suffix=".fasta", mode="w", delete=False
        ) as fh:
            for i, f in enumerate(targets):
                s = _spliced_seq(genome_seq, f.exons, f.strand)
                if s:
                    fh.write(f">{i}\n{s}\n")
            qfa = fh.name
        result = subprocess.run([
            "blastn", "-query", qfa, "-db", trna_db,
            "-outfmt", "6 qseqid sseqid bitscore",
            "-word_size", "7", "-dust", "no", "-max_target_seqs", "5",
        ], capture_output=True, text=True, timeout=60)
    except (subprocess.SubprocessError, OSError):
        return
    finally:
        if qfa and os.path.exists(qfa):
            os.unlink(qfa)

    best = {}  # qid → (bitscore, name)
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        p = line.split("\t")
        if len(p) < 3:
            continue
        qid = int(p[0])
        name = _extract_gene_from_sid(p[1])
        if not name:
            continue
        bs = float(p[2])
        if qid not in best or bs > best[qid][0]:
            best[qid] = (bs, name)

    for i, f in enumerate(targets):
        hit = best.get(i)
        if hit and hit[1] != f.gene_name:
            f.notes.append(f"renamed {f.gene_name}→{hit[1]} by BLAST")
            f.gene_name = hit[1]
