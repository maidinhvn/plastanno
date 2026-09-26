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
import hashlib
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import List, Dict

from Bio.Seq import Seq

from ..core.feature import Feature
from ..core import trna_identity as _TI
from ..core import coords as _coords


# ── HMM-based CDS detection ───────────────────────────────────────────────────

def run_hmm_scan(genome_seq, hmm_dir, threads=4, gene_catalog=None):
    # gene_catalog is accepted for interface symmetry with the other detectors and
    # is not read here: the scan appends the whole sequence rather than a length
    # derived from the longest catalogued gene. See WRAP below.
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
    # The genome is circular, but a six-frame translation of a linear string cuts
    # every ORF that crosses position 1 in two. accD in a genome deposited with the
    # origin inside it came back 255 bp short here while Engine A had it whole, so
    # the two engines disagreed about a gene neither had got wrong. Scanning the
    # sequence with its own head appended makes such an ORF contiguous; the extra
    # copy of a gene that falls in the appended part maps to the same coordinates
    # and is removed by the duplicate check below.
    #
    # The whole sequence is appended, not a fixed tail. A shorter overlap is enough
    # to stop an origin-crossing ORF being split, but it is NOT enough to make the
    # SCORE origin-independent: a gene inside the appended region gets scanned
    # twice, once in each context, and the better of the two wins, while a gene
    # further along is scanned once. psbE scored 0.921 in one presentation of a
    # genome and 0.929 in a rotated copy of the same genome for exactly that
    # reason. Doubling gives every gene the same two contexts, and the score comes
    # out the same wherever the origin falls: measured 0.9210 in both frames, at a
    # cost of about one second on a 155 kb genome.
    WRAP = L
    fwd_ext = fwd + fwd[:WRAP]
    Lx = len(fwd_ext)
    rev = str(Seq(fwd_ext).reverse_complement())
    MIN_AA = 20

    targets = []  # (name, protein_segment)
    for strand_char, seqstr in (("f", fwd_ext), ("r", rev)):
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
            gs, ge, strand = _map_to_genome(frame_id, int(seg_off),
                                            env_from, env_to, Lx)
            if ge <= gs or gs < 0 or (ge - gs) > L:
                continue
            if gs >= L:                 # wholly inside the appended copy
                gs, ge = gs - L, ge - L
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
            arcs = _coords.region_arcs(gs % L, ge % L if ge != L else L, L) \
                if ge > L else [(gs, ge)]
            # Overlap on the circle: the same locus reached once directly and once
            # through the appended copy must collapse to one feature, and two IR
            # copies must not.
            if any(_coords._inter_bp(arcs, k) > 0 for k in kept):
                continue
            kept.append(arcs)
            f = Feature(gene_name=gene, gene_type="CDS",
                        start=gs % L, end=(ge % L if ge != L else L),
                        strand=strand, engine="B", s_model=sm)
            if ge > L:                       # the ORF crosses the origin
                f.exons = arcs
            features.append(f)
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

# The production ARAGORN invocation, as one constant rather than a literal at
# the call site. `-c` is not a preference: a plastome is circular, and in linear
# mode (-l) ARAGORN reports a tRNA crossing the origin with a negative start,
# e.g. "c[-3,71]" for the trnH-GUG of NC_073016.1. `aragorn_coords` refuses such
# a coordinate, so dropping `-c` now fails at the call site instead of quietly
# losing every origin-crossing tRNA.
ARAGORN_ARGS = ("-i", "-t", "-c", "-w")


def check_aragorn_args(args=ARAGORN_ARGS):
    """Raise unless the invocation puts ARAGORN in circular mode.

    Separate from run_aragorn so the lock can be tested without a subprocess.
    """
    args = tuple(args)
    if "-c" not in args:
        raise RuntimeError(
            "ARAGORN must run in circular mode: '-c' missing from %r. In linear "
            "mode an origin-crossing tRNA is reported with a negative start, "
            "which the coordinate parser rejects." % (args,))
    if "-l" in args:
        raise RuntimeError(
            "ARAGORN linear mode ('-l') is not supported: %r. A plastome is "
            "circular; use '-c'." % (args,))
    return args


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


def project_to_circle(linear_exons, glen):
    """Map exons from an UNROLLED transcript axis back onto the circle.

    `linear_exons` are 0-based half-open and may run past `glen`, because a
    feature crossing the origin is unrolled to s1 .. e1+glen before its intron is
    split out. An exon that straddles the origin becomes TWO storage arcs, which
    is what GenBank's join() expresses too; the biological exon count is
    recoverable because the two arcs meet exactly at the seam.

    Returns arcs ascending; transcription order is restored by
    core.coords.transcript_parts, which knows how to read a wrapped feature.
    """
    out = []
    for a, b in linear_exons:
        if b <= glen:
            out.append((a, b))
        elif a >= glen:
            out.append((a - glen, b - glen))
        else:
            out.append((a, glen))          # up to the origin
            out.append((0, b - glen))      # and past it
    return sorted(x for x in out if x[1] > x[0])


def exons_from_aragorn(s1, e1, strand, i_off, i_len, glen):
    """(arcs, exons_or_None, wrapped) for one ARAGORN call.

    THE single decision about how an ARAGORN line becomes internal coordinates.
    Both run_aragorn and the coordinate-translation gate call this, so the gate
    cannot pass by re-implementing the branching it is meant to check.

    An origin-crossing feature is unrolled onto a linear transcript axis before
    its intron is split out, then projected back. The intron does not disappear
    because the genome was rotated.
    """
    arcs, wrapped = aragorn_coords(str(s1), str(e1), glen)
    exons = None
    if i_off is not None:
        if wrapped:
            lin = _aragorn_intron_exons(int(s1), int(e1) + glen, strand,
                                        int(i_off), int(i_len))
            exons = project_to_circle(lin, glen) if lin else None
        else:
            exons = _aragorn_intron_exons(int(s1), int(e1), strand,
                                          int(i_off), int(i_len))
    return arcs, exons, wrapped


def aragorn_coords(s_str, e_str, glen):
    """ARAGORN's 1-based inclusive [start,end] -> (arcs, wrapped).

    `arcs` is a list of 0-based half-open intervals: one for an ordinary feature,
    two for one crossing the origin, which ARAGORN writes as [high, low]. Split out
    so the conversion has a single definition and can be tested directly; getting it
    wrong shifted every non-intron tRNA by one base, and swapping the two numbers of
    a wrapped feature inflated it to nearly the whole genome.
    """
    s1, e1 = int(s_str), int(e_str)
    # Fail loudly rather than translate nonsense. Under -c every coordinate is a
    # real position on the circle, so anything outside 1..glen means the caller
    # is not in circular mode (linear mode writes an origin-crossing feature as
    # "c[-3,71]") or is passing coordinates from a different sequence. Silently
    # accepting a negative start produced an arc at a negative index, which then
    # sliced the sequence from the wrong end without any error.
    if not (1 <= s1 <= glen) or not (1 <= e1 <= glen):
        raise ValueError(
            "ARAGORN coordinate outside 1..%d: [%s,%s]. A zero or negative "
            "position is what linear mode (-l) emits for a feature crossing the "
            "origin; production runs %r." % (glen, s_str, e_str, ARAGORN_ARGS))
    if e1 < s1:
        return [(s1 - 1, glen), (0, e1)], True
    return [(s1 - 1, e1)], False


def parse_aragorn_output(text, glen):
    """Turn ARAGORN's batch (-w) stdout into Features. Split out from the subprocess
    call so the whole parse path — coordinates, strand, exons, intron flag and
    provenance note — can be tested on raw lines rather than only end to end."""
    AA3TO1 = {
        "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
        "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
        "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
        "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
        "Sec": "U",
    }
    line_re = re.compile(
        r"^\s*\d+\s+"                  # index
        r"(?:tRNA|tmRNA)-(\w+?)\d*\s+" # tRNA-<AA3> (drop trailing digit)
        r"(c?)\[(\d+),(\d+)\]\s+"      # [c][start,end]
        r"\d+\s+"                       # length
        r"\(([a-zA-Z]+)\)"              # (anticodon)
        r"(?:i\((\d+),(\d+)\))?"        # optional intron i(offset,length)
    )
    # A line that is shaped like a hit but does not satisfy `line_re` used to be
    # skipped in the same breath as ARAGORN's header and blank lines, so a
    # coordinate form the parser cannot read — the negative start of linear mode
    # among them — cost a feature with no diagnostic. Recognise the shape loosely
    # and refuse it explicitly. Header ('>NC_...'), the 'N genes found' line and
    # blanks do not match this, so they still pass through untouched.
    hit_shaped = re.compile(r"^\s*\d+\s+(?:tRNA|tmRNA)\S*\s+c?\[")
    # ARAGORN writes 'tRNA-???' when it finds a tRNA structure but cannot assign
    # an amino acid to it. That is ordinary output, not corruption, and it was
    # reaching the raise below: one such line aborted the whole genome. Three of
    # sixty genomes died this way in benchmark_v3/intron_glocal/baseline_exits.txt
    # (NC_034893.1, NC_034864.1, NC_077615.1). An unassignable call costs one
    # candidate, not an annotation -- but it is still announced, because the
    # point of the raise was that nothing should vanish silently.
    unknown_aa = re.compile(r"^\s*\d+\s+(?:tRNA|tmRNA)-\?+\s+c?\[")
    features = []
    for line in text.splitlines():
        m = line_re.match(line)
        if not m:
            if unknown_aa.match(line):
                print("        [warn] ARAGORN could not assign an amino acid; "
                      "candidate skipped: %s" % line.strip())
                continue
            if hit_shaped.match(line):
                raise ValueError(
                    "unparsable ARAGORN hit line: %r. Expected "
                    "'<n> tRNA-<AA> [c][start,end] <len> (anticodon)[i(off,len)]' "
                    "with 1-based positions inside the genome." % line.strip())
            continue
        aa3, comp, s_str, e_str, anti, i_off, i_len = m.groups()
        # ARAGORN reports 1-based inclusive [start,end]; this module returns
        # 0-based half-open [start,end). The start therefore needs -1, which was
        # missing: every non-intron tRNA was annotated one base short at its 5'
        # end, verified against raw ARAGORN output and the reference records.
        # (The intron branch below builds its exons separately and was correct.)
        # ARAGORN reports 1-based inclusive positions. Keep them as `s1, e1` and
        # convert once, where each consumer needs it: `_aragorn_intron_exons`
        # takes the 1-based pair and does its own conversion, so handing it a
        # value that has already been decremented subtracts twice.
        s1, e1 = int(s_str), int(e_str)
        arcs, wrapped = aragorn_coords(s_str, e_str, glen)
        start, end = arcs[0][0], arcs[-1][1]
        strand = -1 if comp == "c" else 1
        anti_rna = anti.upper().replace("T", "U")

        if aa3 in ("fMet", "fMeth"):
            code = "fM"
        else:
            code = AA3TO1.get(aa3.capitalize())
        gene_name = ("trn%s-%s" % (code, anti_rna)) if code else ("trn?-%s" % anti_rna)

        # one decision, shared with the coordinate-translation gate
        arcs, exons, wrapped = exons_from_aragorn(
            s1, e1, strand, i_off, i_len, glen)

        # A wrapped feature must keep start > end, which is how coords.is_wrapped
        # recognises it and how transcript_parts splits tail from head. Taking
        # min/max of the projected exons would silently unwrap it.
        if wrapped:
            f_start, f_end = arcs[0][0], arcs[-1][1]
        elif exons:
            f_start, f_end = exons[0][0], exons[-1][1]
        else:
            f_start, f_end = arcs[0][0], arcs[-1][1]
        f = Feature(
            gene_name=gene_name,
            gene_type="tRNA",
            product="tRNA-%s" % aa3,
            start=f_start,
            end=f_end,
            strand=strand,
        )
        f.engine = "B"
        if exons:
            f.exons = exons
            f.has_intron = True
            f.notes.append("ARAGORN intron i(%s,%s)" % (i_off, i_len))
        elif wrapped:
            # both arcs go on the feature; without them a wrapped tRNA is
            # described only by a start/end pair that runs backwards
            f.exons = list(arcs)
            f.notes.append("ARAGORN (crosses the origin)")
        else:
            f.notes.append("ARAGORN")
        features.append(f)
    return features


def run_aragorn(genome_seq, genome_len):
    glen = genome_len
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

        # -c, not -l: a plastome is circular, and in linear mode ARAGORN reports a
        # tRNA crossing the origin with a negative start ("c[-3,71]") that the
        # line parser cannot read, so the feature was dropped without a word.
        # The arguments live in ARAGORN_ARGS and are checked here, so an edit
        # that removes circular mode stops the run instead of thinning the
        # annotation by one tRNA per genome.
        result = subprocess.run(
            ["aragorn"] + list(check_aragorn_args()) + [tmp.name],
            capture_output=True, text=True, timeout=120,
        )

        if result.returncode != 0:
            raise RuntimeError("aragorn exited %d: %s"
                               % (result.returncode, (result.stderr or "")[:200]))
        parsed = parse_aragorn_output(result.stdout, glen)
        # ARAGORN states how many genes it found; if the parser recovered fewer, a
        # line was silently discarded. That is how every origin-crossing tRNA was
        # lost for as long as the tool ran in linear mode, so it must be loud.
        m = re.search(r"^\s*(\d+)\s+genes found", result.stdout, re.M)
        if m and int(m.group(1)) != len(parsed):
            print("      WARNING: ARAGORN reported %s tRNA but %d were parsed"
                  % (m.group(1), len(parsed)))
        features.extend(parsed)
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
        # A timeout is NOT an absence of tRNAs. Returning the partial list in
        # silence makes a slow genome indistinguishable from a clean negative,
        # and the caller has no way to tell that the evidence is incomplete.
        print("      WARNING: tRNA BLAST timed out; %d candidate(s) kept, the "
              "tRNA evidence for this genome is INCOMPLETE" % len(features))
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

# Approximate expected intron length (bp) for the long-intron plastid tRNAs. Used
# ONLY as a conservative rescue in the exon pairing: the default "tightest intron"
# rule can mis-pair a 5' exon with a spurious intermediate 3' exon hit inside the
# real intron, truncating the tRNA (e.g. trnA-UGC annotated ~453 bp vs the real
# ~880 bp). When the tightest pick is clearly SHORT for such a gene, prefer a
# plausible full-length pair. Genes not listed are never rescued -> unchanged.
_INTRON_EXPECTED = {
    "trnA-UGC": 800, "trnI-GAU": 950, "trnG-UCC": 700,
    "trnK-UUU": 2500, "trnV-UAC": 570, "trnL-UAA": 500,
}



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
        # Same reasoning as blast_trna: an intron-bearing tRNA that is missed
        # because the search ran out of time must not look like one that is not
        # there.
        print("      WARNING: intron-tRNA BLAST timed out; %d hit(s) kept, the "
              "intron-tRNA evidence for this genome is INCOMPLETE" % len(raw_hits))
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
        cands = []  # all valid (gap, -combined_pident, left_exon, right_exon) pairs
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
                cands.append((gap, -(a["pident"] + b["pident"]), la, rb))
        if not cands:
            continue
        best = min(cands)              # default: tightest intron (unchanged)
        # Conservative rescue: if the tightest pick is clearly TRUNCATED for a
        # long-intron gene (a spurious intermediate 3' exon was chosen inside the
        # real intron), prefer a plausible full-length pair. Fires ONLY when the
        # pick is < 0.6x the gene's expected intron — so normal, correct picks
        # (>= 0.6x expected) are never touched.
        exp = _INTRON_EXPECTED.get(gene)
        if exp and best[0] < 0.6 * exp:
            full = [c for c in cands if 0.6 * exp <= c[0] <= 1.5 * exp]
            if full:
                # c[2] and c[3] are GENOMIC exon coordinates, which a rotation or
                # a strand flip changes; ordering by them is arbitrary between two
                # presentations of one genome. The exon sequences are not.
                def _pair_digest(c):
                    seq = "".join(genome_seq[x:y] for x, y in (c[2], c[3]))
                    return hashlib.sha256(seq.encode()).hexdigest()
                best = min(full, key=lambda c: (c[1], abs(c[0] - exp),
                                                _pair_digest(c)))
        _, _, la, rb = best
        features.append(Feature(
            gene_name=gene, gene_type="tRNA",
            start=la[0], end=rb[1], strand=strand,
            exons=[la, rb], has_intron=True, engine="B",
            s_model=1.0,
        ))

    # Deduplicate per locus (two IR copies stay separate; they don't overlap).
    # Keyed on actual overlap, not on a 200-bp bucket of the start coordinate: two
    # calls of one tRNA at 399 and 401 fall either side of a bucket edge and both
    # used to survive as separate genes, while two genuinely distinct copies inside
    # one bucket were collapsed into one.
    glen = len(genome_seq) or None
    deduped = []
    for f in sorted(features, key=lambda x: x.s_model, reverse=True):
        if any(k.gene_name == f.gene_name and k.strand == f.strand
               and _coords.locus_overlap_fraction(k, f, glen) > 0
               for k in deduped):
            continue
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
                # The guard above already established this line is hit-shaped, so
                # a parse failure here means a real hit was dropped. The ARAGORN
                # parser was hardened for exactly this case; keep the two
                # consistent rather than skipping in silence.
                print("      WARNING: unreadable tRNAscan-SE hit line skipped: %r"
                      % line.rstrip()[:120])
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
    exon_mode = "aragorn",
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
    hmm_feats = run_hmm_scan(genome_seq, hmm_dir, threads, gene_catalog)
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

    # 2d. tRNA geometry: legacy source-priority merge, or the hybrid rule table.
    #
    # hybrid takes intron-free boundaries from tRNAscan-SE and everything else
    # from ARAGORN/exon-DB. It is fail-closed: a missing or failing tRNAscan-SE
    # raises rather than returning [], because silently falling back would emit
    # legacy coordinates under a hybrid label.
    # tRNA geometry is NOT decided here any more.
    #
    # The hybrid boundary swap used to run at this point, inside Engine B. That
    # placed it UPSTREAM of reconcile._select, which decides the inventory by
    # clustering same-name features that overlap at all. Moving two overlapping
    # duplicates onto their tRNAscan-SE boundaries left them abutting with zero
    # overlap, so they stopped clustering and both survived: a boundary change
    # silently became an inventory change (NC_026958.1 gained a second
    # trnT-GGU). Engine B therefore always produces the legacy inventory, and
    # the swap happens after reconciliation -- see pipeline.run.
    from .trna_hybrid import tag_identity_sources
    # Stamp where each NAME came from before the merge collapses the three
    # lists; afterwards the surviving feature no longer records which caller
    # won. Inert unless hybrid reads it later, so legacy output is unaffected.
    tag_identity_sources(aragorn_feats, blast_feats, intron_feats)

    trnascan_feats = run_trnascan(genome_seq) if use_trnascan else []
    if use_trnascan:
        print(f"        tRNAscan-SE: {len(trnascan_feats)} tRNAs")
    # Merge tRNA: intron > tRNAscan-SE > ARAGORN > BLAST
    trna_feats = _merge_trna_sources(
        aragorn_feats, blast_feats, intron_feats, trnascan_feats,
        glen=len(genome_seq), exon_mode=exon_mode,
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


def _merge_trna_sources(aragorn, blast, intron, trnascan=None, glen=None,
                        exon_mode="aragorn"):
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
    # ARAGORN's own intron call, ranked ABOVE the exon-DB call. See below.
    SRC_ARAGORN_INTRON = 3.5
    SRC_NAMES = {3.5: "ARAGORN-intron", 3: "intron-BLAST", 2.5: "tRNAscan-SE",
                 2: "ARAGORN", 1: "BLAST"}

    def _aragorn_rank(f):
        """Structure places boundaries; similarity detects.

        The exon-DB path has always outranked ARAGORN (3 against 2), so at an
        intron-bearing locus a BLAST similarity call beats a structural one --
        and measured over all 238 such loci in the development set, the
        similarity call is far worse: acceptor stem 4.35/7 with 2 exon pairs
        exactly matching the reference, against ARAGORN's 6.51/7 with 44. The
        reference itself scores 6.58/7, so ARAGORN reaches the quality of the
        annotations it is being judged against while the exon DB sits near the
        2.74/7 random-sequence null.

        Only INTRON-BEARING ARAGORN candidates are promoted. A single-exon
        ARAGORN call keeps its usual rank, so nothing outside this locus class
        is touched. ARAGORN had a call at 238 of 238 intron loci, so the exon DB
        is never needed as a geometry fallback in practice -- but it keeps its
        DETECTION role, and where ARAGORN has no candidate the exon-DB call
        still wins by default.
        """
        if exon_mode != "aragorn":
            return SRC_ARAGORN
        arcs = f.exons if getattr(f, "exons", None) else [(f.start, f.end)]
        return (SRC_ARAGORN_INTRON
                if _coords.biological_exon_count(sorted(arcs), glen) > 1
                else SRC_ARAGORN)

    cand = ([(f, SRC_INTRON)   for f in intron]
            + [(f, SRC_TRNASCAN) for f in (trnascan or [])]
            + [(f, _aragorn_rank(f)) for f in aragorn]
            + [(f, SRC_BLAST)    for f in blast])

    def same_locus(a, b):
        # strand matters: the two inverted-repeat copies of a tRNA face opposite
        # ways, and merging them would delete one. Extents come from core.coords
        # so a feature crossing the origin is measured over its arcs.
        if a.strand != b.strand:
            return False
        inter   = _coords.overlap_bp(a, b)
        shorter = min(_coords.spliced_length(a), _coords.spliced_length(b))
        if shorter <= 0:
            return False
        return shorter > 0 and inter / shorter > 0.5

    def rank(item):
        f, pr = item
        return (pr, round(getattr(f, "s_model", 0.0) or 0.0, 3),
                _coords.spliced_length(f))

    cand.sort(key=rank, reverse=True)   # best source/score/length first
    kept = []
    for f, pr in cand:
        rival = next((k for k in kept if same_locus(f, k)), None)
        if rival is None:
            kept.append(f)
            continue
        # Another source called the same locus. Source priority decides which model
        # is written, but when the two disagree about where the feature begins or
        # ends that choice is not evidence — it is a tie broken by convention. Keep
        # the discarded model in the provenance and mark the survivor for review, so
        # a boundary conflict is visible rather than resolved silently. This is how
        # the wrapped ARAGORN trnH and the linear BLAST trnH differ.
        # Record the discarded model in full. Whether it is worth a reviewer's
        # attention is decided later, by finalize_qc, against a configurable
        # threshold — the merge itself must not lose the alternative.
        if ((rival.start, rival.end) != (f.start, f.end)
                or rival.gene_name != f.gene_name or rival.strand != f.strand):
            rival.alternatives.append({
                "caller":   SRC_NAMES.get(pr, "unknown"),
                "gene_name": f.gene_name,
                "type":     f.gene_type,
                "strand":   f.strand,
                "start":    f.start,
                "end":      f.end,
                "exons":    [list(x) for x in (f.exons or [])],
                "wrapped":  f.start > f.end,
                "score":    round(getattr(f, "s_model", 0.0) or 0.0, 3),
                # circular distance: a 2 bp disagreement across the origin is 2 bp,
                # not the width of the genome
                "distance_bp": _coords.boundary_distance(rival, f, glen),
            })
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
        if not hit:
            continue
        # Route through the single identity API. Assigning gene_name here and
        # leaving product alone is what produced 121 self-contradictory tRNA
        # records (gene=trnG-UCC beside product=tRNA-Ser).
        canon, note = _TI.resolve(_TI.from_gene_name(f.gene_name),
                                  _TI.from_gene_name(hit[1]), _TI.BLAST_TRNA_DB,
                                  detector_symbol=f.gene_name)
        if _TI.gene_name(canon) != f.gene_name or _TI.product(canon) != f.product:
            _TI.apply_trna_identity(f, canon, note)
