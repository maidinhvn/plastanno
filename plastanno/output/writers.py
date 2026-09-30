"""
Output writers for Plastanno v2.
Writes: .gb, .gff3, .faa, .ffn, .frn, .report, and (unless no_plot) a circular
plastome map (_map.png/.pdf/.svg).

Improvements over v1:
- Provenance in GenBank notes
- Confidence flags in GFF3
- Clean FASTA headers
"""
import re
import hashlib
from pathlib import Path
from typing import List, Dict
from datetime import datetime
from ..core.feature import Feature
from ..core import coords as _coords


# Written to /organism when the user gives no --organism. It is a stand-in, not a
# real determination, so write_all() warns whenever it is used.
PLACEHOLDER_ORGANISM = "Viridiplantae"


def _strand_str(strand):
    if strand == 1  or strand == "+": return "+"
    if strand == -1 or strand == "-": return "-"
    return "."


def _sanitize_id(name):
    """Identifier-safe form of a sequence name: keep the part before the first
    comma or whitespace (the FASTA/GenBank convention), and map any remaining
    invalid character to '_'. Idempotent — a clean accession such as
    'NC_053537.1' is returned unchanged. A FASTA header like
    'Centella_asiatica_CA1_plastid, complete genome' becomes
    'Centella_asiatica_CA1_plastid'."""
    base = re.split(r"[,\s]", str(name).strip(), 1)[0]
    base = re.sub(r"[^A-Za-z0-9_.\-]", "_", base).strip("_")
    return base or "seq"


def _locus_name(name, maxlen=16):
    """A spec-valid GenBank LOCUS name: sanitized and at most `maxlen` characters
    (the GenBank locus field is columns 13-28). A long name's distinguishing part
    is often its suffix (e.g. a sample tag CA1/CA2/CA18), so a blind truncation
    would make several genomes collide on one LOCUS name; we therefore keep a
    prefix plus a short hash of the full name to stay short AND unique."""
    base = _sanitize_id(name)
    if len(base) <= maxlen:
        return base
    h = hashlib.md5(str(name).encode()).hexdigest()[:4]
    return f"{base[:maxlen - 5]}_{h}"


# Column order of the tRNA source-merge ledger, shared with output.verify so the
# checker validates the real schema rather than "at least one column".
ALT_COLUMNS = ["kept_gene", "kept_type", "kept_strand", "kept_start", "kept_end",
               "alt_caller", "alt_gene", "alt_strand", "alt_start", "alt_end",
               "alt_exons", "alt_wrapped", "alt_score", "distance_bp"]


def _exon_parts(ann, genome_len=None):
    """Exons of a feature in transcript order, each with its own strand.

    Thin adapter over core.coords.transcript_parts — the single definition of what
    transcription order means, including for a feature crossing the origin, which
    the writers used to order wrongly. Kept in the ((start, end), strand) shape the
    GenBank, GFF3 and feature-table writers already consume, so all three formats
    and the extracted sequence describe one exon structure.
    """
    return [((s, e), st) for s, e, st in _coords.transcript_parts(ann, genome_len)]


def _is_trans_spliced(ann, parts):
    """True for a trans-spliced gene: either the feature says so, or its exons do
    not all lie on one strand."""
    return bool(getattr(ann, "is_trans_spliced", False)) or len({st for _, st in parts}) > 1


def write_genbank(annotations, genome_seq, accession,
                   genome_len, ir_boundaries, relatives,
                   out_path, organism=None):
    """Write GenBank format with provenance notes.

    `organism` is the taxon name (from --organism). It is written consistently to
    the ORGANISM/SOURCE header, the DEFINITION line and the source feature's
    /organism qualifier, so downstream tools (NCBI submission, CPJSdraw, the
    circular map) all see the same taxon. When it is not supplied we fall back to
    the generic placeholder and the caller warns — a silent placeholder previously
    ended up labelling every sample "Viridiplantae".
    """
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio.SeqFeature import (SeqFeature, FeatureLocation, CompoundLocation,
                                BeforePosition, AfterPosition)

    org = (organism or "").strip() or PLACEHOLDER_ORGANISM
    # No trailing '.' — Biopython appends one when writing the DEFINITION line.
    definition = (f"{org} chloroplast, complete genome"
                  if organism else "Annotated by Plastanno v2")

    seq    = Seq(genome_seq)
    record = SeqRecord(
        seq,
        id          = _sanitize_id(accession),
        name        = _locus_name(accession),   # valid GenBank LOCUS (<=16 chars, no comma)
        description = definition,
    )
    record.annotations["molecule_type"]      = "DNA"
    record.annotations["topology"]           = "circular"
    # Division: PLN covers plants / algae / fungi, i.e. the entire scope of a
    # plastome annotator. Without this Biopython writes the invalid default 'UNK'.
    record.annotations["data_file_division"] = "PLN"
    # Without this the SOURCE / ORGANISM header lines are written as a bare '.'
    record.annotations["organism"]           = org
    record.annotations["source"]             = org
    record.annotations["date"]          = \
        datetime.now().strftime("%d-%b-%Y").upper()

    # Source feature. /organelle is required for an organellar GenBank submission.
    record.features.append(SeqFeature(
        FeatureLocation(0, genome_len),
        type = "source",
        qualifiers = {
            "organism" : [org],
            "mol_type" : ["genomic DNA"],
            "organelle": ["plastid:chloroplast"],
        }
    ))

    # Add annotations
    for ann in sorted(annotations, key=lambda x: x.start):
        # Build exon locations. Exons are spliced out, so a multi-exon gene gets a
        # join() and not its outer span; a trans-spliced gene keeps the given
        # transcript order and per-exon strand (do not re-sort).
        parts = _exon_parts(ann, genome_len)
        locs  = [FeatureLocation(s, e, st) for (s, e), st in parts]
        loc   = locs[0] if len(locs) == 1 else CompoundLocation(locs)
        trans = _is_trans_spliced(ann, parts)

        # Provenance notes (v2 improvement)
        note = (f"engine={ann.engine}; "
                f"confidence={ann.confidence:.2f}; "
                f"flag={ann.flag}")
        if ann.notes:
            note += "; " + "; ".join(ann.notes)

        gene_q = {"gene": [ann.gene_name]}
        if trans:
            # A join() whose segments lie on different strands is only valid with
            # /trans_splicing; NCBI rejects it otherwise.
            gene_q["trans_splicing"] = [None]

        if ann.is_pseudogene:
            # NCBI convention for a truncated remnant such as the IR copy of ycf1:
            # a gene feature carrying /pseudo and NO coding feature — a pseudogene
            # has no valid CDS or translation to report.
            gene_q["pseudo"] = [None]
            gene_q["note"]   = [note + f"; pseudogene_reason={ann.pseudogene_reason}"]
            pl = loc if (trans or ann.start > ann.end) else FeatureLocation(
                ann.start, ann.end, ann.strand)
            record.features.append(SeqFeature(pl, type="gene", qualifiers=gene_q))
            continue

        # A tRNA product must name the RNA, never the gene symbol. This fallback
        # is where 112 records acquired /product="trnI-GAU": the intron-tRNA path
        # builds a Feature with no product and the gene name was substituted.
        _prod = ann.product or (
            "tRNA-OTHER" if ann.gene_type == "tRNA" else ann.gene_name)
        qualifiers = {
            "gene"   : [ann.gene_name],
            "product": [_prod],
            "note"   : [note],
        }
        if trans:
            qualifiers["trans_splicing"] = [None]

        # Add protein sequence for CDS
        if getattr(ann, "rna_edited_start", False):
            qualifiers["exception"] = ["RNA editing"]
        if ann.gene_type == "CDS" and ann.protein:
            qualifiers["translation"] = [ann.protein]
            qualifiers["codon_start"] = ["1"]
        elif ann.gene_type == "CDS" and getattr(ann, "orf_incomplete", False):
            # Unresolved boundary: mark both ends partial ('<'/'>') so the length is
            # legally not a multiple of 3, and give codon_start without asserting a
            # translation. The coordinates are the ones we called — nothing moved.
            qualifiers["codon_start"] = ["1"]
            # '<' goes before the lowest coordinate and '>' after the highest,
            # independent of strand (Biopython renders the complement form).
            i_lo = min(range(len(locs)), key=lambda i: int(locs[i].start))
            i_hi = max(range(len(locs)), key=lambda i: int(locs[i].end))
            locs = [FeatureLocation(BeforePosition(int(l.start)) if i == i_lo else l.start,
                                    AfterPosition(int(l.end))   if i == i_hi else l.end,
                                    l.strand)
                    for i, l in enumerate(locs)]
            loc = locs[0] if len(locs) == 1 else CompoundLocation(locs)

        # Gene feature. A gene is ONE interval spanning its introns — only the
        # coding feature is spliced. The reference records are unambiguous
        # (rpoC1/clpP/ndhB/petB/trnK-UUU all have a single-interval gene over a
        # multi-interval CDS), and emitting a join() here draws NCBI's
        # MultiIntervalGene warning. The exception is a trans-spliced gene, whose
        # pieces are transcribed separately: there the reference does use a join,
        # with /trans_splicing (rps12).
        # A feature crossing the origin has start > end, which FeatureLocation
        # rejects; it must keep the two-part location its exons already describe,
        # exactly as a trans-spliced gene does.
        wrapped = ann.start > ann.end
        gene_loc = loc if (trans or wrapped) else FeatureLocation(ann.start, ann.end,
                                                                  ann.strand)
        record.features.append(SeqFeature(
            gene_loc, type="gene", qualifiers=gene_q
        ))

        # CDS/tRNA/rRNA feature
        record.features.append(SeqFeature(
            loc, type=ann.gene_type, qualifiers=qualifiers
        ))

    SeqIO.write(record, out_path, "genbank")
    return record


def write_gff3(annotations, genome_len, accession, out_path):
    """Write GFF3 format with confidence flags."""
    # Feature IDs must be unique — GFF3 parsers treat repeated IDs as one feature.
    # "<gene>_<start>" alone is not: the two trans-spliced rps12 copies share their
    # 5' exon, so both copies produced ID=rps12_<same start> and a parser merged
    # them into a single six-exon gene. Disambiguate a repeat with a numeric
    # suffix; genes that are already unique keep the familiar form.
    used = {}
    def _uid(ann):
        base = f"{ann.gene_name}_{ann.start}"
        used[base] = used.get(base, 0) + 1
        return base if used[base] == 1 else f"{base}.{used[base]}"

    with open(out_path, "w") as f:
        f.write("##gff-version 3\n")
        f.write(f"##sequence-region {accession} 1 {genome_len}\n")

        for ann in sorted(annotations, key=lambda x: x.start):
            strand = _strand_str(ann.strand)
            gid    = _uid(ann)
            attrs  = (f"ID={gid};"
                      f"Name={ann.gene_name};"
                      f"engine={ann.engine};"
                      f"confidence={ann.confidence:.2f};"
                      f"flag={ann.flag}")
            if ann.is_pseudogene:
                attrs += ";pseudo=true"

            # Gene line — the outer span, introns included (this is what a gene
            # feature means in GFF3).
            # A feature crossing the origin cannot be one GFF3 line: the spec has
            # no way to write start > end. Split the outer span at the origin and
            # emit one line per arc, sharing an ID.
            #
            # This used to read ann.exons, which was wrong twice over. A wrapped
            # tRNA with no explicit exons — every intron-free one that
            # --trna-mode hybrid re-seats across the origin — produced an EMPTY
            # list and therefore no gene line at all, so verify.py could not
            # resolve Parent= back to a Name= and read the feature's own ID as its
            # gene name (trnH-GUG_154576 against the .gb's trnH-GUG). And a
            # wrapped multi-exon gene produced one gene line per exon with the
            # introns cut out, which is not what a gene line means.
            if ann.start > ann.end:
                gspans = [(ann.start, genome_len)]
                if ann.end > 0:
                    gspans.append((0, ann.end))
            else:
                gspans = [(ann.start, ann.end)]
            for gs, ge in gspans:
                f.write("\t".join([
                    accession, "Plastanno", "gene",
                    str(gs+1), str(ge),
                    ".", strand, ".", attrs
                ]) + "\n")

            if ann.is_pseudogene:
                # No coding feature for a pseudogene, matching the GenBank writer.
                continue

            # Coding/RNA lines — ONE PER EXON. Emitting a single line spanning
            # ann.start..ann.end made every intron-containing gene look like an
            # uninterrupted coding sequence, and reported the trans-spliced rps12
            # as a ~70 kb CDS (its three exons are split across the LSC and the
            # IR). Column 8 carries the CDS phase, as the GFF3 spec requires.
            parts     = _exon_parts(ann, genome_len)
            feat_type = ann.gene_type
            parent    = gid
            trans     = ";trans_splicing=true" if _is_trans_spliced(ann, parts) else ""
            coding    = 0
            for (es, ee), est in parts:
                phase = str((3 - coding % 3) % 3) if feat_type == "CDS" else "."
                f.write("\t".join([
                    accession, "Plastanno", feat_type,
                    str(es+1), str(ee),
                    f"{ann.confidence:.2f}", _strand_str(est), phase,
                    f"ID={feat_type.lower()}-{parent};"
                    f"Parent={parent};"
                    f"product={ann.product or ann.gene_name}{trans}"
                ]) + "\n")
                coding += ee - es


def write_fasta_files(annotations, genome_seq,
                       accession, out_dir, genome_len=None):
    """Write .faa, .ffn, .frn FASTA files."""
    faa_path = out_dir / f"{accession}.faa"
    ffn_path = out_dir / f"{accession}.ffn"
    frn_path = out_dir / f"{accession}.frn"

    with open(faa_path,"w") as faa, \
         open(ffn_path,"w") as ffn, \
         open(frn_path,"w") as frn:

        for ann in sorted(annotations, key=lambda x: x.start):
            header = (f">{accession}_{ann.gene_name} "
                      f"{ann.start+1}..{ann.end} "
                      f"[{_strand_str(ann.strand)}] "
                      f"engine={ann.engine} "
                      f"flag={ann.flag}")

            if ann.gene_type == "CDS":
                if ann.is_pseudogene:
                    continue          # no coding sequence to report
                if ann.protein:
                    faa.write(f"{header}\n{ann.protein}\n")
                # Extract nucleotide
                seq = _extract_seq(genome_seq, ann, genome_len)
                ffn.write(f"{header}\n{seq}\n")

            elif ann.gene_type in ("tRNA","rRNA"):
                seq = _extract_seq(genome_seq, ann, genome_len)
                frn.write(f"{header}\n{seq}\n")

    return faa_path, ffn_path, frn_path


def _extract_seq(genome_seq, feat, genome_len=None):
    """Coding-orientation sequence of a feature, via the coordinate API.

    This used to concatenate the exons in ascending genomic order and reverse-
    complement the result. That is right for a linear gene and wrong for one
    crossing the origin, whose first transcribed base is not its lowest
    coordinate: the .ffn and .faa records for such a gene began in the middle.
    core.coords.transcript_parts is the one place that decides part order.
    """
    return _coords.extract(genome_seq, feat, genome_len)


# Standard plastome functional gene classification (Sugiura-style groups), used
# for the categorised gene table in the .report. Members are the canonical gene
# set per group; only those actually annotated are shown, and the count reflects
# annotated instances (so IR-duplicated genes contribute twice).
GENE_CATEGORIES = [
    ("Genes for photosynthesis", [
        ("Subunits of ATP synthase", {"atpA", "atpB", "atpE", "atpF", "atpH", "atpI"}),
        ("Subunits of photosystem II", {"psbA", "psbB", "psbC", "psbD", "psbE", "psbF",
            "psbH", "psbI", "psbJ", "psbK", "psbL", "psbM", "psbN", "pbf1", "psbT", "psbZ"}),
        ("Subunits of NADH dehydrogenase", {"ndhA", "ndhB", "ndhC", "ndhD", "ndhE", "ndhF",
            "ndhG", "ndhH", "ndhI", "ndhJ", "ndhK"}),
        ("Subunits of cytochrome b/f complex", {"petA", "petB", "petD", "petG", "petL", "petN"}),
        ("Subunits of photosystem I", {"psaA", "psaB", "psaC", "psaI", "psaJ", "ycf3", "ycf4"}),
        ("Subunit of Rubisco", {"rbcL"}),
    ]),
    ("Self-replication", [
        ("Large subunit of ribosome", {"rpl2", "rpl14", "rpl16", "rpl20", "rpl22", "rpl23",
            "rpl32", "rpl33", "rpl36"}),
        ("DNA-dependent RNA polymerase", {"rpoA", "rpoB", "rpoC1", "rpoC2"}),
        ("Small subunit of ribosome", {"rps2", "rps3", "rps4", "rps7", "rps8", "rps11",
            "rps12", "rps14", "rps15", "rps16", "rps18", "rps19"}),
    ]),
]
OTHER_GENES = {
    "accD": "Subunit of Acetyl-CoA-carboxylase",
    "ccsA": "C-type cytochrome synthesis gene",
    "cemA": "Envelope membrane protein",
    "clpP": "Protease",
    "infA": "Translational initiation factor",
    "matK": "Maturase",
    "ycf1": "Conserved hypothetical ORF",
    "ycf2": "Conserved hypothetical ORF",
    "ycf15": "Conserved hypothetical ORF",
    "lhbA": "Conserved hypothetical ORF",
}


_AA3_BY_LETTER = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln",
    "E": "Glu", "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys",
    "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp",
    "Y": "Tyr", "V": "Val",
}


def _trna_product(gene_name):
    """NCBI product name for a tRNA gene, e.g. trnA-UGC -> 'tRNA-Ala'.

    The gene NAME (trnA-UGC) is not a valid /product: table2asn rejects it with
    MissingTrnaAA because it cannot read the encoded amino acid out of it. The
    product must name the amino acid. trnfM is the formyl-methionine initiator,
    which NCBI writes as tRNA-Met.
    """
    m = re.match(r"^trn([A-Za-z]{1,2})", gene_name or "")
    if not m:
        return None
    tag = m.group(1)
    letter = "M" if tag.lower() in ("fm", "im") else tag[0].upper()
    aa = _AA3_BY_LETTER.get(letter)
    return "tRNA-%s" % aa if aa else None


def write_tbl(annotations, accession, out_path, organism=None, genome_len=None):
    """Write an NCBI 5-column feature table (.tbl).

    A GenBank flatfile is a *distribution* format: NCBI does not accept one as an
    annotation submission. The submission route for an organelle genome is
    table2asn, which takes the sequence (.fsa) plus this feature table and produces
    the .sqn that is actually sent. Without a .tbl the annotation cannot be
    submitted at all, however valid the .gb file is.

    Format: a '>Feature <seqid>' header, then per feature one line of
    'start<TAB>end<TAB>key' followed by one 'start<TAB>end' line per additional
    interval, then '<TAB><TAB><TAB>qualifier<TAB>value' lines. Coordinates are
    1-based and INCLUSIVE, and a minus-strand feature is written with start > end,
    so the intervals are emitted in transcript order. '<' / '>' mark an unresolved
    end, exactly as in the GenBank location.
    """
    lines = [">Feature %s" % _sanitize_id(accession)]
    _seen = {}

    def _locus_tag(ann):
        """A locus_tag unique within the record: PLAS_<gene>[_2] for a second copy."""
        base = re.sub(r"[^A-Za-z0-9]", "", ann.gene_name) or "gene"
        _seen[base] = _seen.get(base, 0) + 1
        n = _seen[base]
        return "PLAS_%s" % base if n == 1 else "PLAS_%s_%d" % (base, n)
    for ann in sorted(annotations, key=lambda x: x.start):
        parts = _exon_parts(ann, genome_len)
        partial = (ann.gene_type == "CDS" and getattr(ann, "orf_incomplete", False))

        def _iv(seg, st, first, last):
            s0, e0 = seg
            a, b = (s0 + 1, e0) if st != -1 else (e0, s0 + 1)
            lo = "<" if (first and partial) else ""
            hi = ">" if (last and partial) else ""
            return "%s%d\t%s%d" % (lo, a, hi, b)

        n = len(parts)
        rows = [_iv(seg, st, i == 0, i == n - 1) for i, (seg, st) in enumerate(parts)]

        # gene feature — the SAME interval structure as the product feature, not the
        # outer span. A trans-spliced rps12 spans ~70 kb end to end, so a single
        # interval made the gene swallow every locus in between (table2asn:
        # "Gene contains 57 other genes"). The GenBank writer already emits a
        # join() here, and the reference records do too.
        trans = _is_trans_spliced(ann, parts)
        if trans or ann.start > ann.end:
            # a trans-spliced gene, or one crossing the origin: the gene feature
            # keeps the same intervals as the product feature. Writing a wrapped
            # gene as a single start..end pair sends it the wrong way round the
            # whole genome.
            lines.append("%s\tgene" % rows[0])
            lines.extend(rows[1:])
        else:
            gs, ge = ann.start + 1, ann.end
            if ann.strand == -1:
                gs, ge = ann.end, ann.start + 1
            lines.append("%s%d\t%s%d\tgene"
                         % ("<" if partial else "", gs, ">" if partial else "", ge))
        lines.append("\t\t\tgene\t%s" % ann.gene_name)
        if trans:
            # a gene whose intervals lie on different strands is only legal with
            # /trans_splicing (table2asn: MixedStrand otherwise)
            lines.append("\t\t\ttrans_splicing")
        # A unique locus_tag is what ties a product feature to ITS gene. Without one
        # the two IR copies of a gene share a name, and table2asn attached each CDS
        # to whichever copy it saw first - on the wrong strand half the time
        # (GeneXrefStrandProblem / CDSgeneRange).
        tag = _locus_tag(ann)
        lines.append("\t\t\tlocus_tag\t%s" % tag)
        if ann.is_pseudogene:
            lines.append("\t\t\tpseudo")
            continue                       # a pseudogene gets no product feature

        lines.append("%s\t%s" % (rows[0], ann.gene_type))
        lines.extend(rows[1:])
        product = ann.product or ann.gene_name
        if ann.gene_type == "tRNA":
            product = _trna_product(ann.gene_name) or product
        lines.append("\t\t\tproduct\t%s" % product)
        lines.append("\t\t\tgene\t%s" % ann.gene_name)
        lines.append("\t\t\tlocus_tag\t%s" % tag)
        if getattr(ann, "rna_edited_start", False):
            # /exception alone is not enough: with no /translation in a .tbl,
            # table2asn translates the CDS itself and an ACG start yields a protein
            # that does not begin with Met (BadProteinStart). /transl_except names
            # the edited codon explicitly, exactly as the reference records do.
            fs, fe = parts[0][0]
            st0 = parts[0][1]
            if st0 == -1:
                pos = "%d..%d" % (fe - 2, fe)
                pos = "complement(%s)" % pos
            else:
                pos = "%d..%d" % (fs + 1, fs + 3)
            # /transl_except ALONE. Pairing it with /exception="RNA editing" is
            # rejected (TranslExceptAndRnaEditing): the reference records can use
            # the bare exception because they carry a /translation, but a feature
            # table has none, so the edited codon must be named explicitly.
            lines.append("\t\t\ttransl_except\t(pos:%s,aa:Met)" % pos)
        if ann.gene_type == "CDS":
            lines.append("\t\t\ttransl_table\t11")
            if partial or ann.protein:
                lines.append("\t\t\tcodon_start\t1")
        if trans:
            lines.append("\t\t\ttrans_splicing")
        if ann.flag == "NEEDS_REVIEW":
            lines.append("\t\t\tnote\tlow-confidence prediction; requires review")

    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return out_path


def write_report(annotations, accession, genome_len,
                  ir_boundaries, relatives, elapsed,
                  out_path):
    """Write QC report with provenance summary."""
    with open(out_path, "w") as f:
        f.write(f"Plastanno v2 Report\n")
        f.write(f"{'='*60}\n")
        f.write(f"Accession  : {accession}\n")
        f.write(f"Length     : {genome_len:,} bp\n")
        f.write(f"Runtime    : {elapsed:.1f}s\n")
        f.write(f"Date       : {datetime.now():%Y-%m-%d}\n\n")

        # IR boundaries
        f.write(f"IR/LSC/SSC Boundaries:\n")
        for region, (s,e) in ir_boundaries.items():
            f.write(f"  {region}: {s:,} - {e:,} ({e-s:,} bp)\n")

        # Warnings for atypical structures (manual curation advised)
        warns = []
        if not any(r in ir_boundaries for r in ("IRa", "IRb")):
            warns.append("No inverted repeat detected; genome annotated as single-copy. "
                         "If a standard quadripartite plastome was expected, inspect the IR "
                         "manually; if this is an IR-lacking or reduced plastome (e.g. a "
                         "heterotrophic/parasitic taxon), features near the former IR may need curation.")
        if genome_len < 80_000:
            warns.append(f"Genome length {genome_len:,} bp is unusually short (< 80 kb); "
                         f"possibly a reduced plastome — manual curation advised.")
        if warns:
            f.write("\nWarnings:\n")
            for w in warns:
                f.write(f"  ! {w}\n")

        # Top relatives
        f.write(f"\nTop relatives:\n")
        for r in relatives[:5]:
            f.write(f"  {r['accession']:<15} "
                    f"{r['pident']:.1f}% "
                    f"{r.get('genus','')}\n")

        # Gene summary
        cds   = [a for a in annotations if a.gene_type=="CDS"]
        rrna  = [a for a in annotations if a.gene_type=="rRNA"]
        trna  = [a for a in annotations if a.gene_type=="tRNA"]
        pseudo= [a for a in annotations if a.is_pseudogene]

        f.write(f"\nGene Summary:\n")
        f.write(f"  CDS       : {len(cds)}\n")
        f.write(f"  rRNA      : {len(rrna)}\n")
        f.write(f"  tRNA      : {len(trna)}\n")
        f.write(f"  Pseudogene: {len(pseudo)}\n")
        f.write(f"  Total     : {len(annotations)}\n")

        # Categorised gene table (functional groups, counting annotated instances
        # so IR-duplicated genes appear twice — e.g. rRNA totals 8).
        import textwrap, re
        from collections import Counter
        cds_inst = Counter(a.gene_name for a in cds)
        natkey = lambda g: [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", g)]

        def cat_row(label, genes, n):
            f.write(f"  {label:<34} {n:>3}\n")
            if genes:
                f.write(textwrap.fill(", ".join(genes), width=66,
                        initial_indent="      ", subsequent_indent="      ") + "\n")

        f.write("\nGene Categories (functional groups, by annotated instances):\n")
        f.write(" rRNA\n");  cat_row("rRNA", sorted({a.gene_name for a in rrna}, key=natkey), len(rrna))
        f.write(" tRNA\n");  cat_row("tRNA", sorted({a.gene_name for a in trna}, key=natkey), len(trna))
        classified = set()
        for category, subs in GENE_CATEGORIES:
            wrote = False
            for sublabel, members in subs:
                present = sorted((g for g in cds_inst if g in members), key=natkey)
                if not present:
                    continue
                if not wrote:
                    f.write(f" {category}\n"); wrote = True
                classified.update(present)
                cat_row(sublabel, present, sum(cds_inst[g] for g in present))
        other = sorted((g for g in cds_inst if g not in classified), key=natkey)
        if other:
            f.write(" Other genes\n")
            for g in other:
                cat_row(OTHER_GENES.get(g, "Other gene"), [g], cds_inst[g])
        f.write(f"  {'Total genes (instances)':<34} "
                f"{len(rrna)+len(trna)+sum(cds_inst.values()):>3}\n")

        # Confidence distribution
        high   = sum(1 for a in annotations if a.flag=="HIGH")
        medium = sum(1 for a in annotations if a.flag=="MEDIUM")
        review = sum(1 for a in annotations if a.flag=="NEEDS_REVIEW")
        f.write(f"\nConfidence:\n")
        f.write(f"  HIGH        : {high}\n")
        f.write(f"  MEDIUM      : {medium}\n")
        f.write(f"  NEEDS_REVIEW: {review}\n")

        # Every CDS NCBI's validator would reject (finalize.submission_check), in one place,
        # so a submitter knows exactly what to fix before running table2asn.
        from ..core.finalize import SUBMIT
        note_of = lambda a: next(str(n) for n in a.notes if str(n).startswith(SUBMIT))
        failing = [a for a in cds if any(str(n).startswith(SUBMIT) for n in a.notes)]
        ncbi = [a for a in failing if "would fail NCBI validation" in note_of(a)]
        repaired = [a for a in failing if a not in ncbi]
        f.write("\nSubmission check (NCBI table2asn CDS rules, genetic code 11):\n")
        if not ncbi:
            f.write("  No CDS would fail NCBI validation.\n")
        else:
            f.write(f"  {len(ncbi)} CDS would fail NCBI validation; "
                    "fix or mark them before submitting:\n")
            for a in sorted(ncbi, key=lambda x: x.start):
                why = note_of(a)[len(SUBMIT):].replace("would fail NCBI validation: ", "")
                f.write(f"  {a.gene_name:<10} {a.start + 1}..{a.end} "
                        f"({'+' if a.strand == 1 else '-'})  {why}\n")
        moved5 = [a for a in repaired if "start codon was moved" in note_of(a)]
        done3 = [a for a in repaired if a not in moved5]
        for group, what in ((done3, "completed their 3' end"),
                            (moved5, "moved their start codon")):
            if not group:
                continue
            f.write(f"  {len(group)} CDS pass only because the pipeline {what}; "
                    "check their start and exon structure:\n")
            for a in sorted(group, key=lambda x: x.start):
                why = note_of(a)[len(SUBMIT):]
                f.write(f"  {a.gene_name:<10} {a.start + 1}..{a.end} "
                        f"({'+' if a.strand == 1 else '-'})  {why}\n")

        # Genes needing review
        if review > 0:
            f.write(f"\nGenes needing review:\n")
            for a in annotations:
                if a.flag == "NEEDS_REVIEW":
                    f.write(f"  {a.gene_name:<15} "
                            f"engine={a.engine} "
                            f"C={a.confidence:.2f}\n")
                    for note in a.notes:
                        f.write(f"    → {note}\n")


def write_all(annotations, genome_seq, accession,
               genome_len, ir_boundaries, relatives,
               out_dir, prefix, no_plot=False, elapsed=0, organism=None):
    """Write all output files."""
    out_dir = Path(out_dir)
    files   = []

    # A missing taxon name is not fatal, but it must not be silent: the generic
    # placeholder otherwise propagates into /organism, the DEFINITION line and any
    # downstream figure, and would mis-assign the taxon in a GenBank submission.
    if not (organism or "").strip():
        print(f"      NOTE: no --organism given; /organism is the placeholder "
              f"'{PLACEHOLDER_ORGANISM}'.\n"
              f"            Set --organism \"Genus species\" before submitting to GenBank.")

    # One clean identifier for every output (LOCUS name, GFF3 seqid, FASTA
    # headers, filenames). `prefix` is the user's --prefix or the input filename
    # stem, which is cleaner than record.id (the raw FASTA header, which may carry
    # a trailing ', complete genome'). Sanitising here keeps all files consistent.
    seq_id = _sanitize_id(prefix)

    # Annotation content is finalised before this point (core.finalize.finalize_qc,
    # called by the pipeline). The writers only serialise; they must not change a
    # feature, or the same annotation would differ between output formats and the
    # summary counts printed earlier would no longer describe what was written.

    # Fail-closed: a tRNA whose /gene and /product name different amino acids is
    # not a tidiness problem, it is an unidentifiable locus. Every rename now
    # goes through core.trna_identity, so reaching this point means a code path
    # was missed -- surface it rather than emit the contradiction.
    from ..core import trna_identity as _TI
    # normalise first, then check. The two are deliberately separate.
    _TI.normalize_display_fields(annotations)
    _bad = []
    for ann in annotations:
        if ann.gene_type != "tRNA":
            continue
        why = _TI.validate(ann.gene_name, ann.product or "",
                           role=getattr(getattr(ann, "trna_identity", None),
                                        "role", None))
        if why:
            _bad.append("%s @%d: %s" % (ann.gene_name, ann.start, why))
    if _bad:
        raise ValueError(
            "refusing to write %d self-contradictory tRNA annotation(s); every "
            "rename must go through core.trna_identity.apply_trna_identity:\n  %s"
            % (len(_bad), "\n  ".join(_bad[:5])))

    # GenBank
    gb_path = out_dir / f"{seq_id}.gb"
    write_genbank(
        annotations, genome_seq, seq_id,
        genome_len, ir_boundaries, relatives, gb_path, organism=organism
    )
    files.append(gb_path)

    # GFF3
    gff_path = out_dir / f"{seq_id}.gff3"
    write_gff3(annotations, genome_len, seq_id, gff_path)
    files.append(gff_path)

    # NCBI feature table — the format table2asn needs for a submission
    tbl_path = out_dir / f"{seq_id}.tbl"
    write_tbl(annotations, seq_id, tbl_path, organism=organism,
              genome_len=genome_len)
    files.append(tbl_path)

    # FASTA files
    faa, ffn, frn = write_fasta_files(
        annotations, genome_seq, seq_id, out_dir, genome_len=genome_len
    )
    files.extend([faa, ffn, frn])

    # tRNA-source alternatives — the models that tRNAscan-SE, ARAGORN, the exon
    # BLAST and the tRNA BLAST proposed for a locus where a different one of them
    # was kept. This is the tRNA source-merge ledger ONLY. The reconciliation layer,
    # the fragment collapse in _select and the special-case handlers all discard
    # candidates too, and none of those decisions is recorded here yet; do not read
    # this file as a complete record of what the pipeline threw away.
    alt_path = out_dir / f"{seq_id}.trna_alternatives.tsv"
    with open(alt_path, "w") as fh:
        fh.write("# tRNA source-merge alternatives only; see writers.write_all\n")
        fh.write("\t".join(ALT_COLUMNS) + "\n")
        for ann in sorted(annotations, key=lambda x: x.start):
            for a in getattr(ann, "alternatives", []) or []:
                fh.write("\t".join(str(x) for x in [
                    ann.gene_name, ann.gene_type, ann.strand, ann.start, ann.end,
                    a.get("caller"), a.get("gene_name"), a.get("strand"),
                    a.get("start"), a.get("end"),
                    ";".join("%d-%d" % (s, e) for s, e in a.get("exons", [])),
                    a.get("wrapped"), a.get("score"), a.get("distance_bp"),
                ]) + "\n")
    files.append(alt_path)

    # Report
    rep_path = out_dir / f"{seq_id}.report"
    write_report(
        annotations, seq_id, genome_len,
        ir_boundaries, relatives, elapsed, rep_path
    )
    files.append(rep_path)

    # Circular plastome map (7th output). Skipped under no_plot (e.g. benchmarks,
    # which is why F1 numbers are unaffected). Reads back the GenBank just written
    # so the map reflects the final annotation. A plotting failure must never
    # break the annotation outputs, so it is fully guarded.
    if not no_plot:
        try:
            from Bio import SeqIO as _SeqIO
            # The in-package module only. A fallback used to load
            # scripts/viz/plastome_circular_map.py on ANY error here -- a stale copy
            # that lacked the title-escaping fix -- so a broken package module would
            # have been replaced by old code without a word. A failure now skips
            # the map with the message below.
            from plastanno.viz import plastome_circular_map as _pcm
            _rec = next(_SeqIO.parse(str(gb_path), "genbank"))
            _ir = _pcm.ir_from_blast(str(_rec.seq).upper()) or _pcm.ir_from_annotation(_rec)
            # Label the map with the real taxon; the generic placeholder is no more
            # informative than the sequence id, so prefer the id in that case.
            _org = (_rec.annotations.get("organism") or "").strip()
            if not _org or _org == PLACEHOLDER_ORGANISM:
                _org = seq_id
            _pcm.draw_map(_rec, _ir, _org, 180, str(out_dir / f"{seq_id}_map"), dpi=300)
            for _ext in (".png", ".pdf", ".svg"):
                _p = out_dir / f"{seq_id}_map{_ext}"
                if _p.exists() and _p not in files:
                    files.append(_p)
        except Exception as _e:
            print(f"      (circular map skipped: {_e})")

    return files
