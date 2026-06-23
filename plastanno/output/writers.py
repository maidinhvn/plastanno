"""
Output writers for Plastanno v2.
Writes: .gb, .gff3, .faa, .ffn, .frn, .report, and (unless no_plot) a circular
plastome map (_map.png/.pdf/.svg).

Improvements over v1:
- Provenance in GenBank notes
- Confidence flags in GFF3
- Clean FASTA headers
"""
from pathlib import Path
from typing import List, Dict
from datetime import datetime
from ..core.feature import Feature


def _strand_str(strand):
    if strand == 1  or strand == "+": return "+"
    if strand == -1 or strand == "-": return "-"
    return "."


def write_genbank(annotations, genome_seq, accession,
                   genome_len, ir_boundaries, relatives,
                   out_path):
    """Write GenBank format with provenance notes."""
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

    seq    = Seq(genome_seq)
    record = SeqRecord(
        seq,
        id          = accession,
        name        = accession,
        description = "Annotated by Plastanno v2",
    )
    record.annotations["molecule_type"] = "DNA"
    record.annotations["topology"]      = "circular"
    record.annotations["date"]          = \
        datetime.now().strftime("%d-%b-%Y").upper()

    # Source feature
    record.features.append(SeqFeature(
        FeatureLocation(0, genome_len),
        type = "source",
        qualifiers = {
            "organism": ["Viridiplantae"],
            "mol_type": ["genomic DNA"],
        }
    ))

    # Add annotations
    for ann in sorted(annotations, key=lambda x: x.start):
        # Build exon locations
        if ann.exon_strands and len(ann.exon_strands) == len(ann.exons):
            # Trans-spliced gene: exons may lie on different strands. Keep the
            # given transcript order and per-exon strand (do not re-sort).
            parts = [FeatureLocation(s, e, st)
                     for (s, e), st in zip(ann.exons, ann.exon_strands)]
            loc = CompoundLocation(parts)
        elif len(ann.exons) > 1:
            parts = [
                FeatureLocation(s, e, ann.strand)
                for s, e in (sorted(ann.exons)[::-1] if ann.strand == -1 else sorted(ann.exons))
            ]
            loc = CompoundLocation(parts)
        else:
            loc = FeatureLocation(
                ann.start, ann.end, ann.strand
            )

        qualifiers = {
            "gene"   : [ann.gene_name],
            "product": [ann.product or ann.gene_name],
        }

        # Add protein sequence for CDS
        if ann.gene_type == "CDS" and ann.protein:
            qualifiers["translation"] = [ann.protein]
            qualifiers["codon_start"] = ["1"]

        # Provenance notes (v2 improvement)
        note = (f"engine={ann.engine}; "
                f"confidence={ann.confidence:.2f}; "
                f"flag={ann.flag}")
        if ann.notes:
            note += "; " + "; ".join(ann.notes)
        qualifiers["note"] = [note]

        if ann.is_pseudogene:
            qualifiers["pseudo"] = [""]
            qualifiers["note"][0] += \
                f"; pseudogene_reason={ann.pseudogene_reason}"

        # Gene feature
        record.features.append(SeqFeature(
            loc, type="gene", qualifiers={"gene":[ann.gene_name]}
        ))

        # CDS/tRNA/rRNA feature
        record.features.append(SeqFeature(
            loc, type=ann.gene_type, qualifiers=qualifiers
        ))

    SeqIO.write(record, out_path, "genbank")
    return record


def write_gff3(annotations, genome_len, accession, out_path):
    """Write GFF3 format with confidence flags."""
    with open(out_path, "w") as f:
        f.write("##gff-version 3\n")
        f.write(f"##sequence-region {accession} 1 {genome_len}\n")

        for ann in sorted(annotations, key=lambda x: x.start):
            strand = _strand_str(ann.strand)
            attrs  = (f"ID={ann.gene_name}_{ann.start};"
                      f"Name={ann.gene_name};"
                      f"engine={ann.engine};"
                      f"confidence={ann.confidence:.2f};"
                      f"flag={ann.flag}")
            if ann.is_pseudogene:
                attrs += ";pseudo=true"

            # Gene line
            f.write("\t".join([
                accession, "Plastanno", "gene",
                str(ann.start+1), str(ann.end),
                ".", strand, ".", attrs
            ]) + "\n")

            # Feature line
            feat_type = ann.gene_type
            f.write("\t".join([
                accession, "Plastanno", feat_type,
                str(ann.start+1), str(ann.end),
                f"{ann.confidence:.2f}", strand, ".",
                f"Parent={ann.gene_name}_{ann.start};"
                f"product={ann.product or ann.gene_name}"
            ]) + "\n")


def write_fasta_files(annotations, genome_seq,
                       accession, out_dir):
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
                if ann.protein:
                    faa.write(f"{header}\n{ann.protein}\n")
                # Extract nucleotide
                seq = _extract_seq(
                    genome_seq, ann.exons, ann.strand, ann.exon_strands
                )
                ffn.write(f"{header}\n{seq}\n")

            elif ann.gene_type in ("tRNA","rRNA"):
                seq = _extract_seq(
                    genome_seq, ann.exons or
                    [(ann.start, ann.end)], ann.strand
                )
                frn.write(f"{header}\n{seq}\n")

    return faa_path, ffn_path, frn_path


def _populate_proteins(annotations, genome_seq):
    """Translate each CDS's spliced coding sequence into Feature.protein when the
    engines have not set it (the usual case). Uses the same coding-sequence
    extraction as the .ffn writer (so trans-spliced rps12 is handled), the plastid
    translation table (11), strips the terminal stop, and renders a recognised
    alternative initiator codon (GTG/ACG/TTG for psbL/ndhD/rps19/rpl2/ycf1/rps12)
    as Met — matching the GenBank /translation convention. Does not alter any
    coordinates, so annotation scoring is unchanged."""
    from Bio.Seq import Seq
    from ..core.reconcile import SPECIAL_START_CODONS
    for ann in annotations:
        if ann.gene_type != "CDS" or ann.protein:
            continue
        exons = ann.exons or [(ann.start, ann.end)]
        seq = _extract_seq(genome_seq, exons, ann.strand, ann.exon_strands)
        seq = seq[: len(seq) // 3 * 3]
        if len(seq) < 3:
            continue
        aa = str(Seq(seq).translate(table=11))
        if aa.endswith("*"):
            aa = aa[:-1]
        starts = ("ATG",) + tuple(SPECIAL_START_CODONS.get(ann.gene_name, ()))
        if aa and seq[:3].upper() in starts and aa[0] != "M":
            aa = "M" + aa[1:]
        ann.protein = aa


def _extract_seq(genome_seq, exons, strand, exon_strands=None):
    """Extract and concatenate exon sequences (coding orientation)."""
    from Bio.Seq import Seq
    if exon_strands and len(exon_strands) == len(exons):
        # Trans-spliced: take each exon in its own coding orientation, in the
        # given transcript order.
        return "".join(
            str(Seq(genome_seq[s:e]).reverse_complement()) if st == -1
            else genome_seq[s:e]
            for (s, e), st in zip(exons, exon_strands)
        )
    seq = "".join(
        genome_seq[s:e] for s, e in exons
    )
    if strand == -1:
        seq = str(Seq(seq).reverse_complement())
    return seq


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
               out_dir, prefix, no_plot=False, elapsed=0):
    """Write all output files."""
    out_dir = Path(out_dir)
    files   = []

    # Fill in CDS protein translations (used by both the GenBank /translation
    # qualifier and the .faa file); the engines leave Feature.protein empty.
    _populate_proteins(annotations, genome_seq)

    # QC: a spliced CDS whose translation carries an internal stop codon is
    # biologically invalid — a mis-placed exon boundary, a pseudogene, or an
    # unrecognised RNA-editing site. Flag it NEEDS_REVIEW with a note so users
    # inspect it rather than trusting it silently. (Coordinates are unchanged.)
    for ann in annotations:
        if ann.gene_type == "CDS" and ann.protein and "*" in ann.protein:
            n_stop = ann.protein.count("*")
            ann.flag = "NEEDS_REVIEW"
            ann.notes.append(
                f"internal stop codon(s) (n={n_stop}): possible pseudogene or RNA-editing site")

    # GenBank
    gb_path = out_dir / f"{prefix}.gb"
    write_genbank(
        annotations, genome_seq, accession,
        genome_len, ir_boundaries, relatives, gb_path
    )
    files.append(gb_path)

    # GFF3
    gff_path = out_dir / f"{prefix}.gff3"
    write_gff3(annotations, genome_len, accession, gff_path)
    files.append(gff_path)

    # FASTA files
    faa, ffn, frn = write_fasta_files(
        annotations, genome_seq, accession, out_dir
    )
    files.extend([faa, ffn, frn])

    # Report
    rep_path = out_dir / f"{prefix}.report"
    write_report(
        annotations, accession, genome_len,
        ir_boundaries, relatives, elapsed, rep_path
    )
    files.append(rep_path)

    # Circular plastome map (7th output). Skipped under no_plot (e.g. benchmarks,
    # which is why F1 numbers are unaffected). Reads back the GenBank just written
    # so the map reflects the final annotation. A plotting failure must never
    # break the annotation outputs, so it is fully guarded.
    if not no_plot:
        try:
            import importlib.util
            from Bio import SeqIO as _SeqIO
            _mp = Path(__file__).resolve().parents[2] / "scripts" / "viz" / "plastome_circular_map.py"
            if _mp.exists():
                _spec = importlib.util.spec_from_file_location("_pcm", _mp)
                _pcm = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(_pcm)
                _rec = next(_SeqIO.parse(str(gb_path), "genbank"))
                _ir = _pcm.ir_from_blast(str(_rec.seq).upper()) or _pcm.ir_from_annotation(_rec)
                _org = (_rec.annotations.get("organism") or "").strip() or accession
                _pcm.draw_map(_rec, _ir, _org, 180, str(out_dir / f"{prefix}_map"), dpi=300)
                for _ext in (".png", ".pdf", ".svg"):
                    _p = out_dir / f"{prefix}_map{_ext}"
                    if _p.exists() and _p not in files:
                        files.append(_p)
        except Exception as _e:
            print(f"      (circular map skipped: {_e})")

    return files
