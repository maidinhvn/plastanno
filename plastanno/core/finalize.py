"""Finalise annotation content before anything is written or counted.

Every check that can change a feature lives here, and it runs once, in the pipeline,
before the summary counts are printed and before any output file is produced. It
used to live in the writer, which meant the annotation was still being modified
while the first format was being serialised: the counts printed to the user
described a different state from the files on disk, and a check added to one writer
would not apply to the others.

Writers may not modify a feature after this point.
"""

from . import coords as _coords

# Every note this module adds carries this prefix, so a second pass can remove the
# notes the previous one left before writing its own. Without it, finalising twice
# appends the same sentence twice, and a note describing a problem that has since
# been fixed stays on the feature for good.
QC = "[qc] "


def _populate_proteins(annotations, genome_seq, genome_len=None):
    """Translate each CDS's spliced coding sequence into Feature.protein.

    Uses the coordinate contract, so an origin-crossing or trans-spliced CDS is read
    in transcription order. Strips the terminal stop and renders a recognised
    alternative initiator as Met, matching the GenBank /translation convention.

    The translation is recomputed from the coordinates every time rather than kept
    when already present. A step between the engines and here may move a boundary —
    ORF completion, the ycf1 handler, the missing-first-exon rescue all do — and a
    protein carried over from the pre-move coordinates then no longer matches the
    CDS it is written next to.
    """
    from Bio.Seq import Seq
    from .reconcile import SPECIAL_START_CODONS
    for ann in annotations:
        if ann.gene_type != "CDS" or ann.is_pseudogene:
            continue
        ann.protein = ""
        seq = _coords.extract(genome_seq, ann, genome_len)
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


def finalize_qc(annotations, genome_seq, boundary_conflict_bp=None, genome_len=None):
    """Run every content-changing check, once, before output.

    Idempotent: running it twice on the same annotations gives exactly the same
    result as running it once, and running it again after a coordinate has moved
    describes the new coordinate. Everything it derives — the translation, the
    partial-CDS mark, the RNA-edited-start mark, its own notes and its effect on
    the flag — is recomputed from the final coordinates on every pass, and the
    previous pass's output is cleared first. That property is what makes it safe to
    call from more than one place, which is how the old writer-embedded version
    produced counts that disagreed with the files.

    `boundary_conflict_bp` decides which recorded alternative is worth a reviewer's
    attention. **The default is None, which flags nothing**: every alternative is
    recorded in `Feature.alternatives` and written to the run's alternatives file,
    but no rule promotes one to NEEDS_REVIEW on its own.

    Two candidate rules were considered and rejected. A fixed 3 bp threshold came
    from counting disagreements on one development genome, where 21 of 26 were one
    or two bases. And "one caller wrapped the origin, the other did not" is not a
    biological difference at all: whether a feature crosses position 1 depends on
    where the submitter placed the origin, and the same feature stops being wrapped
    if the sequence is rotated. Neither belongs in a default.
    """
    from .reconcile import SPECIAL_START_CODONS

    if genome_len is None:
        genome_len = len(genome_seq) or None

    # Clear what a previous pass derived, so nothing survives a coordinate change.
    for ann in annotations:
        if ann.pre_qc_flag is None:
            ann.pre_qc_flag = ann.flag        # first pass: remember reconciliation's verdict
        else:
            ann.flag = ann.pre_qc_flag        # later passes: start from it again
        ann.notes = [n for n in ann.notes if not str(n).startswith(QC)]
        ann.orf_incomplete = False
        ann.rna_edited_start = False
        ann.geometry_invalid = False

    _populate_proteins(annotations, genome_seq, genome_len)

    for ann in annotations:
        # A structurally impossible exon model must not pass silently. The measuring
        # helpers merge overlapping exons so a length stays meaningful, while the
        # translator reads the model as written — so such a feature reports one
        # length to QC and hands a longer sequence to the translation, and neither
        # number is wrong on its own terms. Flag it rather than pick a side.
        problems = _coords.geometry_problems(ann, genome_len)
        if problems:
            # Quarantine, not just a flag. The spliced length and the translated
            # sequence disagree for such a model, so a protein derived from it is
            # not a fact about the genome — withhold it and write the feature as a
            # partial rather than let a wrong translation reach a submission.
            ann.geometry_invalid = True
            ann.orf_incomplete = True
            ann.protein = ""
            ann.flag = "NEEDS_REVIEW"
            for m in problems:
                ann.notes.append(QC + "invalid exon model: " + m)
            ann.notes.append(
                QC + "translation withheld and the feature written as partial: the "
                "spliced length (%d bp) and the sequence the model yields (%d bp) "
                "do not agree"
                % (_coords.spliced_length(ann, genome_len),
                   len(_coords.extract(genome_seq, ann, genome_len))))

        # a locus two callers disagreed about is not settled by source priority
        if getattr(ann, "alternatives", None):
            worst = max((a.get("distance_bp", 0) for a in ann.alternatives), default=0)
            if boundary_conflict_bp is not None and worst > boundary_conflict_bp:
                ann.flag = "NEEDS_REVIEW"
                ann.notes.append(
                    QC + "boundary conflict: another caller placed this locus at %s"
                    % "; ".join("%d-%d (%d bp away%s)"
                                % (a["start"], a["end"], a["distance_bp"], "")
                                for a in ann.alternatives))

        if ann.gene_type != "CDS":
            continue

        if ann.protein and "*" in ann.protein:
            n_stop = ann.protein.count("*")
            ann.flag = "NEEDS_REVIEW"
            ann.notes.append(QC + "internal stop codon(s) (n=%d): possible pseudogene "
                             "or RNA-editing site" % n_stop)

        alts = SPECIAL_START_CODONS.get(ann.gene_name)
        if alts and not ann.is_pseudogene:
            if _coords.extract(genome_seq, ann, genome_len)[:3].upper() in tuple(alts):
                ann.rna_edited_start = True

        # a CDS that is not a whole number of codons cannot be read; keep the
        # coordinates, drop the meaningless translation, mark it partial
        if not ann.is_pseudogene and _coords.spliced_length(ann, genome_len) % 3:
            ann.orf_incomplete = True
            ann.protein = ""
            ann.flag = "NEEDS_REVIEW"
            ann.notes.append(
                QC + "spliced length %d is not a multiple of 3: unresolved boundary; "
                "annotated as a partial CDS, translation withheld"
                % _coords.spliced_length(ann, genome_len))
    return annotations


# Notes written by submission_check carry this prefix, so a second pass can clear them.
SUBMIT = "[submission] "
# The terminal-stop step (annotate/refine_splice) marks a CDS whose 3' end it completed or
# trimmed beyond what 3.0.1 did, so submission_check keeps it NEEDS_REVIEW even when the ORF
# is now valid: a clean 3' end says nothing about the rest of the gene model.
COMPLETED_3P = "[3' completed] "
# The start-codon step (annotate/refine_splice) marks a CDS whose first codon was not a start
# codon and which it moved to one nearby; submission_check keeps it NEEDS_REVIEW for the same
# reason: a valid start codon does not make it the right one.
MOVED_5P = "[start moved] "
# Engine A marks a CDS it found in the inverted repeat because its catalog region held no usable
# hit (an IR that has expanded over the region). submission_check keeps it NEEDS_REVIEW: a gene
# found outside its usual region is worth a look even when its ORF is valid.
OUT_OF_REGION = "[outside its region] "
# Initiation codons of NCBI genetic code 11, the ones table2asn accepts as a CDS start.
TABLE11_STARTS = ("ATG", "GTG", "TTG", "CTG", "ATT", "ATC", "ATA")
STOP_CODONS = ("TAA", "TAG", "TGA")


def submission_check(annotations, genome_seq, genome_len=None):
    """Flag every CDS that NCBI's validator would reject. Changes flags and notes only.

    Until this check, most CDS that NCBI's validator rejects reached the user flagged MEDIUM
    or HIGH: the pipeline knew the codons and said nothing.

    Three checks are applied, those of table2asn for genetic code 11, to each CDS as the
    writers will write it:

      no valid start  the first codon is not a table-11 initiation codon, and the CDS is not
                      written with the RNA-editing exception (rna_edited_start)
      no stop         the 3' end is not a stop codon
      internal stop   a stop codon occurs before the last codon

    A CDS written as partial (orf_incomplete: '<' and '>' on both ends) is exempt from the
    first two, exactly as it is for the validator, but not from the third. Pseudogenes are
    written without a CDS and are skipped.

    A failing CDS becomes NEEDS_REVIEW, with a note naming each problem and the flag it had.
    So does a CDS whose 3' end the terminal-stop step completed or trimmed (a COMPLETED_3P note),
    or whose start the start-codon step moved (a MOVED_5P note), even when its ORF is now valid:
    the repair made it pass the validator, and a clean ORF says nothing about whether its start
    or its splice sites are the right ones. So does a CDS Engine A found in the inverted repeat,
    outside its catalog region (an OUT_OF_REGION note).
    Nothing else changes: no coordinate, no feature, no confidence. That is why this runs
    after revoke_implausible_rescues. Inside finalize_qc its flags would feed that gate and
    remove rescued CDS, which is a change to the annotation, not a warning about it.

    Idempotent. The notes this function wrote are cleared first, and a flag it lowered is
    lowered again from the same evidence, with the same note: the flag it found is read back
    from its own earlier note, so a second call does not report "NEEDS_REVIEW" as the flag it
    lowered. It never raises a flag back up. Restoring flags is finalize_qc's job, which is
    why this must run after it (and after nothing that moves a coordinate).
    Returns [(feature, [problem, ...])].
    """
    import re
    if genome_len is None:
        genome_len = len(genome_seq) or None
    failing = []
    for ann in annotations:
        lowered_from = None
        for n in ann.notes:
            m = re.search(r"flag lowered from (\w+)", str(n)) if str(n).startswith(SUBMIT) else None
            if m:
                lowered_from = m.group(1)
        ann.notes = [n for n in ann.notes if not str(n).startswith(SUBMIT)]
        if ann.gene_type != "CDS" or ann.is_pseudogene:
            continue
        seq = _coords.extract(genome_seq, ann, genome_len).upper()
        whole = len(seq) - len(seq) % 3
        codons = [seq[i:i + 3] for i in range(0, whole, 3)]
        problems = []
        if not getattr(ann, "orf_incomplete", False):
            if codons and codons[0] not in TABLE11_STARTS and \
                    not getattr(ann, "rna_edited_start", False):
                problems.append("no valid start codon (%s)" % codons[0])
            if len(seq) % 3 or not codons or codons[-1] not in STOP_CODONS:
                problems.append("no stop codon at the 3' end (%s)" % (seq[-3:] or "-"))
        # the last codon is the terminal stop when there is one, and is examined above
        internal = sum(1 for c in codons[:-1] if c in STOP_CODONS)
        if internal:
            problems.append("%d internal stop codon(s)" % internal)
        completed = [str(n)[len(COMPLETED_3P):] for n in ann.notes if str(n).startswith(COMPLETED_3P)]
        moved = [str(n)[len(MOVED_5P):] for n in ann.notes if str(n).startswith(MOVED_5P)]
        outside = [str(n)[len(OUT_OF_REGION):] for n in ann.notes if str(n).startswith(OUT_OF_REGION)]
        if not problems and not completed and not moved and not outside:
            continue
        was = ann.flag
        if was == "NEEDS_REVIEW" and lowered_from:     # lowered by this function's earlier pass
            was = lowered_from
        ann.flag = "NEEDS_REVIEW"
        repairs = []
        if moved:
            repairs.append("its start codon was moved by the pipeline (%s)" % "; ".join(moved))
        if completed:
            repairs.append("its 3' end was completed by the pipeline (%s)" % "; ".join(completed))
        if outside:
            repairs.append("it was %s" % "; ".join(outside))
        if problems:
            text = "; ".join(["would fail NCBI validation: " + "; ".join(problems)] + repairs)
        else:
            text = "; ".join(repairs) + ": check its start and exon structure before submitting"
        ann.notes.append(SUBMIT + text
                         + ("" if was == "NEEDS_REVIEW" else "; flag lowered from " + was))
        failing.append((ann, problems
                        or ["start codon moved by the pipeline"] * bool(moved)
                        + ["3' end completed by the pipeline"] * bool(completed)
                        + ["found outside its catalog region"] * bool(outside)))
    return failing


def assign_cds_products(annotations, gene_catalog):
    """Give each CDS its protein name from the catalog, once, before output.

    Until now no CDS but rps12 reached the writers with a product, and each writer
    substituted the gene symbol: /product="psbA" where GenBank practice is
    "photosystem II protein D1", in the .gb, the .gff3 and the .tbl meant for
    submission (parked/product_qualifier/FINDINGS.md). The catalog held the
    protein names all along; nothing read them.

    A product already set is kept (rps12's, every tRNA's and rRNA's). A pseudogene
    is written without a product feature, so it is skipped. A CDS whose name the
    catalog cannot name keeps the symbol fallback and says so in its notes.

    Returns (named, unnamed).
    """
    named = unnamed = 0
    for ann in annotations:
        if ann.gene_type != "CDS" or ann.product or ann.is_pseudogene:
            continue
        product = (gene_catalog.get(ann.gene_name) or {}).get("product")
        if product and product != ann.gene_name:
            ann.product = product
            named += 1
        else:
            ann.notes.append(QC + "no protein name in the catalog; /product "
                             "falls back to the gene symbol")
            unnamed += 1
    return named, unnamed
