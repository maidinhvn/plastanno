"""One source of truth for tRNA identity.

THE DEFECT THIS EXISTS TO PREVENT
    Two places renamed a tRNA by assigning `gene_name` and leaving `product`
    alone. On 30 development genomes that produced 121 features (10.8% of all
    predicted tRNA) asserting two different amino acids at once --
    `gene=trnG-UCC` beside `product=tRNA-Ser`. The reference rate is 0.5%. Such a
    record is not merely untidy: `gene_naming.classify_trna` reports
    "family disagreement" and the locus becomes unidentifiable to any consumer.

FOUR THINGS, KEPT APART
    raw detector identity   what ARAGORN or the exon DB called it
    rename evidence         what BLAST, the anticodon, or IR position says
    canonical identity      the decision
    resolution status       how firm that decision is

`gene_name` and `product` are BOTH derived from the canonical identity, by
construction, so they cannot disagree. The detector's original call is never
deleted -- it is recorded as provenance, because an ARAGORN/BLAST disagreement is
a real observation about the locus, not stray text.

WHAT THIS DOES NOT DO
    Making the two qualifiers agree does not establish that the chosen name is
    biologically correct. Precedence here is a documented policy, not evidence.
    Adjudicating the identity is the tRNA audit, separately.
"""
from dataclasses import dataclass, field, replace
import re

RESOLVED = "RESOLVED"
FAMILY_ONLY = "FAMILY_ONLY"
UNRESOLVED = "UNRESOLVED"
CONFLICT = "CONFLICT"
#: a display identity chosen by BLAST over a disagreeing detector call. It is
#: PROVISIONAL: no biological audit has shown BLAST is right here, only that the
#: rename exists because ARAGORN reads intron-tRNA anticodons unreliably. It must
#: not be treated as a confirmed identity and must not raise any support score.
BLAST_SELECTED = "BLAST_SELECTED"

#: evidence sources. BLAST against the tRNA DB SELECTS THE DISPLAY IDENTITY when
#: it disagrees with the structural caller -- that is the historical behaviour
#: this refactor preserves, and the documented reason the rename exists. It is
#: NOT called authoritative: no biological audit has established that BLAST is
#: correct here. It never touches coordinates.
BLAST_TRNA_DB = "BLAST_TRNA_DB"
ARAGORN_STRUCTURE = "ARAGORN_STRUCTURE"
POSITIONAL_IR = "POSITIONAL_IR"
_SELECTS_DISPLAY_IDENTITY = {BLAST_TRNA_DB}

AA1_TO_3 = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
            "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
            "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
            "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
            "U": "Sec"}
_AA3_TO_1 = {v.lower(): k for k, v in AA1_TO_3.items()}
_AA3_TO_1["fmet"] = "M"
_AA3_TO_1["imet"] = "M"

INITIATOR, ELONGATOR = "INITIATOR", "ELONGATOR"


@dataclass(frozen=True)
class TRNAIdentity:
    family: str = None          #: one-letter amino acid
    anticodon: str = None       #: RNA triplet, e.g. "UCC"
    role: str = None            #: INITIATOR / ELONGATOR -- methionine only
    status: str = UNRESOLVED

    def same_family(self, other):
        return (self.family is not None and other.family is not None
                and self.family == other.family)

    def disagrees_with(self, other):
        """True only where BOTH sides resolve an axis and the values differ."""
        for a in ("family", "anticodon", "role"):
            x, y = getattr(self, a), getattr(other, a)
            if x is not None and y is not None and x != y:
                return True
        return False


def from_gene_name(name):
    """Parse a plastid tRNA gene name. Never guesses: an unparseable name gives
    UNRESOLVED rather than a fabricated family."""
    n = (name or "").strip()
    m = re.match(r"^trn(fM|f?[A-Za-z])(?:-([ACGUTacgut]{3}))?$", n)
    if not m:
        return TRNAIdentity(status=UNRESOLVED)
    tag, tri = m.group(1), m.group(2)
    role = None
    if tag.lower() in ("fm",):
        fam, role = "M", INITIATOR
    else:
        fam = tag.upper()[-1]
    if fam not in AA1_TO_3:
        return TRNAIdentity(status=UNRESOLVED)
    anti = tri.upper().replace("T", "U") if tri else None
    if fam == "M" and role is None and anti:
        role = ELONGATOR          # plain trnM is the elongator by convention
    return TRNAIdentity(family=fam, anticodon=anti, role=role,
                        status=RESOLVED if anti else FAMILY_ONLY)


def gene_name(idn):
    """The /gene value. Derived from the identity, never set independently."""
    if idn.family is None:
        return "trnX"
    tag = "fM" if (idn.family == "M" and idn.role == INITIATOR) else idn.family
    return "trn%s-%s" % (tag, idn.anticodon) if idn.anticodon else "trn%s" % tag


def product(idn):
    """The /product value, from the SAME identity as gene_name.

    INSDC wants the RNA product named by its amino acid, not the gene symbol.
    /product="trnI-GAU" is a gene symbol and semantically wrong; the value must
    be "tRNA-Ile". An unresolved family gives tRNA-OTHER rather than an invented
    amino acid -- and never tRNA-Stop or tRNA-Pyl, which are not plastid tRNA
    products and reached the output from a raw ARAGORN anticodon call.

    trnfM carries product tRNA-Met: INSDC has no controlled value "fMet", and the
    initiator role is carried by the trnfM symbol and by `role`, not by /product.
    """
    if idn.family is None:
        return "tRNA-OTHER"
    return "tRNA-%s" % AA1_TO_3[idn.family]


def anticodon_triplet(idn):
    """The anticodon as the identity holds it, or None.

    NOT emitted as an /anticodon qualifier: INSDC requires
    (pos:<range>,aa:<AA>,seq:<nnn>) and this pipeline does not compute the
    anticodon POSITION, so emitting the qualifier would be invalid. The triplet
    is produced here so the validator can check the third leg of
    gene family <-> product amino acid <-> anticodon, which the gene symbol
    already encodes. Adding a real /anticodon qualifier needs position
    calculation and would be a new output field (annotation_schema 2).
    """
    return idn.anticodon


def display_fields(idn):
    """THE single factory. All three display values from one identity.

    Every caller must take all three from here. Setting one of them
    independently is what produced 233 inconsistent tRNA records across 30
    genomes: 117 where /gene and /product named different amino acids, 112 where
    /product was the gene symbol, and 4 where /product was tRNA-Stop or
    tRNA-Pyl.
    """
    return gene_name(idn), product(idn), anticodon_triplet(idn)


def resolve(detector, evidence, source, detector_symbol=None):
    """(canonical identity, provenance note).

    Agreement merges detail. Disagreement is resolved ONLY by a source that is
    authoritative on identity, and the superseded call is recorded. A
    non-authoritative source that disagrees yields CONFLICT: the detector's call
    is kept and the locus is marked, rather than forcing the qualifiers to agree
    and manufacturing false confidence.
    """
    if evidence is None or evidence.family is None:
        return detector, ""
    if not detector.disagrees_with(evidence):
        # No contradiction -- but the detector may have resolved NOTHING, as
        # ARAGORN does when it reports trn?-UUA. That raw call must still be
        # recorded: it is why the product read tRNA-Stop, and losing it would
        # erase the only trace of what the structural caller actually said.
        raw_note = ""
        if detector.family is None and detector_symbol:
            raw_note = ("raw_detector_symbol=%s raw_detector_family=NA "
                        "identity_selection=%s coordinates_unchanged=true"
                        % (detector_symbol, source))
        merged = TRNAIdentity(
            family=evidence.family or detector.family,
            anticodon=evidence.anticodon or detector.anticodon,
            role=evidence.role or detector.role,
            status=UNRESOLVED)
        st = RESOLVED if merged.anticodon else (
            FAMILY_ONLY if merged.family else UNRESOLVED)
        if merged.family == "M" and merged.role is None:
            st = FAMILY_ONLY      # the two methionines are different genes
        # a merge over an unresolved detector is SELECTED, not independently
        # resolved: only the evidence had an opinion
        if raw_note:
            st = BLAST_SELECTED if source in _SELECTS_DISPLAY_IDENTITY else st
        return replace(merged, status=st), raw_note
    if source in _SELECTS_DISPLAY_IDENTITY:
        # The superseded call is recorded in a form that CANNOT be re-parsed as a
        # gene name. Writing "trnS-CGA" here put a second family into /note, and
        # gene_naming.classify_trna -- which reads /note -- then reported a family
        # disagreement on a record whose qualifiers were already consistent. The
        # provenance must not corrupt the identity it documents.
        # NAMESPACED so provenance cannot be mistaken for an official qualifier
        # and cannot be re-read as identity evidence by any classifier.
        note = ("raw_detector_family=%s raw_detector_anticodon=%s "
                "identity_selection=%s coordinates_unchanged=true"
                % (detector.family or "NA", detector.anticodon or "NA", source))
        # PROVISIONAL, not RESOLVED: the detector and BLAST genuinely disagree
        # and only one was displayed. Both are kept in provenance.
        return replace(evidence, status=BLAST_SELECTED), note
    note = ("raw_detector_family=%s raw_detector_anticodon=%s "
            "rename_evidence_family=%s rename_evidence_anticodon=%s "
            "identity_selection=NONE identity_status=CONFLICT source=%s"
            % (detector.family or "NA", detector.anticodon or "NA",
               evidence.family or "NA", evidence.anticodon or "NA", source))
    return replace(detector, status=CONFLICT), note


def apply_trna_identity(feature, idn, note=""):
    """Set gene_name AND product together. The only sanctioned way to rename."""
    feature.gene_name = gene_name(idn)
    feature.product = product(idn)
    feature.trna_identity = idn
    if note:
        feature.notes.append(note)
    return feature


#: the only product values this pipeline may emit for a tRNA
_VALID_PRODUCTS = frozenset(["tRNA-%s" % v for v in AA1_TO_3.values()] +
                            ["tRNA-OTHER"])


def validate(gene, prod, anticodon=None, role=None):
    """Report inconsistency. NEVER repairs -- that is the factory's job.

    Checks the whole tuple, not just a pair:
        gene family <-> product amino acid <-> anticodon amino acid <-> CAU role

    Returns a reason string, or None when consistent. A malformed product is
    reported explicitly rather than declining to judge: an earlier version only
    compared well-formed products and silently passed /product="trnI-GAU" and
    /product="tRNA-Stop", so 116 invalid records went unflagged.
    """
    g = from_gene_name(gene)
    p = (prod or "").strip()

    if not p:
        return "no product; a tRNA must name its RNA product"
    if p not in _VALID_PRODUCTS:
        if re.match(r"^trn[A-Za-z]{1,3}(-[ACGUTacgut]{3})?$", p):
            return ("product %r is a gene symbol, not an RNA product name "
                    "(expected tRNA-<AminoAcid>)" % p)
        return ("product %r is not a valid tRNA product; permitted values are "
                "tRNA-<one of the 20 amino acids or Sec> or tRNA-OTHER" % p)

    if g.family is None:
        # an unidentifiable gene may only carry the non-committal product
        return (None if p == "tRNA-OTHER" else
                "gene %r resolves no amino acid, so product must be tRNA-OTHER, "
                "not %r" % (gene, p))
    want = "tRNA-%s" % AA1_TO_3[g.family]
    if p != want:
        return ("gene %s is %s but product is %s"
                % (gene, AA1_TO_3[g.family], p))

    # third leg: the anticodon, where one is supplied
    if anticodon and g.anticodon and anticodon.upper().replace("T", "U") != g.anticodon:
        return ("gene %s encodes anticodon %s but %s was supplied"
                % (gene, g.anticodon, anticodon))

    # fourth leg: the CAU role. trnI-CAU has no functional role; the two
    # methionines are distinguished by it and must not be merged.
    if g.family == "M":
        want_role = INITIATOR if gene.startswith("trnfM") else ELONGATOR
        if role is not None and role != want_role:
            return ("gene %s implies role %s but %s was supplied"
                    % (gene, want_role, role))
    elif role not in (None, "NA"):
        return "gene %s is not methionine, so it has no functional role" % gene
    return None


def normalize_display_fields(annotations):
    """Derive /gene and /product for every tRNA from ONE identity.

    THE NORMALIZER. Separate from `validate`, which only reports: a checker that
    repairs cannot also be trusted to detect. Every tRNA leaves here with display
    fields produced by `display_fields`, so the three classes of inconsistency
    measured on 30 development genomes -- 117 gene/product disagreements, 112
    products holding a gene symbol, 4 products reading tRNA-Stop or tRNA-Pyl --
    cannot be emitted.

    The identity is taken from `feature.trna_identity` when a resolver already
    set one, otherwise parsed from the gene symbol. Coordinates, exons, strand,
    scores and flags are untouched.

    Returns the number of features whose display fields changed.
    """
    n = 0
    for f in annotations:
        if getattr(f, "gene_type", None) != "tRNA":
            continue
        idn = getattr(f, "trna_identity", None) or from_gene_name(f.gene_name)
        g, pr, _anti = display_fields(idn)
        if (f.gene_name, f.product) != (g, pr):
            f.gene_name, f.product = g, pr
            n += 1
    return n
