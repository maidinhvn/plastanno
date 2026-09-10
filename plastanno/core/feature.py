"""
Feature data structure with provenance tracking.
Every annotation carries source engine + confidence score.
"""
from dataclasses import dataclass, field
from typing import List, Tuple, Optional

@dataclass
class Exon:
    start  : int
    end    : int
    strand : int  # 1 or -1

@dataclass
class Feature:
    """Single annotated gene with full provenance."""

    # Identity
    gene_name  : str
    gene_type  : str          # CDS / tRNA / rRNA
    product    : str = ""

    # Coordinates
    start      : int = 0
    end        : int = 0
    strand     : int = 1
    exons      : List[Tuple[int,int]] = field(default_factory=list)
    # Per-exon strand, in transcript order, parallel to `exons`. Only set for
    # trans-spliced genes whose exons lie on different strands (e.g. rps12). When
    # empty, every exon is assumed to share `strand` (the normal cis-spliced case).
    exon_strands : List[int] = field(default_factory=list)
    has_intron : bool = False
    # True for a gene whose exons are transcribed separately and spliced in trans
    # (rps12). Set explicitly rather than inferred from `exon_strands`, because the
    # IRb copy of rps12 happens to have all three exons on the same strand and
    # would otherwise look cis-spliced to the writers.
    is_trans_spliced : bool = False

    # Sequences
    protein    : str = ""     # AA sequence (CDS only)
    sequence   : str = ""     # nucleotide

    # Provenance (NEW in v2)
    engine     : str = ""     # "A", "B", or "AB"
    confidence : float = 0.0  # C score 0-1

    # Score components
    s_overlap  : float = 0.0  # Jaccard overlap A∩B/A∪B
    s_ref      : float = 0.0  # Exonerate pident/100
    s_model    : float = 0.0  # normalized HMM bitscore
    s_orf      : float = 0.0  # ORF validity

    # Flags
    is_pseudogene   : bool = False
    pseudogene_reason: str = ""
    # True when the spliced CDS length is not a multiple of 3, i.e. the boundary
    # could not be resolved into whole codons. The feature is kept exactly as
    # called; it is written as a partial CDS with no /translation.
    orf_incomplete  : bool = False
    # Number of frameshift indels Exonerate had to model in the alignment that
    # produced this feature. Non-zero means the aligned extent is not a whole
    # number of codons, so the reading frame is broken.
    frameshifts     : int = 0
    # True when the CDS begins at a recognised NON-ATG initiator (psbL/ndhD/...),
    # which in plastids is created by C-to-U RNA editing. GenBank marks such a CDS
    # with /exception="RNA editing"; without it the start codon is simply illegal.
    rna_edited_start: bool = False
    flag            : str = "HIGH"  # HIGH/MEDIUM/NEEDS_REVIEW
    notes           : List[str] = field(default_factory=list)

    def compute_confidence(self, available=None, base_weights=None):
        """
        Confidence = weighted average over the AVAILABLE signal components,
        renormalised so the weights sum to 1. `available` is the set of signal
        names meaningful for this feature (a reference-only CDS has {'ref','orf'}
        but not 'overlap'/'model'). Renormalising prevents a single-engine
        feature from being unfairly capped (the old fixed-weight formula limited
        such features to 0.4 and then forced NEEDS_REVIEW).
        """
        if not isinstance(base_weights, dict):
            base_weights = None
        if not isinstance(available, (set, list, tuple, frozenset)):
            available = None  # tolerate legacy positional calls; infer instead
        bw = base_weights or {"overlap": 0.4, "ref": 0.2, "model": 0.2, "orf": 0.2}
        vals = {"overlap": self.s_overlap, "ref": self.s_ref,
                "model": self.s_model, "orf": self.s_orf}
        if available is None:
            available = {k for k, v in vals.items() if v > 0}
            if self.gene_type == "CDS":
                available.add("orf")
        available = {k for k in available if k in bw}
        total_w = sum(bw[k] for k in available)
        self.confidence = (sum(bw[k] * vals[k] for k in available) / total_w) if total_w > 0 else 0.0
        if self.confidence >= 0.8:
            self.flag = "HIGH"
        elif self.confidence >= 0.5:
            self.flag = "MEDIUM"
        else:
            self.flag = "NEEDS_REVIEW"
        return self.confidence

    def to_dict(self):
        return {
            "gene_name"       : self.gene_name,
            "gene_type"       : self.gene_type,
            "start"           : self.start,
            "end"             : self.end,
            "strand"          : self.strand,
            "exons"           : self.exons,
            "has_intron"      : self.has_intron,
            "is_trans_spliced": self.is_trans_spliced,
            "protein"         : self.protein,
            "engine"          : self.engine,
            "confidence"      : round(self.confidence, 3),
            "s_overlap"       : round(self.s_overlap, 3),
            "s_ref"           : round(self.s_ref, 3),
            "s_model"         : round(self.s_model, 3),
            "s_orf"           : round(self.s_orf, 3),
            "flag"            : self.flag,
            "is_pseudogene"   : self.is_pseudogene,
            "orf_incomplete"  : self.orf_incomplete,
            "frameshifts"     : self.frameshifts,
            "rna_edited_start": self.rna_edited_start,
            "pseudogene_reason": self.pseudogene_reason,
            "notes"           : self.notes,
        }
