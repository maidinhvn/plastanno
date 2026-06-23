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
            "protein"         : self.protein,
            "engine"          : self.engine,
            "confidence"      : round(self.confidence, 3),
            "s_overlap"       : round(self.s_overlap, 3),
            "s_ref"           : round(self.s_ref, 3),
            "s_model"         : round(self.s_model, 3),
            "s_orf"           : round(self.s_orf, 3),
            "flag"            : self.flag,
            "is_pseudogene"   : self.is_pseudogene,
            "pseudogene_reason": self.pseudogene_reason,
            "notes"           : self.notes,
        }
