"""tRNA identity classification, policy version 2.

WHY A SECOND POLICY EXISTS
    `benchmark/gene_naming.classify_trna` scans /gene, /product, /standard_name,
    /gene_synonym AND /note with equal standing. Free-text provenance therefore
    participates in deciding identity: a note reading "renamed trnS-CGA -> trnG-UCC"
    put a second amino acid into the evidence and the locus came back
    "family disagreement" although its structured qualifiers were consistent.

    Rewording the note only made the old text stop matching the regex. The raw
    disagreement was still recorded, and any other phrasing containing a gene name
    would reintroduce the problem. The contract is what has to change.

POLICY v2
    1. STRUCTURED QUALIFIERS DECIDE. /gene, /product and /anticodon are the
       evidence. When any of them resolves the family, /note cannot change the
       answer -- not the family, not the anticodon, not the role.
    2. /note is consulted ONLY when no structured qualifier resolves anything.
       A result reached that way is marked NOTE_DERIVED and is never FULL.
    3. Nothing is invented. An unresolved axis stays unresolved.

    gene_naming.classify_trna is UNTOUCHED: the frozen H4 scorer depends on it and
    its outputs are part of a frozen record.
"""
import re

RESOLVED = "RESOLVED"
FAMILY_ONLY = "FAMILY_ONLY"
UNRESOLVED = "UNRESOLVED"
CONFLICT = "CONFLICT"
NOTE_DERIVED = "NOTE_DERIVED"

POLICY_VERSION = 2

_STRUCTURED = ("gene", "product", "standard_name", "gene_synonym")

AA3 = {"ala": "A", "arg": "R", "asn": "N", "asp": "D", "cys": "C", "gln": "Q",
       "glu": "E", "gly": "G", "his": "H", "ile": "I", "leu": "L", "lys": "K",
       "met": "M", "phe": "F", "pro": "P", "ser": "S", "thr": "T", "trp": "W",
       "tyr": "Y", "val": "V", "sec": "U"}

_RE_TRN = re.compile(r"\btrn(fm|f?[a-z])(?:[-_( ]\s*([acgtu]{3})\b)?", re.I)
_RE_AA = re.compile(r"\bt?rna[- ]?(f?[a-z]{3})\b", re.I)
_RE_ANTI_SEQ = re.compile(r"seq:\s*([acgtu]{3})", re.I)
_RE_ANTI_AA = re.compile(r"aa:\s*([A-Za-z]{3})", re.I)
_INIT = re.compile(r"\bf-?met\b|\bformyl|\binitiator\b|\btrnfm\b", re.I)
_ELONG = re.compile(r"\belongator\b|\btrnm[-_( ]", re.I)
_RE_ATTACHED = re.compile(
    r"(?:trn[a-z]{1,3}|rna-[a-z]{3})[-_( ]+\s*([acgtu]{3})\b", re.I)


def _scan(values):
    """(families, triplets, initiator, elongator) from a list of strings.

    PRECEDENCE MATTERS. "tRNA-Gly" contains the substring "tRN" followed by "A",
    so a bare trn-pattern reads it as family A -- the trnA / tRNA-Ala collision
    the original classifier documents. An amino-acid LABEL is therefore tried
    first, and only a value with no such label is read as a trnX name.
    """
    fams, tris, init, elong = set(), set(), False, False
    for v in values:
        if not v:
            continue
        init = init or bool(_INIT.search(v))
        elong = elong or bool(_ELONG.search(v))
        fam_here = None
        m = _RE_AA.search(v)
        if m:
            tag = m.group(1).lower()
            base = tag[1:] if tag.startswith("f") and tag[1:] in AA3 else tag
            if base in AA3:                     # a real amino-acid label
                fam_here = AA3[base]
                if tag.startswith("f"):
                    init = True
        if fam_here is None:
            for m2 in _RE_TRN.finditer(v):
                tag = m2.group(1).lower()
                if tag == "fm":
                    fam_here, init = "M", True
                elif len(tag) == 1 and tag.upper() in AA3.values():
                    fam_here = tag.upper()
                if m2.group(2):
                    tris.add(m2.group(2).upper().replace("T", "U"))
        else:
            # the label fixed the family; a triplet may still be attached to it
            for m3 in _RE_ATTACHED.finditer(v):
                tris.add(m3.group(1).upper().replace("T", "U"))
        if fam_here:
            fams.add(fam_here)
    return fams, tris, init, elong


def classify(quals):
    """Policy-v2 identity of a tRNA feature.

    Returns a dict with family, anticodon, role, status and evidence_source.
    """
    out = {"family": None, "anticodon": None, "role": None,
           "status": UNRESOLVED, "evidence_source": None, "conflicts": []}

    structured = []
    for k in _STRUCTURED:
        structured.extend(quals.get(k, []) or [])
    fams, tris, init, elong = _scan(structured)

    # /anticodon is structured and authoritative for the triplet when it parses
    for v in quals.get("anticodon", []) or []:
        m = _RE_ANTI_SEQ.search(v or "")
        if m:
            tris = {m.group(1).upper().replace("T", "U")}
        m = _RE_ANTI_AA.search(v or "")
        if m and m.group(1).lower() in AA3:
            fams = {AA3[m.group(1).lower()]}

    source = "STRUCTURED"
    if not fams and not tris:
        # 2. only now may /note speak, and the answer is marked as such
        fams, tris, init, elong = _scan(quals.get("note", []) or [])
        source = "NOTE" if (fams or tris) else None

    if len(fams) > 1:
        out["conflicts"].append("family disagreement: " + ", ".join(sorted(fams)))
        out["status"] = CONFLICT
        out["evidence_source"] = source
        return out
    if fams:
        out["family"] = next(iter(fams))
    if len(tris) > 1:
        out["conflicts"].append("multiple anticodons: " + ", ".join(sorted(tris)))
    elif tris:
        out["anticodon"] = next(iter(tris))

    if out["family"] == "M":
        if init and elong:
            out["conflicts"].append("both initiator and elongator markers")
        elif init:
            out["role"] = "INITIATOR"
        elif elong:
            out["role"] = "ELONGATOR"

    if out["family"] is None:
        out["status"] = UNRESOLVED
    elif source == "NOTE":
        out["status"] = NOTE_DERIVED
    elif out["anticodon"] is None or (out["family"] == "M" and out["role"] is None):
        out["status"] = FAMILY_ONLY
    else:
        out["status"] = RESOLVED
    out["evidence_source"] = source
    return out
