"""Post-write verification: the run only succeeded if the files did.

The pipeline used to report success on the basis that write_all() returned without
raising. It once returned five files instead of seven — two writer calls had been
deleted — and still exited 0, so a batch of genomes was recorded as annotated with
no GenBank file anywhere. Existence is not enough either: a writer that fails
part-way leaves a truncated file that a downstream parser rejects. So each expected
output is checked for presence, for content, and for being readable in its own
format.
"""

from pathlib import Path


class OutputError(RuntimeError):
    """Raised when the outputs of a run are missing, empty or unparseable."""


# suffix -> (human name, whether it must parse as its own format)
REQUIRED = [
    (".gb",                    "GenBank flatfile"),
    (".gff3",                  "GFF3"),
    (".tbl",                   "NCBI feature table"),
    (".faa",                   "protein FASTA"),
    (".ffn",                   "CDS nucleotide FASTA"),
    (".frn",                   "RNA nucleotide FASTA"),
    (".report",                "text report"),
    (".trna_alternatives.tsv", "tRNA source-merge ledger"),
]

# A genome may legitimately contain no feature of a given kind, so an empty FASTA
# is not by itself a failure — but the file must still exist.
MAY_BE_EMPTY = {".faa", ".ffn", ".frn", ".trna_alternatives.tsv"}


def _check_genbank(path):
    from Bio import SeqIO
    # SeqIO.parse given a path leaves the handle open until the iterator is
    # exhausted or garbage-collected; opening it here keeps the checker from
    # leaking a descriptor per genome in a batch.
    with open(path) as fh:
        recs = list(SeqIO.parse(fh, "genbank"))
    if not recs:
        raise OutputError("%s: no GenBank record could be read" % path.name)
    if len(recs[0].seq) == 0:
        raise OutputError("%s: the record carries no sequence" % path.name)
    return "%d record(s), %d feature(s)" % (len(recs), len(recs[0].features))


def _check_gff3(path):
    """Structural check of every data line: nine tab-separated columns, integer
    start <= end, a legal strand and phase. Written out rather than delegated
    because a GFF3 library is not a dependency of this package."""
    n = 0
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            col = line.split("\t")
            if len(col) != 9:
                raise OutputError("%s:%d: %d columns, expected 9"
                                  % (path.name, lineno, len(col)))
            try:
                st, en = int(col[3]), int(col[4])
            except ValueError:
                raise OutputError("%s:%d: non-integer coordinates %r %r"
                                  % (path.name, lineno, col[3], col[4]))
            if st < 1 or en < st:
                raise OutputError("%s:%d: bad interval %d..%d (GFF3 is 1-based, "
                                  "ascending, inclusive)" % (path.name, lineno, st, en))
            if col[6] not in ("+", "-", ".", "?"):
                raise OutputError("%s:%d: bad strand %r" % (path.name, lineno, col[6]))
            if col[7] not in (".", "0", "1", "2"):
                raise OutputError("%s:%d: bad phase %r" % (path.name, lineno, col[7]))
            n += 1
    if n == 0:
        raise OutputError("%s: no feature lines" % path.name)
    return "%d feature line(s)" % n


def _check_tbl(path):
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh]
    if not lines or not lines[0].startswith(">Feature "):
        raise OutputError("%s: missing the '>Feature <seqid>' header" % path.name)
    n = 0
    for lineno, line in enumerate(lines[1:], 2):
        if not line:
            continue
        col = line.split("\t")
        if line.startswith("\t\t\t"):
            continue                       # qualifier line
        if len(col) < 2:
            raise OutputError("%s:%d: %r is neither an interval nor a qualifier"
                              % (path.name, lineno, line))
        for v in col[:2]:
            if not v.lstrip("<>").isdigit():
                raise OutputError("%s:%d: non-numeric coordinate %r"
                                  % (path.name, lineno, v))
        if len(col) >= 3 and col[2]:
            n += 1
    if n == 0:
        raise OutputError("%s: no features" % path.name)
    return "%d feature(s)" % n


def _check_fasta(path):
    n = 0
    with open(path) as fh:
        first = fh.readline()
        if first and not first.startswith(">"):
            raise OutputError("%s: does not begin with a FASTA header" % path.name)
        if first:
            n = 1 + sum(1 for l in fh if l.startswith(">"))
    return "%d record(s)" % n


def _check_tsv(path):
    """The ledger must match its declared schema exactly.

    Checking only that every row has as many columns as the header accepted a
    two-column file as valid while the real schema has fourteen — a writer that
    silently stopped emitting most of the evidence would have passed.
    """
    from .writers import ALT_COLUMNS
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh]
    rows = [l for l in lines if not l.startswith("#")]
    if not rows:
        raise OutputError("%s: no header row" % path.name)
    head = rows[0].split("\t")
    if head != ALT_COLUMNS:
        missing = [c for c in ALT_COLUMNS if c not in head]
        extra = [c for c in head if c not in ALT_COLUMNS]
        raise OutputError(
            "%s: header does not match the declared schema (%d columns, expected "
            "%d%s%s)" % (path.name, len(head), len(ALT_COLUMNS),
                         "; missing " + ",".join(missing) if missing else "",
                         "; unexpected " + ",".join(extra) if extra else ""))
    for lineno, line in enumerate(rows[1:], 2):
        if line and len(line.split("\t")) != len(ALT_COLUMNS):
            raise OutputError("%s: row %d has %d columns, expected %d"
                              % (path.name, lineno, len(line.split("\t")),
                                 len(ALT_COLUMNS)))
    return "%d row(s), %d columns" % (len(rows) - 1, len(ALT_COLUMNS))


def _check_text(path):
    with open(path) as fh:
        return "%d line(s)" % sum(1 for _ in fh)


PARSERS = {
    ".gb": _check_genbank, ".gff3": _check_gff3, ".tbl": _check_tbl,
    ".faa": _check_fasta, ".ffn": _check_fasta, ".frn": _check_fasta,
    ".report": _check_text, ".trna_alternatives.tsv": _check_tsv,
}


def _feature_intervals(out_dir, seq_id):
    """The same annotation as each format states it, for cross-checking.

    Returns {format: {(gene, type): sorted tuple of 0-based half-open intervals}}.
    Every format is read with its own rules — GenBank via Biopython, GFF3 and the
    feature table by parsing the text — so agreement between them is evidence, not
    a tautology.
    """
    from Bio import SeqIO
    out = {}

    with open(out_dir / (seq_id + ".gb")) as fh:
        rec = next(SeqIO.parse(fh, "genbank"))
    gb = {}
    for f in rec.features:
        if f.type in ("CDS", "tRNA", "rRNA"):
            key = (f.qualifiers.get("gene", ["?"])[0], f.type)
            gb.setdefault(key, []).append(
                tuple(sorted((int(p.start), int(p.end)) for p in f.location.parts)))
    out["gb"] = {k: sorted(v) for k, v in gb.items()}

    # A GFF3 product line carries only Parent=<gene id>; the gene name lives on the
    # gene line it points at, so the gene lines are read first.
    gff, cur, gene_name = {}, {}, {}
    with open(out_dir / (seq_id + ".gff3")) as fh:
        lines = [l.rstrip("\n").split("\t") for l in fh if not l.startswith("#")]
    for c in lines:
        if len(c) == 9 and c[2] == "gene":
            at = dict(kv.split("=", 1) for kv in c[8].split(";") if "=" in kv)
            if "ID" in at:
                gene_name[at["ID"]] = at.get("Name", at["ID"])
    for c in lines:
        if len(c) != 9 or c[2] not in ("CDS", "tRNA", "rRNA"):
            continue
        at = dict(kv.split("=", 1) for kv in c[8].split(";") if "=" in kv)
        ident = at.get("Parent") or at.get("ID", "?")
        name = gene_name.get(ident, at.get("Name", ident))
        cur.setdefault((name, c[2], ident), []).append((int(c[3]) - 1, int(c[4])))
    for (name, typ, _), ivs in cur.items():
        gff.setdefault((name, typ), []).append(tuple(sorted(ivs)))
    out["gff3"] = {k: sorted(v) for k, v in gff.items()}

    tbl, gene, key = {}, None, None
    pending = []
    with open(out_dir / (seq_id + ".tbl")) as fh:
        tbl_lines = fh.readlines()
    for line in tbl_lines:
        line = line.rstrip("\n")
        if line.startswith(">Feature") or not line:
            continue
        c = line.split("\t")
        if line.startswith("\t\t\t"):
            if len(c) >= 5 and c[3] == "gene":
                gene = c[4]
            continue
        if len(c) >= 3 and c[2]:
            if key and pending:
                tbl.setdefault(key, []).append(tuple(sorted(pending)))
            pending = []
            key = None
            if c[2] in ("CDS", "tRNA", "rRNA"):
                key = (gene, c[2])
        if key is not None and len(c) >= 2:
            a = int(c[0].lstrip("<>")); b = int(c[1].lstrip("<>"))
            lo, hi = (a, b) if a <= b else (b, a)
            pending.append((lo - 1, hi))
    if key and pending:
        tbl.setdefault(key, []).append(tuple(sorted(pending)))
    out["tbl"] = {k: sorted(v) for k, v in tbl.items()}
    return out


def cross_check_formats(out_dir, seq_id):
    """Every format must describe the same features at the same coordinates.

    Each writer builds its own location from the same Feature, so a bug in one of
    them shows up here and nowhere else: the run still produces three complete,
    individually valid files that disagree with each other.
    """
    out_dir = Path(out_dir)
    got = _feature_intervals(out_dir, seq_id)
    problems = []
    ref = got["gb"]
    for fmt in ("gff3", "tbl"):
        other = got[fmt]
        for key in sorted(set(ref) | set(other), key=lambda k: (k[0], k[1])):
            a, b = ref.get(key, []), other.get(key, [])
            if a != b:
                problems.append(
                    "%s %s: .gb has %d copy(ies) %s, .%s has %d %s"
                    % (key[0], key[1], len(a), a[:2], fmt, len(b), b[:2]))
    if problems:
        raise OutputError("%d disagreement(s) between output formats for %s:\n  - %s"
                          % (len(problems), seq_id, "\n  - ".join(problems[:40])))
    return "%d locus group(s) agree across .gb, .gff3 and .tbl" % len(ref)


def verify_outputs(out_dir, seq_id, verbose=True):
    """Check every required output of one genome. Returns a {suffix: detail} map.

    Raises OutputError listing every problem found, not just the first — a run that
    lost two writers should say so in one message.
    """
    out_dir = Path(out_dir)
    problems, detail = [], {}
    for suffix, human in REQUIRED:
        path = out_dir / ("%s%s" % (seq_id, suffix))
        if not path.exists():
            problems.append("%s (%s) was not written" % (path.name, human))
            continue
        size = path.stat().st_size
        if size == 0 and suffix not in MAY_BE_EMPTY:
            problems.append("%s (%s) is empty" % (path.name, human))
            continue
        try:
            detail[suffix] = PARSERS[suffix](path) if size else "empty (allowed)"
        except OutputError as e:
            problems.append(str(e))
        except Exception as e:
            problems.append("%s (%s) could not be parsed: %s"
                            % (path.name, human, e))
    if problems:
        raise OutputError("%d output problem(s) for %s:\n  - %s"
                          % (len(problems), seq_id, "\n  - ".join(problems)))
    detail["cross-format"] = cross_check_formats(out_dir, seq_id)
    if verbose:
        for suffix, _ in REQUIRED:
            print("      verified %-24s %s" % (seq_id + suffix, detail.get(suffix, "")))
        print("      verified %-24s %s" % ("cross-format", detail["cross-format"]))
    return detail
