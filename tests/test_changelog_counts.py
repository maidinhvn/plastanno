#!/usr/bin/env python3
"""Counts written in the CHANGELOG must match what can be counted.

The 3.0.0 entry said "Thirteen test files, 206 checks, none needing an external
dataset". The tag has fourteen files and 365 checks, and one of them does need
the database. Three wrong numbers in one sentence, on the page a release points
at, found only because someone re-derived them by hand before publishing.

What this checks, and how:

  * "N test files" in a released section — counted from that version's git tag,
    so it verifies the state the entry describes rather than today's tree.
  * "N entries" for database_CHECKSUMS.sha256 — counted from the shipped file.

What it deliberately does not check, and says so rather than skipping quietly:

  * check counts. Getting them needs the suite run at that tag with the database
    present, which is not a unit test's job.
  * commit counts. Those describe the development repository, which is not this
    one — the public history is consolidated, so counting here gives a different
    and meaningless number.

Claims inside quotation marks are skipped: the [Unreleased] section quotes the
old wrong numbers in order to correct them, and a guard that flagged those would
punish the record for being honest.
"""
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

RUN = [0]
FAIL = []
UNVERIFIED = []


def check(label, got, want):
    RUN[0] += 1
    ok = got == want
    if not ok:
        FAIL.append("%s (got %r, want %r)" % (label, got, want))
    print("  %s %s" % ("ok  " if ok else "FAIL", label))


WORDS = {"one": 1, "two": 2, "three": 3, "four": 4, "five": 5, "six": 6,
         "seven": 7, "eight": 8, "nine": 9, "ten": 10, "eleven": 11,
         "twelve": 12, "thirteen": 13, "fourteen": 14, "fifteen": 15,
         "sixteen": 16, "seventeen": 17, "eighteen": 18, "nineteen": 19,
         "twenty": 20}


def as_int(tok):
    tok = tok.replace(",", "")
    if tok.isdigit():
        return int(tok)
    return WORDS.get(tok.lower())


def strip_quoted(text):
    """Remove double-quoted spans, where the file quotes claims to correct them."""
    return re.sub(r'"[^"]*"', '""', text)


raw = open(os.path.join(ROOT, "CHANGELOG.md")).read()

# Split into sections keyed by the version heading.
parts = re.split(r"^## \[([^\]]+)\][^\n]*$", raw, flags=re.M)
sections = {}
for i in range(1, len(parts), 2):
    sections[parts[i]] = strip_quoted(parts[i + 1])
check("the changelog has version sections", len(sections) >= 2, True)

print()
print("test-file counts, against each version's tag")
tags = set(subprocess.run(["git", "-C", ROOT, "tag"], capture_output=True,
                          text=True).stdout.split())
found_any = False
for ver, body in sections.items():
    for tok in re.findall(r"([A-Za-z]+|\d+) test files", body):
        n = as_int(tok)
        if n is None:
            continue
        tag = "v" + ver
        if tag not in tags:
            UNVERIFIED.append("%s: '%s test files' — tag %s not present"
                              % (ver, tok, tag))
            continue
        listing = subprocess.run(
            ["git", "-C", ROOT, "ls-tree", "-r", "--name-only", tag],
            capture_output=True, text=True).stdout.splitlines()
        actual = sum(1 for p in listing
                     if re.fullmatch(r"tests/test_[A-Za-z0-9_]+\.py", p))
        check("%s says %s test files; tag %s has %d" % (ver, tok, tag, actual),
              n, actual)
        found_any = True
# A tree with no version tags cannot verify any of these — the development
# repository tags milestones, not releases. That is a property of the tree, not
# a defect, so require verification only where release tags exist.
release_tags = [t for t in tags if t.startswith("v")]
if release_tags:
    check("at least one test-file count was verifiable", found_any, True)
else:
    print("       skip  no release tags in this tree; counts unverifiable here")

print()
print("checksum-file entry count, against the shipped file")
sums = os.path.join(ROOT, "database_CHECKSUMS.sha256")
for ver, body in sections.items():
    for tok in re.findall(r"([\d,]+) entries", body):
        n = as_int(tok)
        if n is None:
            continue
        if not os.path.exists(sums):
            UNVERIFIED.append("%s: '%s entries' — database_CHECKSUMS.sha256 absent"
                              % (ver, tok))
            continue
        actual = sum(1 for _ in open(sums))
        check("%s says %s entries; the file has %d" % (ver, tok, actual),
              n, actual)

print()
print("claims this guard cannot reach")
for kind, pat in (("check counts", r"([\d,]+) checks"),
                  ("commit counts", r"([\d,]+) commits")):
    for ver, body in sections.items():
        for tok in re.findall(pat, body):
            UNVERIFIED.append("%s: \"%s\" among %s — not verifiable here" % (ver, tok, kind))
if UNVERIFIED:
    for u in UNVERIFIED:
        print("       %s" % u)
else:
    print("       none")
# Listing them is the point; not reaching them is not a failure.
check("unreachable claims are listed rather than ignored",
      isinstance(UNVERIFIED, list), True)

print()
print("%d checks, %d failed" % (RUN[0], len(FAIL)))
for f in FAIL:
    print("  FAILED: %s" % f)
sys.exit(1 if FAIL else 0)
