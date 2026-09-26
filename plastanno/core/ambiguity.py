"""One policy for "the evidence does not choose".

Several places pick a single winner from competing candidates: locus selection,
the rRNA cluster choice, the rps12 anchor, ORF start selection, the intron-tRNA
pair. Each used max() over a list, so when the evidence tied the answer came from
the order the list happened to be in — and that order reverses when a genome is
presented on the other strand.

The rule here, in order:

  1. rank by BIOLOGICAL evidence alone;
  2. collapse candidates that describe the same model — they are not a
     disagreement, they are one answer reached twice;
  3. if several distinct models remain at the same evidence rank, that is an
     ambiguity, and it is reported as one:
       * a stable primary is chosen by a digest — the digest is a canonical
         ORDERING device, not evidence, and carries no biological meaning;
       * the other models are kept in full as alternatives;
       * every member shares one ambiguity_id, so a locus where the evidence
         could not decide is ONE review locus and not several features.

That last point matters beyond tidiness: counting each overlapping candidate as
its own NEEDS_REVIEW would inflate the review rate, which is the number H4 is
about.
"""

import hashlib


def resolve(candidates, evidence, digest, model=None):
    """Pick a representation for a set of competing candidates.

    `evidence` returns the biological ranking key, higher being better.
    `digest`   returns a presentation-independent string. It ORDERS and it builds
               the group id; it never contributes to the ranking.
    `model`    returns what makes two candidates the same answer; defaults to the
               digest, i.e. identical sequence.

    Returns (primary, alternatives, ambiguity_id, orderable).

      ambiguity_id is derived from the UNORDERED SET of the group's digests, so it
      is the same whatever order the candidates arrive in and whatever rotation or
      strand the genome was presented on.

      orderable is False when two distinct models share a digest. Then there is no
      semantic primary: nothing presentation-independent can put one first, and a
      caller that emitted only `primary` would emit coordinates that change with
      the input order. Such a caller must emit EVERY member, all carrying the same
      ambiguity_id — which is what makes them one review locus rather than several.
    """
    if not candidates:
        return None, [], None, True
    model = model or digest

    top = max(evidence(c) for c in candidates)
    best = [c for c in candidates if evidence(c) == top]

    # the same model reached twice is one answer, not a disagreement
    seen, distinct = set(), []
    for c in sorted(best, key=digest):
        m = model(c)
        if m in seen:
            continue
        seen.add(m)
        distinct.append(c)

    if len(distinct) == 1:
        return distinct[0], [], None, True

    digests = [digest(c) for c in distinct]
    amb = "AMB" + hashlib.sha256(
        "|".join(sorted(digests)).encode()).hexdigest()[:10]
    orderable = len(set(digests)) == len(digests)
    return distinct[0], distinct[1:], amb, orderable


def describe(alternatives, as_dict):
    """The alternatives in the form Feature.alternatives expects."""
    return [as_dict(a) for a in alternatives]
