# Updating the reference database on Zenodo

The database is not shipped in the wheel. `plastanno fetch-db` downloads it from
Zenodo, checks its MD5 and extracts it. So a database change needs **three**
things kept in step, and getting any one wrong is silent until a user runs the
tool:

1. a new Zenodo version of the deposit,
2. `URL` in `plastanno/fetch_db.py` pointing at that version's file,
3. `MD5` in the same file matching the new tarball.

If the MD5 is not updated, `fetch-db` aborts with a checksum mismatch. If only
the MD5 is updated, users download the old bundle and the check fails. Both must
change together.

## What changed in this update

Eight genes that the previous bundle lacked:

    chlB  chlL  chlN  lhbA  pafI  pafII  rpl21  ycf66

Each contributes two files — an HMM profile in `hmm_db/hmm_profiles/` and a
protein FASTA in `protein_db/` — plus the regenerated `hmm_db/all_profiles.hmm`
that concatenates them. Nothing else in the database differs: `blast_db`,
`trna_db`, `rrna_db`, `exon_db` and `boundary_db` are byte-identical to the
previous bundle.

Verified before packaging: the tree matches the frozen snapshot
`db_snapshots/v2026-09-24_final` on **110 of 110** manifest entries, and
`protein_db` has the same 403 files as the tested tree. A single-genome run
produces output byte-identical to the tested tree across all seven output files.

`database_CHECKSUMS.sha256` was regenerated: 5634 entries, all verifying. It now
also covers `boundary_db/*` and `exon_templates.json`, which ship in the
repository but were missing from the previous checksum list. The six `*.bak`
entries were dropped — those are development snapshots and are not part of the
package.

## Step 1 — the tarball

Already built at `/data06/users/vutrinh/plastanno-database.tar.gz`. Its archive
root contains `database/`, which is what `fetch_db` requires: it extracts into
the platform data directory and then expects `<parent>/database/blast_db` to
exist.

To rebuild it from scratch:

    cd <repo root>          # the directory that CONTAINS database/
    tar czf /tmp/plastanno-database.tar.gz database
    md5sum /tmp/plastanno-database.tar.gz

Build it from the directory above `database/`, never from inside it — a tarball
whose root is `blast_db/`, `hmm_db/`, … extracts to the wrong place and
`fetch-db` then reports the database as absent.

## Step 2 — publish a new version on Zenodo

The deposit has a **concept DOI** (`10.5281/zenodo.20807994`) that always points
at the newest version, and a per-version record (`20807995` for the current one).
Publishing a new version keeps the concept DOI and mints a new record.

1. Open the concept DOI, or the current record `20807995`.
2. Choose **New version** (not a new upload — a new upload would create a
   separate deposit with its own concept DOI and break the link to earlier
   versions).
3. Remove the old `plastanno-database.tar.gz` from the file list and upload the
   new one. Zenodo does not allow replacing a file in place.
4. In the description, note what changed: the eight genes above.
5. **Publish.** Note the new record number from the URL —
   `https://zenodo.org/records/<NEW_RECORD>`.

## Step 3 — point the tool at it

Two constants in `plastanno/fetch_db.py`:

    URL = "https://zenodo.org/records/<NEW_RECORD>/files/plastanno-database.tar.gz"
    MD5 = "<md5 of the new tarball>"

Then verify end to end, on a machine where the database is not already present,
or with `--force`:

    plastanno fetch-db --force

It must print the download, pass the checksum step, and finish with
`Database installed at …`. A checksum failure here means the MD5 and the uploaded
file disagree; do not edit the MD5 to match whatever was downloaded without first
confirming the upload is the tarball you built.

## Step 4 — confirm the installed copy is the tested one

After `fetch-db`, from inside the installed `database/`:

    sha256sum -c ../database_CHECKSUMS.sha256

Expect 5634 lines OK. `database_CHECKSUMS.sha256` is tracked in the repository,
so a user can run this check themselves; that is the point of shipping it.

## Why this matters for the release

Without step 2, the code and the database go out of step: the tool would ship
expecting eight genes the downloadable bundle does not contain. The failure is
not loud — those genes would simply never be annotated — which is exactly the
kind of discrepancy that is hard to notice and hard to explain later.
