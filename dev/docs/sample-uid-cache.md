# Per-sample caching via derived sample UIDs

**Status:** proposal — not implemented. All design choices are settled (reviewed in
#69); §11 records the three that were open and how they were decided.
**Depends on:** the Genome/name split (#61, merged; closed #39)

**Decisions taken:**

- The UID tree is scoped to `output_dir`; no shared cache this round (§2).
- SRA is keyed on the accession, not on ENA's free MD5 (§4.2).
- mtime is a re-verification trigger, not a key component (§4.2).
- Per-sample paths carry parameter digests alongside the UID (§5).
- Registration is implicit at Snakefile parse time (§7).
- The transition layer is a relative-symlink farm (§8).
- Tier-2 content hashing runs by default, with a documented opt-out (§11.1).
- A short-prefix UID collision is a hard error (§11.2).
- Nothing is evicted automatically; deletion is manual and documented (§11.3).

## 1. Problem

`Sample_ID` is a sample's identity throughout the pipeline, and it conflates four
different kinds of thing:

| Concern | Nature | Example |
|---|---|---|
| Genome | analysis parameter | `ColCEN` |
| Levels | mutable user metadata | `genotype:WT,tissue:root` |
| Read_layout | inherent library property | `PE` |
| Replicate_ID | inherent, but free text | `rep1` |

Only the sequencing data is invariant. Because outputs are keyed on names derived
from the mutable members of that list, changing the analytical design — adding a
factor, renaming a level, regrouping samples — invalidates per-sample work whose
inputs did not change.

That work dominates. On the hg38 chr21 integration case (33 samples, 428 jobs,
2 h wall), downloads alone account for 3 h 1 m of CPU — 151 % of wall time — and
the single longest job in the run is a 1 h 27 m bismark alignment. The test data
is a chromosome subset; on real libraries the per-sample share grows.

## 2. Goal and non-goals

**Goal.** Per-sample outputs survive any change to Factors, Levels, or analysis
grouping — none of those describe the sample or how it was processed.

Genome selection is deliberately not in that list, because the two halves of the
per-sample phase answer to it differently, along the same seam #61 established:

- **Download, trim and QC are genome-independent** and survive a genome change
  completely. One download and one trim serve every reference the sample is
  mapped to.
- **Alignment is per-genome by definition** and must *not* survive one. Changing
  the reference has to produce a new alignment — that is why `<genome>` is a path
  component and the reference identity is inside `align-digest` (§5.2).

What the cache buys for genomes is therefore additive rather than invariant:
adding a reference to a sample costs exactly one new alignment, with no
re-download and no re-trim, and leaves the alignments to every other reference
untouched.

**Non-goals this round:**

- A lab-shared or cross-project cache. **The UID tree is scoped to `output_dir`.**
  The path is configurable so a shared cache becomes a config change rather than a
  redesign.
- Content-verified identity as the cache *key* (see §4.2).
- Caching of analysis-phase outputs.

## 3. Why not Snakemake's `--cache`

Evaluated against Snakemake 9.21, empirically as well as from source.

Its key is a SHA-256 chain over the raw shell command template, all `params:`
values, full-content SHA-256 of every input not produced by another job, the conda
environment, and recursively the hashes of upstream jobs. Outputs are *moved* into
`$SNAKEMAKE_OUTPUT_CACHE/<sha256><ext>` and replaced with an absolute symlink.

It gets the core idea right — with no sample-identifying params, renaming a sample
is a cache hit and changing input bytes is a miss. Five things rule it out:

1. **Params poison the key.** 53 of 122 rules pass a sample name through `params:`.
   Verified: identical reads, `sample_name` in params, renamed sample → a second
   cache entry. The rename is precisely the case the refactor exists to serve.
2. **The chain bottoms out at a string.** `get_fastq_pe` declares no file input at
   all; the source is `params.seq_id` / `params.fastq_path`. Every provenance chain
   is rooted in an accession or path string, not content.
3. **Silent staleness.** Verified: with the source as a param, editing the source
   file's contents still served the old cached output, with no warning. EPICC's
   local-FASTQ path has exactly this shape.
4. **`temp()` is silently defeated.** Verified: a `temp()` output is cleaned from
   `results/` but retained permanently in the cache. We `temp()` raw FASTQs because
   they are the largest files in the pipeline.
5. **No management surface.** A flat `<sha256>.ext` directory: no per-sample
   navigability, no manifest, no prune command, and entries written `0666` by
   design, since the cache is meant to be shared.

**Verdict:** `--cache` caches on *how a file was made*, rooted at an accession
string. We need caching on *what the data is*, rooted at content. Those coincide
only when the source is a local file declared as a rule input and no param carries
a sample name — neither holds here.

Its key composition is still a good model, and its multi-output constraint
(`multiext()` or all-named outputs) is a useful precedent.

## 4. Sample identity

### 4.1 Probes

Identity is established **before DAG construction** from cheap, protocol-specific
metadata. No source is downloaded to compute its UID.

| Source | Probe | Strength |
|---|---|---|
| SRA accession | ENA `filereport` API → `fastq_md5`, `fastq_bytes`, `read_count` | true content MD5, no download |
| Public S3 | HEAD → ETag (content MD5, or `md5-of-md5s-N` for multipart), `x-amz-version-id` | content-derived |
| Generic HTTPS | HEAD → `Content-Length`, `Last-Modified`, `ETag` | size + mtime proxy |
| Local file | `stat` → realpath, size, mtime | size + mtime proxy |

Verified against real inputs from our own test sheets. `SRR20678305` returns both
mates' MD5s and byte counts in a single GET. ENA and `transfer.lemna.org` return
Apache-style `inode-size-mtime` ETags — **not** content hashes, so generic-HTTP
ETags must be treated as opaque tokens. Public S3 ETags are genuine MD5s.

Servers that omit `Content-Length` (chunked encoding) or reject HEAD (405) need a
`Range: bytes=0-0` GET fallback, which returns `Content-Range: bytes 0-0/TOTAL`.

### 4.2 What is keyed, and what is only recorded

The UID is derived from the probe. Everything else the probe returns is recorded
in the manifest for verification and drift detection, not for keying.

**SRA is keyed on the accession, not the MD5.** This is deliberate and
counterintuitive, since the MD5 is free. `get_fastq_pe` tries ENA first and falls
back to `fasterq-dump`, and the two routes produce **different bytes for the same
accession** — `fasterq-dump` reconstructs from SRA-normalised format and read names
differ. Keying on `fastq_md5` would flip the UID exactly when a download was
hardest. Key on the accession; record the MD5, byte count, and *the route actually
taken*. Accessions are stable in practice, and the recorded MD5 lets registration
detect the rare case where they are not.

**mtime is a re-verification trigger, not a key component.** Keying on
`path + size + mtime` trades silent staleness for spurious invalidation — copying a
project between filesystems, or an `rsync` that does not preserve times, would
orphan the entire cache. Instead: key on `realpath + size`; if mtime has moved
since registration, hash the file once and compare against the recorded content
hash. Match → keep the UID, update the recorded mtime. Mismatch → new UID. A hash
is paid only when something genuinely looks different, and neither failure mode
survives. The same rule applies to `Last-Modified` on HTTP sources.

**Tier-2 content hash.** After a file lands, compute a BLAKE3 digest once and
record it. The data was just written, so the page cache is warm and the cost is
negligible against download plus alignment. This buys a real fingerprint for
provenance, an integrity check for cached artefacts, and the ability to detect two
UIDs holding identical data. This runs by default; see §11.1 for the opt-out.

### 4.3 Canonicalisation

- **PE pairs:** order-dependent. R1/R2 is semantic.
- **`+`-merged components:** order-**in**dependent. Probe each component, sort the
  digests, then combine — the same data listed in a different order must give the
  same UID.
- **Local paths:** `realpath`, absolute. A relative path would make the UID depend
  on the working directory.
- **Missing local file:** hard-fail at registration. Earlier and clearer than
  today's failure.
- **Probe method is recorded per file**, so a future stronger probe can be
  introduced without silently reinterpreting existing manifest entries.

### 4.4 UID form

URL-safe base64 of the digest, truncated to 8 characters = 48 bits. At 10 000
samples the collision probability is ~2 × 10⁻⁷. The full digest is always stored in
the manifest. A short-prefix collision is a hard error; see §11.2 for the remedy.

## 5. Identity is not the cache key

The UID identifies a *sample*. Per-sample outputs also depend on the *computation*,
so paths need both:

```
<output_dir>/.epicc/samples/<UID>/<trim-digest>/trim__R1.fastq.gz
<output_dir>/.epicc/samples/<UID>/<trim-digest>/<genome>/<align-digest>/final.bam
```

Keyed on UID alone, changing `adapter1` or `quality_threshold` would silently reuse
stale trimmed reads, and switching `chip_aligner` would reuse the wrong BAM.

This is the mirror image of what poisons `--cache`: Snakemake hashes *all* params
including the sample name, so nothing survives a rename. We hash params too — but
scoped to those that genuinely affect the step, with sample identity factored out
entirely. **Separating identity from computation, and hashing each explicitly, is
the central principle of this design.**

### 5.1 Membership rule

A parameter belongs in a digest **iff changing it changes the bytes of that step's
output.** Sample names, log and banner text, thread counts, and resource requests do
not qualify — they are exactly what Snakemake's `--cache` gets wrong by hashing
`params:` wholesale.

Per-env config lookups resolve to a *value* before hashing: what enters the digest is
`adapter1=AGATCGGAAGAGC`, not `config['adapter1'][env]`. Two samples in different
envs that resolve to identical trimming settings therefore share a digest, which is
the desired behaviour.

### 5.2 Digest membership

**trim-digest**

- `adapter1`, `adapter2`
- `quality_threshold`
- `min_read_length`
- `trim_front`
- `trimmed_fastqs` — it changes what "raw" means, so it changes the output bytes
- `Read_layout` (SE/PE) — selects a different rule and a different fastp invocation

**align-digest**

- the aligner (`chip_aligner` / `atac_aligner`, or STAR / bismark / ShortStack by env)
- the mapping strategy (`chip_mapping_strategy` / `atac_mapping_strategy`), including
  the documented override that forces bowtie2 for `repeat` / `repeatall`
- the reference identity

The reference is probed with **the same machinery as read sources** (§4.1) — genome
config fields already accept URLs, so `fasta_file` and the index-build inputs go
through the same protocol-specific probes. A reference swapped in place under an
unchanged path is then caught by the same mtime re-verification rule.

`<genome>` stays a *separate* path component rather than being folded into the
digest, so the tree remains navigable by reference name — the readability property
#61 established.

### 5.3 Form and reversibility

Digests use the same construction as the UID: sorted `key=value` pairs, SHA-256,
URL-safe base64, truncated to 8 characters, with the full digest recorded.

**Every digest is recorded in the manifest along with the parameter set that produced
it.** Without that reverse lookup an 8-character directory name is undebuggable, and
"why did this re-run?" — the question this design will be asked most often — becomes
unanswerable. Registration should be able to print the diff between two digests.

## 6. The manifest

Lives at `<output_dir>/.epicc/manifest.json` — under `output_dir`, never under
`repo_folder`, which may be read-only for conda installs where the repo sits at
`$PREFIX/share/epicc/`.

Per entry: UID (short and full), the declared source string, probe method, probe
results (sizes, ETags, MD5s, mtimes, resolved URLs), the download route taken, and
the tier-2 content hash once known.

Written atomically (temp file + `os.replace`), since parse-time registration runs
before Snakemake takes its directory lock. Schema-versioned.

> The one prior parse-time write into `REPO_FOLDER/config/` —
> `analysis_samplefile_<name>.tsv` — has been removed rather than relocated: nothing
> ever read it back. So the manifest will be the only state registration writes.

## 7. Registration

**Implicit, at Snakefile parse time.** Top-level code runs to completion before any
rule is defined or the DAG is built, and it covers raw `snakemake ...` invocations,
which we have committed to keeping functional. There is precedent: the top level
already reads and validates the sheet and writes a TSV as a parse-time side effect.

**Registration cannot be a rule.** `ncbi_taxid` looks like a counterexample — it is
auto-computed over the network — but it is a rule producing a JSON, which works only
because no other rule's *output path* depends on it. A UID determines paths, so it
must exist before rules are declared. A checkpoint would work and would make the
entire downstream tree checkpoint-dependent.

**The constraint that shapes this: every cluster job re-parses the whole Snakefile.**
The SLURM executor spawns `python -m snakemake --snakefile <SF> ... --mode remote
--target-jobs ...` per job. Naive parse-time registration would mean N network probes
times every job in the run. Therefore:

1. **Manifest-first and incremental.** Probe only sources that are absent from the
   manifest or whose declared source string changed. Steady state is one JSON read,
   so dry-runs and re-runs cost nothing.
2. **Read-only under `--mode remote`** (detectable in `sys.argv`). Remote jobs
   require the manifest, never probe, never write — which also removes any write race
   between concurrent jobs.
3. **Atomic write** on the head node.

**Invocation paths:**

- `--unlock` skips registration entirely — extend the existing `UNLOCKING` guard.
- `epicc validate` registers. Catching an unreachable URL before a six-hour run is
  worth the network round-trip; provide a flag to skip.

**`epicc register` remains available** as an explicit pre-flight — not a required
step. This mirrors `epicc validate --build-envs`, which exists because conda env
creation fails inside SLURM allocations at some sites. Same class of problem:
network reachability can differ between login and compute nodes, and some sites
firewall outbound HTTP from compute nodes.

**Failure-mode change to accept.** A network hiccup becomes a startup failure rather
than a job failure; `get_fastq_pe` currently carries `retries: 3`. The incremental
rule contains this — only genuinely new samples touch the network, and those cannot
proceed anyway. Add per-probe timeouts and retries, and an offline mode that requires
an existing manifest.

## 8. Layout and the transition layer

The per-sample phase writes into the UID tree. A transition rule maps those onto the
semantic names used today, so `results/<env>/` keeps its current readable structure.

**The links are relative.** Since the UID tree lives under `output_dir` (§2) and the
semantic names live under `output_dir` too, every link is expressible as a relative
path within one directory tree — so the whole results directory relocates, gets
copied to another filesystem, or moves between a scratch and an archive mount without
breaking. Absolute links would be simpler to generate and would forfeit that; this is
also the one thing Snakemake's `--cache` gets wrong here, since its cache lives at an
arbitrary `$SNAKEMAKE_OUTPUT_CACHE` outside the results tree and its links are
therefore absolute by necessity.

The alternative — downstream rules reading the UID tree directly, with semantic names
reserved for final user-facing outputs — was considered and rejected: it trades a
readable `results/` for fewer links, and readable intermediates are worth more than
the links cost.

Consequences to handle:

- `tar` / `rsync` need `-L` / `--copy-links`, or dangling links get archived. This
  belongs in the docs, and in whatever we say about moving completed runs.
- Snakemake resolves symlink mtimes to the target; combined with `temp()`, reclaiming
  a target leaves a dangling link that still looks present to a casual `ls`.
- Deleting the UID tree guts `results/`. **The UID tree is the real data** — this
  needs a prominent documented warning, and `epicc clean` must not offer to remove it
  as though it were scratch.

## 9. Retention

Retention stays a policy knob; nothing is retained unconditionally. The current
ladder already does the right thing — raw FASTQs are `temp()`, trimmed FASTQs are
`temp()` unless `keep_trimmed_fastqs`, BAMs are kept. So the cache naturally holds
BAMs plus QC and stats, and a re-download happens only when a re-trim is genuinely
needed. This is the outcome `--cache` cannot offer, since caching a rule there
retains its output permanently regardless of `temp()`.

## 10. Sequencing

1. **Genome as a post-alignment path token** — done, #61.
2. **Registration pass, manifest, UIDs — no path changes.** Establish and validate
   identity while everything still writes where it writes today. Independently
   shippable and independently revertible.
3. **Move per-sample processing into the UID tree**, add the transition layer.

§11.1 and §11.2 govern step 2; §11.3 governs step 3.

## 11. Decided choices

These three were the open questions at review. Each is recorded with the options and
the trade-off that produced the decision, so the reasoning survives the decision.

### 11.1 Is tier-2 content hashing on by default, or opt-in?

**Decided: on by default, with a documented opt-out flag.**

After a file lands, we can compute a BLAKE3 digest once and record it (§4.2). The
question is whether that happens automatically.

- **On by default.** Every cached artefact gets a real fingerprint. Enables
  `epicc cache verify`, makes the mtime re-verification rule (§4.2) usable — that
  rule *needs* a recorded content hash to compare against, so without tier-2 the
  trigger degrades to "assume changed" — and lets registration detect two UIDs
  holding identical data. Costs one pass over data just written, with a warm page
  cache; negligible against download plus alignment, but not free for a 50 GB WGBS
  pair on a busy shared filesystem.
- **Opt-in.** Zero added I/O by default. The mtime trigger has nothing to compare
  against, so a moved mtime must either be ignored (silent staleness returns) or
  force a new UID (spurious invalidation returns) — reintroducing exactly the
  failure modes §4.2 was designed to eliminate.

The mtime rule is load-bearing and does not work without a recorded hash to compare
against, so opt-in would quietly reintroduce the failure modes §4.2 exists to
eliminate. The opt-out is for sites where the filesystem is the bottleneck. Note the
cost is bulk sequential read, not metadata churn — the cheaper of the two on a shared
parallel filesystem.

### 11.2 What happens on a short-prefix collision?

**Decided: hard error only.** No automatic migration.

8 characters is 48 bits (§4.4). Collisions are vanishingly unlikely at lab scale but
must not corrupt data if they happen.

- **Hard error only.** Registration fails, names both samples, and documents a knob
  to lengthen the prefix. Simple, obvious, and forces a human decision. Cost: a run
  is blocked until someone intervenes, and lengthening the prefix renames every path
  in the UID tree.
- **Automatic prefix lengthening.** Registration extends the prefix for the colliding
  pair (or globally) and records the new length in the manifest. Never blocks. Cost:
  path length becomes a function of manifest history, so two projects with the same
  samples can have different paths, and the migration has to be implemented and
  tested for something that may never fire.

This matches git's short-SHA behaviour. Automatic migration is real complexity
guarding a ~2 × 10⁻⁷ event, and the failure is loud rather than silent either way.

**What the user does about it.** A hard error is only acceptable if the way out is
obvious, so the error must carry the whole remedy: both `Sample_ID`s, both source
strings, and both *full* digests — which also demonstrates the data really does
differ, rather than this being a duplicate row. Two genuinely identical inputs listed
twice are a different problem, already caught by the existing cross-row duplicate
`Read_files` check.

Recovery is a mechanical rename, not a re-run. Raise `uid_prefix_length` (default 8)
in the options file; because the manifest stores full digests, registration can
compute both the old and new prefix for every entry and rename the existing UID
directories in place. No re-download, no re-alignment, no recomputed hashes. That is
what makes the strict choice cheap: the loud failure costs a config edit and a
directory rename, and the implementation is a rename loop rather than a migration
framework.

### 11.3 How does the UID tree get evicted?

**Decided: nothing is evicted automatically.** Manual deletion, documented.

The tree grows monotonically. Nothing in this design removes anything, and the whole
point is that artefacts outlive the analysis that produced them.

- **No eviction.** Simplest. The tree is the durable record; users manage disk
  themselves. Cost: unbounded growth on shared storage, and no way to answer "what is
  safe to delete?"
- **UID-aware `epicc clean`.** A mode that reports what is in the tree and removes
  entries no longer referenced by the current sample sheet. `epicc clean --tmp` now
  gives us a precedent for the shape: two independent safety gates, refuse when
  liveness is unknowable rather than assuming dead, and report what was skipped and
  why. Cost: "unreferenced" is a dangerous predicate — a user switching sheets between
  projects in one `output_dir` would see their cache deleted, which is exactly the
  work this design exists to preserve.
- **Age-based reclamation.** Remove entries untouched for N days regardless of
  references. Predictable, but silently deletes the artefacts of a long-running
  project between analysis rounds — again the case the design targets.

Every automatic rule considered deletes exactly the work this design exists to
preserve, so the tree keeps all records. What ships alongside it is reporting only:
an `epicc cache status` giving total size and a per-UID breakdown, so a user deleting
by hand can see what they are removing and what it cost to produce. The docs must
also carry the §8 warning that the UID tree is the real data and `results/` is
symlinks into it.

Revisit once there is real usage data on how large the tree gets. Deleting cached
work is far worse than keeping too much of it, and a reporting-only mode costs almost
nothing to ship.

## 12. Risks

- **Startup-time network dependency** (§7) — mitigated by incrementality, but it is a
  genuine change in failure character.
- **Silent staleness on weak probes.** Generic HTTP and local files are keyed on a
  size proxy. The mtime re-verification rule (§4.2) covers the realistic cases, and
  relies on the tier-2 hash of §11.1; a same-size in-place edit with a preserved mtime would not be
  caught.
- **ENA/`fasterq-dump` route divergence** (§4.2) — recorded, not prevented.
- **Scope creep into a shared cache.** Explicitly out of scope; the permission and
  concurrency model that would require is visible in what `--cache` had to do.
