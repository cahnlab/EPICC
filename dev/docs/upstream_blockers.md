# Upstream Blockers

Workarounds in this repo that exist solely because of an upstream constraint
(broken release, hardcoded behavior, naming policy, missing feature). Each
entry should be cheap to revisit: it lists what to remove, where it lives,
and a one-line check that will start failing once the upstream catches up.

When you fix or revisit one of these, **delete the entry** — this file is
not a changelog. If a check turns out to be ambiguous or expensive, simplify
it before merging the entry; otherwise the file will rot.

## How to use

1. When you upgrade a pinned package or notice a related upstream release,
   skim this file and run the **Check** for any entry that touches what
   you changed.
2. If the check passes, lift the workaround and remove the entry.
3. If you add a new workaround, append an entry below following the same
   shape (Constraint / Workaround / Upstream / Check / Remove when).

---

## `setuptools<71` in `epibutton_chip.yaml`

- **Constraint:** `epic2` (broad-peak caller for CUT&RUN) imports
  `pkg_resources` at startup. setuptools 71 unbundled `pkg_resources` into
  a separate PyPI distribution (`pkg_resources` package), and conda-forge
  has not yet shipped it, so any setuptools >=71 in the env breaks `epic2`.
- **Workaround:** `setuptools<71` pin in `workflow/envs/epibutton_chip.yaml`.
- **Upstream:**
  - <https://github.com/pypa/setuptools/issues/4519> (deprecation/removal of
    bundled `pkg_resources`)
  - <https://github.com/biocore-ntnu/epic2> (no release that drops the
    `pkg_resources` import yet)
- **Check:** in a fresh chip env without the pin,
  `python -c "import epic2"` succeeds.
- **Remove when:** epic2 cuts a release that imports `importlib.resources`
  / `importlib.metadata` instead, OR conda-forge ships the standalone
  `pkg_resources` package.

## `numpy<=1.19` in `epibutton_idr.yaml`

- **Constraint:** the `idr` package (Irreproducible Discovery Rate) uses
  `np.int` / `np.float` / `np.bool` aliases that numpy 1.20 deprecated and
  1.24 removed. Without the pin, `idr` imports fail with
  `AttributeError: module 'numpy' has no attribute 'int'`.
- **Workaround:** `numpy<=1.19` pin in `workflow/envs/epibutton_idr.yaml`.
  The env is single-purpose (IDR + bedtools + meme) so the old numpy is
  only loaded for IDR-calling rules, not the rest of the pipeline.
- **Upstream:** <https://github.com/nboley/idr> — repo is dormant; the
  scalar-alias deprecations have been raised in issues but no release.
- **Check:** in a fresh idr env without the pin,
  `python -c "import idr"` succeeds.
- **Remove when:** a maintained fork of IDR ships a numpy-1.24-compatible
  release, or we replace `idr` with an alternative caller.

## DMRcaller `parallel=TRUE, BPPARAM=...` plumbing

- **Constraint:** DMRcaller >= 1.42 added a new `parallel=FALSE` default
  to `computeDMRs` that hardcodes `BPPARAM = SerialParam()` regardless
  of `cores=` or any globally-registered BiocParallel backend. The old
  `cores=` parameter is silently ignored on those versions, so a job that
  requests 16 CPUs runs single-threaded.
- **Workaround:** `workflow/scripts/R_call_DMRs.R` registers
  `MulticoreParam(workers=threads)` AND passes
  `parallel=TRUE, BPPARAM=MulticoreParam(workers=threads)` explicitly when
  those formals are present (introspected via `formals()` so the script
  still works on DMRcaller 1.38 and earlier, which used `parallel::mclapply`
  directly via `cores=`).
- **Upstream:**
  - GitHub issue tracking the regression in this repo:
    <https://github.com/cahnlab/EPICC/issues/23>
  - DMRcaller upstream:
    <https://bioconductor.org/packages/DMRcaller/> (no flag has been added
    to make registered BPPARAM the default again).
- **Check:** with DMRcaller pinned to a future release, run a `dmC` test
  with `--cores 8` and confirm `call_DMRs_pairwise` consumes >1 CPU
  (e.g. `top` shows ~700% during the rule, or `time` shows
  `user >> real`).
- **Remove when:** DMRcaller defaults `BPPARAM` to `bpparam()` again (i.e.
  honors the registered backend), at which point the explicit pair can
  collapse back to a one-line `register(MulticoreParam(...))` call.

## Per-genome `.libPaths()` for AnnotationForge GO packages

- **Constraint:** `AnnotationForge::makeOrgPackage()` enforces a strict
  package name format `org.<G><species>.eg.db` (one initial of the genus
  + species). For projects that use multiple reference genomes of the
  *same species* (e.g. ColCEN and TAIR10 — both `org.Athaliana.eg.db`)
  this collides in the conda env's shared R library: the second build
  silently overwrites the first.
- **Workaround:** `workflow/scripts/R_build_GO_database.R` and
  `workflow/scripts/R_GO_analysis.R` install/load each genome's GO
  package into its own per-genome library directory under
  `genomes/<refgenome>/GO/`, with `.libPaths(c(getwd(), .libPaths()))`
  prepended so `library(dbname, character.only=TRUE)` resolves to the
  right one.
- **Upstream:** AnnotationForge has not relaxed the dbname constraint;
  the request would be a `name=` argument to `makeOrgPackage()` that
  bypasses the auto-derived `org.<G><species>.eg.db`.
- **Check:** in a fresh dev env, build GO databases for two same-species
  genomes back-to-back without the per-genome lib dance and confirm both
  `library(org.Athaliana.eg.db)` invocations return the genome-specific
  GAF table.
- **Remove when:** `makeOrgPackage()` exposes a `package_name=` (or
  similar) argument that lets us namespace by reference genome, OR we
  switch off AnnotationForge to a tool with no name-collision policy.

## Bismark 3.x aligner `--nucleotide_coverage` not implemented

- **Constraint:** Bismark 3.x (the Rust rewrite of the suite) recognises the
  aligner's `--nucleotide_coverage` flag but no-ops it — the log prints
  "these options are recognised but not yet active in this build (wired in a
  later phase)" and no `*.nucleotide_stats.txt` is written. Passing the
  resulting (nonexistent) file to `bismark2report --nucleotide_report` then
  crashes it ("os error 2"). Perl Bismark ≤0.25 produced this file at align
  time.
- **Workaround:** `workflow/rules/mC.smk` drops `--nucleotide_coverage` from
  the `bismark` align call and instead runs a separate `bam2nuc` step on the
  deduplicated BAM (bam2nuc *is* implemented in 3.x). Output is named
  `{PE,SE}__<sample>.deduplicated.nucleotide_stats.txt` (after the BAM, not the
  align-time `..._bismark_bt2_pe.nucleotide_stats.txt`), declared as a rule
  output and fed to `bismark2report --nucleotide_report`.
- **Upstream:** <https://github.com/FelixKrueger/Bismark> (Rust suite; the flag
  is staged "for a later phase").
- **Check:** with a future Bismark, run `bismark ... --nucleotide_coverage`
  on a small mC sample and confirm it writes `*_bismark_bt2_pe.nucleotide_stats.txt`
  at align time (no "not yet active" note in the log).
- **Remove when:** the 3.x aligner implements `--nucleotide_coverage`, at which
  point the separate `bam2nuc` step can be dropped and the align-time flag +
  original nucleotide-report path restored.

## `snakemake-executor-plugin-slurm-jobstep<0.6.1`

- **Constraint:** jobstep 0.6.1 sanitizes the environment it hands the nested
  `srun` (`__init__.py`, `popen_env`): every `SLURM_*` variable not on a
  nine-name `keep_slurm` allow-list is dropped. `SLURM_CONF` is not on the
  list. On clusters where `SLURM_CONF` is the only route to `slurm.conf` —
  no `/etc/slurm/slurm.conf` on the compute nodes, path supplied by an
  environment module, as on Bright Cluster Manager sites — the `srun` falls
  through to a configless DNS SRV lookup and every rule dies:

  ```
  srun: error: resolve_ctls_from_dns_srv: res_nsearch error: Unknown host
  srun: error: fetch_config: DNS SRV lookup failed
  srun: error: _establish_config_source: failed to fetch config
  srun: fatal: Could not establish a configuration source
  ```

  This is a *separate* bug from the `delete_slurm_environment()` entry below,
  and worse: that one only fires when snakemake itself runs inside an
  allocation, this one fires on every run regardless of where it was launched,
  and the profile's `precommand` cannot rescue it — the strip happens inside
  the job, downstream of anything the job shell exports. The same allow-list
  drops `SLURM_MEM_PER_CPU`, which makes steps fall back to `DefMemPerCPU` and
  get OOM-killed (upstream issue #53).
- **Workaround:** pin `snakemake-executor-plugin-slurm-jobstep<0.6.1` in
  `config/epicc-env.txt`, `pyproject.toml` and `conda-recipe/meta.yaml`.
  0.6.0 has the same API surface (`pass_command_as_script`, `array_execs`)
  and satisfies the submit plugin's `>=0.6.0,<1` requirement, so 2.8.0 pairs
  with it unchanged.
- **Upstream:**
  - Regression: <https://github.com/snakemake/snakemake-executor-plugin-slurm-jobstep/pull/46>
    (released as 0.6.1, 2026-05-27)
  - Fix: <https://github.com/snakemake/snakemake-executor-plugin-slurm-jobstep/pull/52>
    — one-line addition of `SLURM_CONF` to `keep_slurm`, open since 2026-07-09
  - Memory fallout: <https://github.com/snakemake/snakemake-executor-plugin-slurm-jobstep/issues/53>
- **Check:** `python -c "import snakemake_executor_plugin_slurm_jobstep as j, inspect; s = inspect.getsource(j.Executor.run_job); print('SLURM_CONF' in s or 'popen_env' not in s)"`
  — prints `True` once the allow-list keeps `SLURM_CONF` (or the sanitization
  is gone).
- **Remove when:** PR #52 (or equivalent) ships in a release, at which point
  the pin becomes `>=<that version>` or drops entirely.

## `delete_slurm_environment()` strips `SLURM_CONF`

- **Constraint:** when Snakemake is launched from inside a SLURM allocation,
  `snakemake-executor-plugin-slurm` calls `delete_slurm_environment()`
  (`utils.py`), which deletes *every* variable starting with `SLURM_` — with no
  exclusions. `SLURM_CONF` is one of them. On clusters where `SLURM_CONF` is the
  only route to `slurm.conf` (no `/etc/slurm/slurm.conf`; the path comes from an
  environment module, as on Bright Cluster Manager sites), the wipe propagates
  through `os.environ` into every submitted job, and the inner `srun` the
  jobstep executor uses falls through to a DNS SRV lookup and dies:

  ```
  srun: error: resolve_ctls_from_dns_srv: res_nsearch error: Unknown host
  srun: error: fetch_config: DNS SRV lookup failed
  srun: fatal: Could not establish a configuration source
  ```

  It hits the first rule scheduled, so it reads as a rule-specific bug when it
  is actually cluster-wide. Note every rule runs under an inner `srun`: the
  submit plugin hardcodes `--executor slurm-jobstep`, which builds
  `srun -n1 --cpu-bind=q ...`. There is no option to disable it.
- **Workaround:** run `epicc run --profile <slurm>` from a login/dev node rather
  than wrapping it in `sbatch` — the plugin warns against the nested case
  itself. Where an sbatch-wrapped launch is unavoidable, re-export the path in
  the profile's `precommand`, which runs inside each job downstream of the wipe:
  `precommand: 'export SLURM_CONF=/path/to/slurm.conf && ...'`. Both are
  documented in `profiles/slurm/config.yaml`.
- **Upstream:** <https://github.com/snakemake/snakemake-executor-plugin-slurm/issues/164>
  — reported 2024-11-04 and **closed without the fix landing**; the report names
  `SLURM_CONFIG`, which is not a real SLURM variable (the real one is
  `SLURM_CONF`), so the actual bug was likely never reproduced. Verified against
  plugin 2.7.0 and upstream `main` on 2026-08-19: the function still has no
  exclusion list.
- **Check:** `python -c "import snakemake_executor_plugin_slurm.utils as u, inspect; print('SLURM_CONF' in inspect.getsource(u.delete_slurm_environment))"`
  — prints `True` once the function grows an exclusion.
- **Remove when:** the plugin preserves `SLURM_CONF` (or offers an opt-out), at
  which point the `precommand` note in the SLURM profile can go.

## `--exclude=None assigned` in `snakemake-executor-plugin-slurm`

- **Constraint:** On a `NODE_FAIL`, the plugin runs
  `sacct -j <id> -n -X -o nodelist` and adds the result to the set of nodes it
  excludes from later submissions. When the job never reached a node, SLURM
  returns the placeholder `None assigned`, which the plugin stores verbatim and
  then interpolates unquoted into every subsequent `sbatch` as
  `--exclude=None assigned`. The space makes `assigned` a positional argument,
  so sbatch reports `Unable to open file assigned` and no further job is ever
  submitted. One flaky node ends the run.
- **Workaround:** none available — the tracking is unconditional and
  `--slurm-exclude-failed-nodes` only *seeds* the same set. The poisoned list is
  per-process, so re-running resumes the workflow.
- **Upstream:** introduced in
  snakemake/snakemake-executor-plugin-slurm#411; present through 2.8.0. The
  value needs the placeholder filtered out and `safe_quote()` applied, as
  `--qos` and `--reservation` already do.
- **Check:** `grep -A6 'NODE_FAIL' $(python -c 'import
  snakemake_executor_plugin_slurm as m; print(m.__file__)') | grep -q
  'None assigned'` — starts succeeding once upstream filters the placeholder.
- **Remove when:** a release filters it; then drop README known-issue 6.

