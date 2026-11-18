# Moving EPICC's plots to R, with one theme

Status: groundwork landed (shared theme module `workflow/scripts/plot_theme.R`).
Phase 1 done (env, `matrix_io.R`, the per-mark scales computed in R). Phase 2
done (`Profile__` and `Profile_pergroup__` drawn in R, the profiles part of
the flag translation back, the phase-1 findings settled). Phases 3–7 are still
to do.

The goal is for every plot EPICC draws to take its look from one editable theme
and to be drawn in R. Today most of the combined-analysis plots come from
deeptools 4, whose defaults have shifted under the pipeline more than once (see
the deeptools notes in `CLAUDE.md`), and the existing R plots each set their
own fonts, palettes and sizes.

Decisions taken with the user before the groundwork:

- **Heatmaps: ComplexHeatmap**, not ggplot2. It provides one `Heatmap` per
  sample joined with `+`, each with its own colour function and legend (our
  `--zMin`/`--zMax` per column), `row_split` for region groups, `anno_lines` for
  the profile above each column, and rasterisation for large matrices. The
  theme reaches it through `epicc_gpar()`.
- **Theme config: named theme + YAML key overrides.** `plot_theme: default`
  picks a built-in theme. `plot_theme_overrides:` replaces individual keys.
  Unknown keys and wrong value types are errors.
- **`heatmaps_plot_params` / `profiles_plot_params`: translate, then
  deprecate.** The known deeptools flags (`--colorMap`, `--missingDataColor`,
  `--interpolationMethod`, `--plotType`) are mapped onto theme keys with a
  deprecation warning. Any other flag stops the run with a message naming the
  theme keys.

---

## 1. Inventory

### Plots deeptools draws

| Output | Rule | Tool | Inputs |
|---|---|---|---|
| `combined/plots/Heatmap__{mp}__{env}__{an}__{G}__{target}.pdf` + `combined/matrix/Heatmap__…_sorted_regions.bed` | `plotting_heatmap_on_targetfile` | `plotHeatmap` | `final_matrix_*.gz`, `params_regions_*` (`--regionsLabel`), `params_heatmap_*` (`--zMin/--zMax`), `params_profile_*` (`--yMin/--yMax` for the profile on top), display table |
| `combined/plots/Heatmap_sorted__{mp}__mC__…pdf` | `plotting_sorted_heatmap_on_targetfile` | `plotHeatmap --sortRegions keep` | `sorted_final_matrix_*__mC__*.gz` (from `sort_heatmap`), mC `params_heatmap`/`params_profile`, display table |
| `combined/plots/Profile__{mp}__{env}__…pdf` | `plotting_profile_on_targetfile` (plot1) | `plotProfile` | matrix, `params_regions`, `params_profile`, display table (labels + `--colors`) |
| `combined/plots/Profile_pergroup__{mp}__{env}__…pdf`, env ≠ mC | same rule (plot2) | `plotProfile --perGroup` | as above; one shared `--yMin`/`--yMax` (the min/max of the per-column ones), `--plotHeight` grown with the sample count |
| `combined/plots/PCA__{env}__{an}__{G}.pdf` | `plot_PCA_correlation` | `plotPCA --transpose` | `pca_matrix_*.npz` (`multiBamSummary` / `multiBigwigSummary`), `--labels`, `--colors`; an empty placeholder when PCA fails |
| `{ChIP,ATAC}/plots/Fingerprint__final__{sample}.png` | `make_fingerprint_plot` | `plotFingerprint` | IP and control BAMs; runs in the `epibutton_chip` env, which has no R |

`{mp}` is `regions`/`tss`/`tes`; `{env}` is `ChIP`, `ATAC`, `RNA`, `sRNA`,
`mC`, `TF` or the cross-type `all`/`most`.

What the R versions must reproduce:

- **Region groups.** `merging_matrix` relabels the plus/minus halves of a
  stranded target to the target name (one group) and the `nostrand` pass to
  `{target}_nostrand` (its own group). The R heatmap splits rows by group in
  header order; profiles draw one line (or one facet row) per group.
- **`--regionsLabel`.** Comes from `matrix_region_labels.py`, which reads the
  group labels and sizes from the matrix header (`matrix_header.py`). In R,
  read the same header JSON directly. Don't rebuild it from the BED: whether a
  `nostrand` pass ran, and `--skipZeros` (applied by `computeMatrix`, not the
  plot), decide both the group count and the region counts.
- **Sorting** (`heatmaps_sort_options`): `keep`, or `descend` by `mean`/`median`
  of each row across all columns, within each group. deeptools sorts with
  numpy's `argsort` (ascending, NaN last) and reverses it for `descend`, so an
  all-NaN row goes to the **top** of its group and tied rows come out in
  reverse order. `sort_rows()` reproduces this (checked in phase 1, deeptools
  3.5.6 and 4.0.0 alike).
- **`--outFileSortedRegions`.** The `most` heatmap's sorted BED feeds
  `sort_heatmap`, which reorders the mC matrix with `computeMatrixOperations
  sort -R` for `Heatmap_sorted__`. The R heatmap must write the same BED
  (deeptools' six columns plus its header line, rows in plotted order, groups
  in order).
- **Scales.** `computing_matrix_scales` writes one `--zMin`/`--zMax` and one
  `--yMin`/`--yMax` value per matrix column, in column order (`heatmaps_scales:
  type|sample`). With `default` the files are empty, and deeptools uses one
  scale for everything. Phase 1 ported the rule to R (`plot_scales.R`, see
  below), and the R plots compute their limits with it from the matrix they
  draw.
- **Colour maps.** `heatmaps_plot_params[env]` currently holds `--colorMap
  seismic` for everything but mC (`Oranges`, `--missingDataColor 0.9`). These
  are now the theme's `heatmap_colors` defaults (`seismic` → `.EPICC_SEISMIC`,
  `Oranges` → `.EPICC_ORANGES`), so the look carries over.
- **Interpolation.** deeptools draws the heatmap as an image with `bilinear`
  (most envs) or `nearest` (mC) interpolation. R pre-bins rows to the output
  resolution (mean of consecutive rows), then rasterises. That is close to
  `nearest` and slightly blurrier than `bilinear`; accept it as a visible
  change.
- **Tick labels.** `--refPointLabel start|end` for `tss`/`tes`;
  `--startLabel start --endLabel end` for `regions`; distances from the
  header's `upstream`/`downstream`/`body` and `bin size`.
- **Average type.** `Profile_pergroup__` mC and the scale computation use
  `profiles_scale` (mean/median); `Profile__` and the heatmap's top profile
  currently use the deeptools default (mean). Make the R versions use
  `profiles_scale` throughout, and note it in the changelog.
- **Labels and colours** come from the display table (`key, label, colour,
  track, series`) via `matrix_sample_labels.py`, in matrix column order.
  Nothing changes there.

### Plots already drawn in R

| Script | Output(s) | Engine | Notes for the theme |
|---|---|---|---|
| `R_profile.R` + `plot_profile.R` (phase 2) | `Profile__*.pdf`, `Profile_pergroup__*.pdf` (mC faceted on context) | ggplot2 | **Uses the theme.** Replaced `plotProfile` and `R_mc_profile_by_context.R`. |
| `R_browser_plot.R` | `Browser_full_chromosomes__*.pdf`, `single_browser__*.pdf` | Gviz | Fonts/sizes via Gviz `displayPars` (`fontsize`, `cex.*`), which can be mapped from the theme. Track colours already come from the pipeline palette. |
| `R_Upset_plot_{peaks,clusters,TSS}.R` + `upset_common.R` | `Upset_combined_*.pdf` | ComplexUpset (ggplot2) | `themes = upset_modify_themes(...)` takes ggplot themes; `colorlist` from all qualitative Brewer palettes should become the theme palette. Fixed category colours (`Distal_downstream`, …) become a theme key. |
| `R_mapping_stats.R` | `mapping_stats_*.pdf` | ggplot2 + cowplot | `brewer.pal(4, "Paired")` → a theme key for category fills. |
| `R_peak_stats.R` | `peak_stats_*.pdf` | ggplot2 | Same. |
| `R_sizes_stats.R` | `srna_sizes_stats_*.pdf` (full + zoom) | ggplot2 | Fixed fill per read category. |
| `R_plot_expression_level.R` | `RNA/plots/plot_expression__*.pdf` | ggplot2 | |
| `R_call_DEGs.R`, `R_call_srna_DEGs.R` | `MDS_*`, `BCV_*`, `Heatmap_{RNAseq,sRNA}_{cpm,zscore}__*.pdf`, DEG bar plots | base graphics (`plotMDS`, `heatmap.2`), ggplot2 | The MDS plots could become ggplot2 on `plotMDS(..., plot = FALSE)` coordinates; `heatmap.2` → ComplexHeatmap. Sample colours use all Brewer qualitative palettes, not the tab20-on-`levels_label` assignment, so a genotype gets a different colour here than in every other plot. |
| `R_GO_analysis.R` | `RNA/plots/*_treemap.pdf` | rrvgo `treemapPlot` (treemap/grid) | Only fonts can be themed; leave as is. |

Not R and not deeptools: `idr --plot` writes `{ChIP,ATAC}/plots/idr_*.png`.
Out of scope.

---

## 2. Theme design (implemented)

**Where it lives.** `workflow/scripts/plot_theme.R` holds the built-in themes
(`EPICC_THEMES`: `default`, `presentation`) as plain R lists, plus helpers:

| Helper | Use |
|---|---|
| `epicc_theme(theme_file)` | Load the theme a script should use: the JSON Snakemake wrote, or the default when `NULL`/`""`. |
| `epicc_ggtheme(th)` | ggplot2 `theme_bw` + strips, spacing, legend, minor grid. |
| `epicc_gpar(th, rel, ...)` | `grid::gpar` at the theme's font, for ComplexHeatmap and Gviz. |
| `epicc_dashes(th, colours)` | Dash pattern per series; only series that share a colour are dashed. |
| `epicc_legend_ncol/_width/_height` | Legend wrapping and the room it needs, for PDF sizing. |
| `epicc_heatmap_colors(th, type)` | Colour ramp per data type with fallbacks (`mCHG` → `mC`, `ChIP_broad` → `ChIP`, unknown → `default`). |
| `epicc_save(plot, file, w, h)` | One place to pick the device later (e.g. `cairo_pdf` for font embedding). |

**Theme keys** (default theme): `base_size` 11, `base_family` "",
`line_width` 0.6, `profile_type` lines (lines/fill/se/std), `dash_cycle`, `grid_minor` FALSE, `panel_spacing` 1.4,
`strip_fill` grey92, `strip_face` bold, `legend_position` right,
`legend_key_width` 2.2, `legend_max_rows` 25, `legend_row_height` 0.22,
`legend_char_width` 0.09, `palette` (tab20), `heatmap_colors` {default, ChIP,
TF, ATAC, RNA, sRNA → seismic; mC → Oranges; an entry may also be one matplotlib map name}, `heatmap_missing_color` grey90,
`panel_width`/`panel_height` 3.2 in, `heatmap_column_width` 1.6 in,
`heatmap_height` 9 in, `profile_strip_height` 1.2 in. The last three are for
phase 3 and unused so far.

**Sample colours stay out of the theme.** They arrive per plot in the display
table (tab20 keyed on `levels_label`, `assign_colors` in
`combined_analysis.smk`), so a genotype keeps its colour across every plot.
`palette` is the theme's fallback for plots without a display table. A unit
test keeps it identical to `_PALETTES["tab20"]`.

**User selection.** In the options YAML:

```yaml
plot_theme: "default"          # or "presentation"
plot_theme_overrides:          # optional; keys of the chosen theme
  base_size: 12
  heatmap_colors:
    RNA: ["white", "darkgreen"]   # heatmap_colors may add new data types
```

`plot_theme_spec()` (`config_validation.py`) checks the shape at startup and
serialises it to stable JSON (`PLOT_THEME_JSON`, defined in the Snakefile).
Plot rules take it as `params.theme`, so a theme change reruns plots only and
never data steps. The rule writes it to `$TMPDIR/plot_theme.json` and passes
the path as the script's last argument. `plot_theme.R` owns the key list and
refuses unknown keys and type changes when the plot is drawn. A startup-time
key check would mean running R from the Snakefile, so it is not done.

**How scripts load it**: the same self-locating `source()` the UpSet scripts
use for `upset_common.R`:

```r
script_dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
source(file.path(script_dir, "plot_theme.R"))
th <- epicc_theme(theme_file)
```

---

## 3. Data path

**Option A (chosen): read the `computeMatrix` `.gz` directly.** Line 1 is `@`
followed by the JSON header (group/sample labels and boundaries, bin size,
upstream/downstream/body, reference point). The rest is a TSV: six BED columns,
then one value per bin per sample, `nan` for missing. In R:
`jsonlite::fromJSON(sub("^@", "", readLines(gzcon(file(m)), n = 1)))` and
`data.table::fread(cmd = "gzip -dc m", skip = 1, header = FALSE, na.strings = "nan")`.

Benchmark (synthetic, ColCEN-TE scale: 43,000 regions × 12 tracks × 180 bins =
93 M values, 302 MB gzipped, epibutton env, R 4.5): `fread` took 5.4 s wall, the
numeric matrix is 709 MB, and the peak RSS of the whole script was 2.2 GB. Row
means for sorting plus column means for the profile took 0.7 s. deeptools also
loads the whole matrix into memory, so peak memory is in the same range as
today's `plotHeatmap`. Scaling: about 60 MB of R memory per 1 M values, so a
30-track `most` matrix on TEs (~230 M values) needs about 5–6 GB. Mitigations,
if needed: read one sample block at a time with `fread(select = ...)`, which
re-parses the file per block but holds one block in memory; and pre-bin rows to
the output resolution (≤ ~2,000 rows per group) straight after sorting, so
ComplexHeatmap never sees the full matrix. Give the heatmap rule `mem_mb`
proportional to the matrix size.

**Option B (rejected as the main path): deeptools text dumps.**
`plotProfile --outFileNameData` gives only the per-group averages, which is
enough for profiles but not for heatmaps. `plotHeatmap --outFileNameMatrix`
gives the sorted matrix, but only by running `plotHeatmap`, the very step being
replaced. `computeMatrixOperations` has no text dump. So B would keep two
deeptools invocations per plot and still need A for heatmaps. The current mC
profile uses B; phase 2 moves it to A so every profile reads the matrix the
same way.

**PCA and fingerprint keep deeptools for the maths.** `plotPCA
--outFileNameData` writes the PC coordinates and eigenvalues (after deeptools'
log2, top-variance row selection and `--transpose`). R then draws the scatter
and scree panels. Reading `.npz` in R would need reticulate or RcppCNPy, and
redoing the PCA in R would change the numbers. Fingerprint: `plotFingerprint
--outRawCounts` (or `--outQualityMetrics`) in the chip env, then a separate
epibutton-env rule draws the cumulative curves in R. That also makes the plot
PDF-capable.

**Packages.** Already declared in `epibutton.yaml`: r-ggplot2, r-cowplot,
r-rcolorbrewer, r-dplyr, r-tidyr, r-purrr, r-complexupset. Present in the
resolved env only as transitive dependencies, so they should be declared now
that we depend on them directly: **r-jsonlite** (theme, matrix header),
**r-data.table** (matrix reading), r-scales. To add for phase 3:
**bioconductor-complexheatmap** (pulls in r-circlize, r-getoptlong,
r-globaloptions, r-clue). Optional: r-patchwork (already transitive) for
composing ggplot panels (PCA scatter + scree). The theme module itself needs
only ggplot2 and jsonlite; phase 1 put the new packages in a separate env,
`epibutton_plots.yaml`, so `epibutton.yaml` stays untouched and nothing that
shares it is rerun.

---

## 4. Phased plan (one PR per plot family)

### Phase 1: env + matrix reader (done)
- A new env, `workflow/envs/epibutton_plots.yaml` (python, gzip, R with
  data.table, jsonlite, scales, ggplot2 and ComplexHeatmap), serves the R
  plotting and plot-scale rules. In phase 1 that is only
  `computing_matrix_scales`, which calls nothing but `python3`, `Rscript` and
  `gzip`. `epibutton.yaml` is unchanged, so the 38 rules that share it (genome
  prep, downloads, trimming, combined analysis) keep their software env and
  are not rerun. A dry-run against the existing pombe run gave the same rerun
  reason for every job under #113's workflow and this one, except
  `computing_matrix_scales`. That run is stale and reruns end to end either
  way, so the env claim rests on `epibutton.yaml` being byte-identical. The
  env solves and builds (R 4.5, data.table 1.18, ComplexHeatmap 2.26). Each
  later phase moves its plot rule to this env when R starts drawing that
  plot; the plot rules that still call deeptools stay where they are.
- `workflow/scripts/matrix_io.R`: `read_matrix_header()`, `read_matrix()`
  (regions, values with `nan` → `NA`, sample and group ranges, per-sample bin
  layout), `sample_columns()`, `row_summary()`/`col_summary()` (NaN-aware like
  numpy's `nan*`), `group_means()`, `sort_rows()`, `bin_rows()` and
  `matrix_ticks()` (deeptools' `getProfileTicks`: the upstream tick reads
  `-2.0`, the downstream one `2.0Kb`).
- `tests/unit/test_matrix_io.py` writes small matrices in deeptools' format (a
  stranded and a `nostrand` group, NaNs, an all-NaN row, tied rows, an empty
  group). Where deeptools is on PATH it also compares `group_means` with
  `plotProfile --outFileNameData` and `sort_rows` with `plotHeatmap
  --outFileSortedRegions`. These pass against deeptools 4.0.0, and the same
  comparison on a matrix built by `computeMatrix` from bigwigs with gaps
  matched to 4e-16. The comparison runs on generated matrices, not the pombe
  run, so it needs no completed pipeline run. Under `--sortUsing max`/`min`,
  rows with exactly equal summaries can come out in a different order, because
  numpy's float sort is not stable. EPICC sorts only by mean or median.
- **Flag translation: deferred to phase 2 (back in phase 2 for profiles).** `legacy_plot_params()`
  translates `heatmaps_plot_params`/`profiles_plot_params` onto theme keys:
  - `--colorMap NAME` → `heatmap_colors.<env>: NAME`;
  - `--colorList` → a list;
  - `--missingDataColor` → `heatmap_missing_color`;
  - `--interpolationMethod` is dropped;
  - `--plotType` → `profile_type`.

  Each translated setting gets a deprecation warning, and any other flag stops
  the run. It is gated per family by a Snakefile `R_DRAWN_PLOT_FAMILIES`.
  Review of #114 noted that it would never run in phase 1, so it was taken
  out. The code and its tests are parked on branch
  `claude/r-plot-legacy-flags` (commit `fc1130c`, in `config_validation.py`,
  `workflow/Snakefile` and `tests/unit/test_plot_theme.py`). Phase 2 brings it
  back with `R_DRAWN_PLOT_FAMILIES = ("profiles",)`, and phase 3 adds
  `"heatmaps"`. Once a family is listed, drop its key from
  `config/epicc-options.yaml` and the builder. Phase 2 did this for
  `profiles_plot_params` (restored from `fc1130c`).
- A `heatmap_colors` entry may now be one matplotlib colour-map name
  (`viridis`, `RdBu_r`), resolved through `workflow/scripts/plot_colormaps.json`
  (all 91 matplotlib maps, 17 stops each; regenerate with
  `dev/gen_plot_colormaps.py`). `plot_theme.R` checks names against it, and so
  will the translation when it returns.
- Open for phase 3: `heatmaps_plot_params` has `all`/`most` entries, which
  translate to `heatmap_colors.all`/`.most`. Decide whether the cross-type
  heatmaps use those or colour each column by its data type.
- **Per-type and per-mark scaling: `workflow/scripts/plot_scales.R`** is a
  port of `computing_matrix_scales` (bash, `dataRange`, `plotProfile
  --outFileNameData`). It was made a phase 1 acceptance requirement so that
  `heatmaps_scales` survives the switch. It provides:
  - `matrix_scales(m, mode, average, mc_zmax)`: one row per matrix column, in
    column order;
  - `scale_for_column()`: the lookup every R plot must use, by its own column
    key;
  - `shared_ylim()`: the `Profile_pergroup__` collapse to one min/max;
  - `default_scales()`: the single scale deeptools falls back to;
  - `scale_flags()`: the rule's exact flag strings.

  What it reproduces, all read off the rule rather than assumed:
  - **zMin/zMax** are `dataRange`'s **10th/90th percentiles**, not min/max.
    numpy propagates NaN, so a column holding any NaN has no range. Every mC
    column has NaNs, so mC heatmaps always land on `heat_mcg`/`heat_mchg`/
    `heat_mchh`. Any other mark whose range is 0..0 gets 0..0.005.
  - **`type`**: one range per mark. The mark is the last `__` field of the
    key, matched exactly. **`sample`**: one range per column.
  - **y**: the column's (or mark's) average profile over all region groups,
    padded by a tenth of the range. When the range is 0 the pad is |max|/10,
    then 0.01. So an all-zero profile gets −0.01..0.01, and only a profile
    with no value at all gets 0..0.01.
  - y is computed from `group_means()`, the same averaging the R profiles will
    draw with. That function is checked against `plotProfile`, so the bounds
    and the curves cannot disagree.

  `tests/unit/test_plot_scales.py` covers:
  - each mode;
  - interleaved stranded columns;
  - H3K9me next to H3K9me2;
  - marks two orders of magnitude apart;
  - an all-zero mC context, an all-zero other mark and an all-NaN column;
  - lookups by key;
  - with deeptools on PATH, the bash rule itself (`dev/scales_rule.py` runs its
    shell block from the `.smk`) against R in `type`/`sample` × mean/median.

  The `scales` comparison family draws a `type`-mode `most` heatmap with the
  rule's limits and with R's. The two are pixel-identical.

  **The rule now runs it.** At the user's call, all scaling moved to R in
  phase 1: `computing_matrix_scales` runs `R_matrix_scales.R`, which writes
  the same flag files, and no longer calls `dataRange` or `plotProfile`.
  Rule connections changed in two ways:
  - its `temp_values`, `temp_profile` and `temp_profile_values` outputs are
    gone (nothing else read them);
  - its code changed, so Snakemake reruns it, and the plots that read its
    files, once. No matrix is rebuilt.

  The old bash is frozen in
  `dev/reference/computing_matrix_scales_deeptools.sh` and run beside the new
  rule by `test_plot_scales.py` and the `scales` comparison family. The
  comparison covered three fixtures plus two matrices `computeMatrix` built
  from bigwigs (one with NaNs, one `--missingDataAsZero --skipZeros`
  reference-point), in every mode with mean and median. Results:
  - **z:** every value is the same double. Where the text differs, it is the
    spelling: `dataRange` prints `0.0` where R prints `0`, and Python may use
    16 digits where R uses 17.
  - **y:** text-identical everywhere. There was no deeptools-average
    difference to explain, because `group_means()` equals `plotProfile`'s
    averages (checked in `test_matrix_io.py`).
  - **Rendered check:** a heatmap drawn with either set of limits is
    pixel-identical.

  Phase 2's R profile takes its limits from the same functions, either by
  reading these files or by calling `matrix_scales()` directly.
- **Behaviours of today's rule found while porting**, settled in phase 2:
  - `Profile__` is drawn with deeptools' default mean even when
    `profiles_scale: median`, while its y-limits come from the median. Its
    curves can therefore leave their axis. **Fixed:** the R profile uses
    `profiles_scale` for both (`test_profile.py::test_lines_are_profiles_scale`).
  - With `heatmaps_scales: default`, the per-group collapse turns the empty
    scale file into `--yMin 99999 --yMax -99999`, and the old
    `Profile_pergroup__` came out blank (`most_default_mean` comparison).
    **Fixed:** R uses `default_scales()`.
  - Because of the NaN propagation above, mC heatmap scales come from config
    (`heat_m*`), never from the data. **Kept** for the heatmaps (phase 3
    decides); it does not touch the profiles, whose y range comes from the
    NaN-aware averages.
  - In `type` mode, the awk search for a mark's range starts from the
    sentinels 999999 and −999999. A mark whose 10th percentile is ≥ 999999 (or
    whose 90th is ≤ −999999) therefore gets 0 on that side. Large raw counts
    can hit this. **Fixed** in `plot_scales.R`: the mark keeps its range. It
    changes the heatmap `--zMin` of such a mark (the only output affected);
    the scale tests apply this one fix to the old rule's output before
    comparing (`old_rule_fixed`).
  - y bounds are written as awk printed them: integral values in full
    (`900000`), anything else with `%.6g`. **Kept** in the flag files, which
    only the deeptools heatmaps still read; the R profiles take the doubles
    straight from `matrix_scales()`, so the rounding no longer reaches a plot
    drawn in R. The flag files go when phase 3 moves the heatmaps.
  - An all-zero profile gets −0.01..0.01, not 0..0.01. **Kept:** a flat line
    at zero drawn mid-panel is right for signed tracks (`bamCompare` log2
    ratios), and 0..0.01 would draw it on the panel's bottom edge.

### Phase 2: profiles (`Profile__`, `Profile_pergroup__`) + scales (done)

What landed:
- `workflow/scripts/plot_profile.R` holds the reusable parts:
  `profile_data()` (the `group_means()` curves plus an se/std band: numpy's
  population `nanstd` per bin, divided by √n for se), `profile_limits()` (y
  range per column: `matrix_scales()` for `type`/`sample`, `default_scales()`
  for `default`), `profile_ggplot()` (lines, fill or a band per
  `profile_type`; per-panel limits drawn as invisible points with no
  expansion, so a panel spans exactly its column's range) and the three plots,
  `plot_profile_per_column()` (`Profile__`), `plot_profile_per_group()`
  (`Profile_pergroup__`) and `plot_profile_by_context()` (mC). Phase 3's
  heatmap top profile should call `profile_data()` and `profile_limits()`
  rather than recompute them.
- `R_profile.R` is the CLI. It reads the matrix and the display table
  (`matrix_sample_labels.py`), works out the limits itself (it does not read
  `params_profile`, so the profile rule no longer depends on
  `computing_matrix_scales`), and writes both PDFs.
  `plotting_profile_on_targetfile` runs it in `CONDA_ENV_PLOTS`.
  `R_mc_profile_by_context.R` and its `plotProfile --outFileNameData` dump are
  gone.
- Visible changes, all in the comparison images: `Profile__` draws each panel
  in its column's display colour with region groups told apart by dash
  pattern (deeptools coloured the groups with the first sample colours of the
  `--colors` list, an accident of reusing it); 4 panels per row instead of 8;
  a y-axis label (`mean signal`, `median methylation`); the mC context panels
  share one y range per context across region groups (previously autoscaled
  per panel; still autoscaled under `default`).
- `profiles_plot_params` is translated (`--plotType` → `profile_type`, with
  a deprecation warning; anything else stops the run) and is gone from the
  shipped options file, the README and the builder.

As planned:
- `R_profile.R`, reading the matrix (option A) and the `params_regions` /
  `params_profile` files. `Profile__`: one panel per column, a line per region
  group, y-limits per column. `Profile_pergroup__` non-mC: one panel per
  region group, a line per sample. `Profile_pergroup__` mC: fold
  `R_mc_profile_by_context.R` in, reading the matrix instead of the
  `plotProfile` dump.
- Drops two `plotProfile` calls per profile job. The `--plotHeight`
  workaround goes away; legend room comes from `epicc_legend_*`.
- Limits come from `plot_scales.R` on the matrix being drawn. The y bounds
  must come from the same `group_means()` call that produces the curves, so
  the plot cannot run past its own limits.
- Scales are already in R (phase 1). Each R profile must take its y-limits
  from `matrix_scales()` with the same `profiles_scale` it draws with.
- Tests must check that the R profile gives matrix column *i* the y-limits of
  column *i* (look up by key with `scale_for_column`, never by position within
  a mark), on an interleaved stranded matrix. Phase 3 adds the same check for
  zMin/zMax in the heatmap.

### Phase 3: heatmaps (`Heatmap__`, `Heatmap_sorted__`)
- `R_heatmap.R` on ComplexHeatmap: one `Heatmap` per column with
  `circlize::colorRamp2` over its `--zMin`/`--zMax` and the theme's ramp for
  the column's data type, `row_split` by group, `anno_lines` top profile with
  per-column y-limits, `use_raster = TRUE`, rows pre-binned.
- Writes `_sorted_regions.bed` in deeptools' format for `sort_heatmap`. Keep
  `sort_heatmap` (deeptools) as is at first. Once R sorts, it could reorder
  the mC matrix itself and drop `sort_heatmap`.

### Phase 4: PCA
- `plotPCA --outFileNameData` in the existing rule, then `R_pca.R` (ggplot2
  scatter + scree, patchwork), colours/labels from the display table. Keep the
  empty-placeholder behaviour when PCA fails.

### Phase 5: scales into R (done in phase 1)
- `computing_matrix_scales` runs `R_matrix_scales.R`. The aim is to take
  plotting off deeptools, not to remove it. deeptools stays for
  `computeMatrix`, `computeMatrixOperations` (rbind, sort, relabel),
  `bamCoverage`/`bamCompare` and the like.

### Phase 6: existing R plots adopt the theme
- UpSet, mapping/peak/size stats, expression, MDS (ggplot2 from `plotMDS(...,
  plot = FALSE)`), DEG heatmaps (ComplexHeatmap), browser (Gviz `displayPars`
  from the theme). Give the DEG/MDS plots the tab20-on-`levels_label` colours.
  Each script gets the theme JSON as a new last argument.

### Phase 7: fingerprint
- `plotFingerprint --outRawCounts` in the chip-env rule; a new epibutton-env
  rule draws the PDF in R. The `.png` output name changes to `.pdf`, so update
  `define_*` outputs and docs.

### Before/after comparison for every conversion

Each conversion PR comes with a side-by-side render of the old plot and the new
one, drawn from the same input and shown to the user before they sign off on
it. The renderer is `dev/plot_compare.py`:

```bash
python dev/plot_compare.py --out dev/plot_compare_out --old-ref <PR base branch> [--family F]
```

For each family it writes `old.pdf`, `new.pdf` (plus any variants such as
another theme), `compare.png` (the renders in one labelled row) and
`report.txt` (whether each render is pixel-identical to the old one) under
`dev/plot_compare_out/<family>/`, which is gitignored. The old render uses
`workflow/scripts` as of `--old-ref`, extracted with `git archive`. Each phase
adds its family to `FAMILIES`. Where the old plot came from deeptools, the
family runs the old deeptools command from the base branch's rule on a
fixture matrix. Small fixtures are generated in the script; real matrices
(pombe, ColCEN) are used when they are already on disk.

| Family | Plots compared | Old | New | Comparison images |
|---|---|---|---|---|
| `mc_profile` (groundwork and phase 1; removed in phase 2) | `Profile_pergroup__*__mC__*` | `R_mc_profile_by_context.R` at the #112 head (phase 1: at the #113 head) | the same script on `plot_theme.R`, default theme + `presentation` | the default theme was pixel-identical both times. Superseded by `profiles`. |
| `profiles` (phase 2) | `Profile__`, `Profile_pergroup__` (non-mC and mC) | `plotting_profile_on_targetfile`'s own shell block at `--old-ref` (`plotProfile`, `--perGroup`, the mC dump script), with that ref's scales rule | the same rule block now (`R_profile.R`) | `profiles/<variant>/{Profile,Profile_pergroup}/compare.png` for `most` (stranded + nostrand groups) under type/mean, sample/median, default/mean; mC under type/mean and default/median; a one-group reference-point ChIP matrix (tss) under type/mean; plus `profile_type: se` and the `presentation` theme. Needs deeptools on PATH |
| `heatmap` (phase 3) | `Heatmap__`, `Heatmap_sorted__` | `plotHeatmap` | `R_heatmap.R` | `heatmap/compare.png`; also diff the `_sorted_regions.bed` |
| `scales` (phase 1) | `heatmaps_scales: type` limits on a `most` heatmap (H3K9me2 0–60 next to H3K4me3 0–0.5, H3K9me, interleaved RNA plus/minus) | `computing_matrix_scales` as bash + `dataRange`/`plotProfile` (`dev/reference/`), drawn by `plotHeatmap` | the rule as it stands (`R_matrix_scales.R`), drawn by the same `plotHeatmap` | `scales/compare.png` + `scales/flags.txt` (both flag sets, largest difference per flag): pixel-identical. Needs deeptools on PATH |
| `pca` (phase 4) | `PCA__` | `plotPCA` | `R_pca.R` | `pca/compare.png` |
| `stats` / `upset` / `deg` / `browser` (phase 6) | existing R plots | current script | themed script | one `compare.png` per family |
| `fingerprint` (phase 7) | `Fingerprint__` | `plotFingerprint` PNG | R PDF | `fingerprint/compare.png` |

---

## 5. Risks and behaviours to match

- **Output names.** `tests/integration/test_pombe_postrun.py` checks
  `Heatmap__{tss,tes,regions}__most__{RUN}__all_genes.pdf`,
  `Profile__…__most__…pdf` and `Browser_full_chromosomes__all__{RUN}.pdf`, and
  `define_final_combined_output` lists every PDF. Keep each name. The only
  planned rename is the fingerprint `.png` → `.pdf` (phase 7).
- **Region groups.** Stranded targets give `target` + `target_nostrand`.
  `--skipZeros` drops regions before the matrix is written (non-mC envs only:
  `base_mc` doesn't set it), so the mC and `most` matrices of one target can
  have different region sets. `sort_heatmap` relies on `computeMatrixOperations
  sort -R` with the `most` order; verify what it does with mC regions missing
  from that BED before R takes over the sorting.
- **`--regionsLabel` vs header labels.** Today `merging_matrix` sets the group
  labels and `matrix_region_labels.py` re-emits them with counts. R should show
  the same text (`label (n)`, if that is what deeptools draws — check a current
  PDF).
- **Sort semantics.** deeptools `descend` + `mean` uses the mean over all
  columns of a row, NaN-aware, within each group. Phase 1 matched
  `--outFileSortedRegions` on matrices with NaNs and an all-NaN row (it goes
  first under `descend`). Recheck on a real mC matrix in phase 3.
- **Scale files with `heatmaps_scales: default`** are empty. R must then
  compute one shared scale, as deeptools does.
- **Colour fidelity.** `seismic`/`Oranges` are copied as hex stops. Spacing
  between stops is linear in both. Check a side-by-side PDF.
- **Memory.** Large `most` matrices on TEs (see section 3). Size `mem_mb` from
  the matrix's value count, and fall back to blockwise reading if a real run
  needs it.
- **Rerun triggers.** Theme JSON and display tables are params of plot rules
  only. `test_cosmetic_plot_labels.py` already guards that data rules take no
  label; extend it to the theme.
- **deeptools stays** for `computeMatrix`, `computeMatrixOperations`,
  `multiBam/BigwigSummary`, `bamCoverage`/`bamCompare` and (for the maths)
  `plotPCA`/`plotFingerprint`, so the deeptools 4 pins and notes in `CLAUDE.md`
  remain relevant.
- **The builder** (`tools/epicc-builder.html`) doesn't know `plot_theme` yet.
  Exported options files omit it and get the default, which is harmless. Add
  it when the theme has user-visible choices worth a UI.
