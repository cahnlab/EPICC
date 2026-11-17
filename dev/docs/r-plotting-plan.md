# Moving EPICC's plots to R, with one theme

Status: groundwork landed (shared theme module `workflow/scripts/plot_theme.R`;
the mC per-context metaplot uses it). Everything below "Phased plan" is still
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
  of each row across all columns, within each group. Rows whose summary is NaN
  must land where deeptools puts them (the end); check this against a real
  matrix.
- **`--outFileSortedRegions`.** The `most` heatmap's sorted BED feeds
  `sort_heatmap`, which reorders the mC matrix with `computeMatrixOperations
  sort -R` for `Heatmap_sorted__`. The R heatmap must write the same BED
  (deeptools' six columns plus its header line, rows in plotted order, groups
  in order).
- **Scales.** `computing_matrix_scales` writes one `--zMin`/`--zMax` and one
  `--yMin`/`--yMax` value per matrix column, in column order (`heatmaps_scales:
  type|sample`). With `default` the files are empty, and deeptools uses one
  scale for everything. R reads these files rather than recomputing them,
  until phase 5 moves the scale computation itself.
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
| `R_mc_profile_by_context.R` | `Profile_pergroup__*__mC__*.pdf` | ggplot2 | **Uses the theme now.** |
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
`line_width` 0.6, `dash_cycle`, `grid_minor` FALSE, `panel_spacing` 1.4,
`strip_fill` grey92, `strip_face` bold, `legend_position` right,
`legend_key_width` 2.2, `legend_max_rows` 25, `legend_row_height` 0.22,
`legend_char_width` 0.09, `palette` (tab20), `heatmap_colors` {default, ChIP,
TF, ATAC, RNA, sRNA → seismic; mC → Oranges}, `heatmap_missing_color` grey90,
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
only ggplot2 and jsonlite; the env change is part of phase 1 so the groundwork
PR leaves `epibutton.yaml` untouched and no env is rebuilt.

---

## 4. Phased plan (one PR per plot family)

### Phase 1: env + matrix reader
- Declare r-jsonlite, r-data.table, r-scales, bioconductor-complexheatmap in
  `epibutton.yaml`. Every rule in that env gets a one-time env rebuild.
- `workflow/scripts/matrix_io.R`: `read_matrix_header()`, `read_matrix()`
  (header + values, sample/group boundaries, bin layout, tick positions/labels
  from the header), `group_means()`, `sort_rows()`, `bin_rows()`. Unit tests on
  small synthetic matrices, including a stranded + `nostrand` merge and a
  matrix with `nan`s. Cross-check `group_means` against `plotProfile
  --outFileNameData` on the pombe fixture.
- Translate-then-deprecate `heatmaps_plot_params`/`profiles_plot_params`: a
  Python helper maps known flags to theme overrides, warns, and refuses other
  flags.

### Phase 2: profiles (`Profile__`, `Profile_pergroup__`)
- `R_profile.R`, reading the matrix (option A) and the `params_regions` /
  `params_profile` files. `Profile__`: one panel per column, a line per region
  group, y-limits per column. `Profile_pergroup__` non-mC: one panel per
  region group, a line per sample. `Profile_pergroup__` mC: fold
  `R_mc_profile_by_context.R` in, reading the matrix instead of the
  `plotProfile` dump.
- Drops two `plotProfile` calls per profile job. The `--plotHeight`
  workaround goes away; legend room comes from `epicc_legend_*`.

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

### Phase 5: scales into R (optional)
- Move `computing_matrix_scales`' awk/bash into R (same per-mark /
  per-sample logic and range padding), which drops its `plotProfile` and
  `dataRange` passes. The per-column scale files keep their format so the
  phase 2–3 scripts don't change.

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
| `mc_profile` (groundwork) | `Profile_pergroup__*__mC__*` | `R_mc_profile_by_context.R` at the #112 head | the same script on `plot_theme.R`, default theme + `presentation` | `mc_profile/compare.png`: the default theme is pixel-identical |
| `profile` (phase 2) | `Profile__`, `Profile_pergroup__` (non-mC and mC) | `plotProfile` / `--perGroup` | `R_profile.R` | `profile/compare.png`, one row per plot type |
| `heatmap` (phase 3) | `Heatmap__`, `Heatmap_sorted__` | `plotHeatmap` | `R_heatmap.R` | `heatmap/compare.png`; also diff the `_sorted_regions.bed` |
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
  columns of a row, NaN-aware, within each group. Check on a matrix with NaNs
  (mC) that the R order matches `--outFileSortedRegions` from deeptools before
  switching.
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
