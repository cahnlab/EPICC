// Behavioral tests for commented samples in tools/epicc-builder.html (#65).
//
// Runs the shipped importTSV / exportTSV against a stub table: a sheet with
// '#'-prefixed sample rows, free-text comments and blank lines must import with
// the commented rows flagged, and export back to the identical text.
//
// Run via tests/unit/test_builder_commented_rows.py, or directly:
//   node tests/unit/js/commented_rows_test.js [sheet.tsv]
// With a path, prints the builder's view of that sheet as JSON and exits
// (used by the pytest wrapper to compare against read_sample_sheet).

const fs = require("fs");
const path = require("path");
const vm = require("vm");

const BUILDER = process.env.EPICC_BUILDER_HTML ||
  path.join(__dirname, "..", "..", "..", "tools", "epicc-builder.html");
const src = fs.readFileSync(BUILDER, "utf8");

function extract(name) {
  const m = src.match(new RegExp("^function " + name + "\\([\\s\\S]*?^}", "m"));
  if (!m) throw new Error("could not extract " + name + " from " + BUILDER);
  return m[0];
}
function extractConst(name) {
  const m = src.match(new RegExp("^(?:const|var) " + name + " = [\\s\\S]*?;$", "m"));
  if (!m) throw new Error("could not extract " + name + " from " + BUILDER);
  return m[0].replace(/^const /, "var ");
}

function makeContext() {
  const ctx = {
    factors: [],
    alerts: [],
    downloads: [],
    tableRows: [],
    setTimeout: () => {},
    console,
  };
  ctx.alert = (m) => ctx.alerts.push(m);
  ctx.table = {
    clearSort() {},
    setData(rows) { ctx.tableRows = rows; },
    getData() { return ctx.tableRows; },
  };
  vm.createContext(ctx);
  vm.runInContext([
    extractConst("COLUMNS"),
    extractConst("COMMENT_PREFIX_RE"),
    "var sheetPreamble = [], sheetTrailer = [];",
    "function factorField(n) { return '_factor_' + n; }",
    "function suggestSampleId() { return ''; }",
    "function rebuildTableColumns() {}",
    "function renderFactorTiles() {}",
    "function refreshUI() {}",
    "function sheetFilename() { return 'sheet.tsv'; }",
    "function downloadBlob(content) { downloads.push(content); }",
    ...["isCommentLine", "parseCommentedRow", "importTSV", "parseLevels",
        "populateFactorsFromLevels", "buildLevelsFromFactors", "sanitizeCell",
        "buildExportMatrix", "commentedRowKeepsImportedLevels", "exportTSV"].map(extract),
  ].join("\n"), ctx);
  return ctx;
}

function builderView(text) {
  const ctx = makeContext();
  ctx.importTSV(text);
  return ctx;
}

if (process.argv[2]) {
  const ctx = builderView(fs.readFileSync(process.argv[2], "utf8"));
  console.log(JSON.stringify(ctx.tableRows.map((r) => ({
    Sample_ID: r.Sample_ID, commented: !!r._commented,
  }))));
  process.exit(0);
}

let passed = 0, failed = 0;
function check(desc, cond) {
  if (cond) { passed++; } else { failed++; console.log("FAIL: " + desc); }
}

const H = ["Sample_ID", "Assay", "Genome", "Levels", "Replicate_ID", "Read_files",
  "Read_layout", "IP_target", "Control", "Peak_type", "Comments"].join("\t");
const row = (...c) => c.join("\t");
const SHEET = [
  "# preamble note",
  "",
  H,
  row("WT_Input_rep1", "ChIP", "ColCEN", "genotype:WT,tissue:leaf", "rep1", "SRR1", "SE", "Input", "", "", ""),
  row("WT_K9_rep1", "ChIP", "ColCEN", "genotype:WT,tissue:leaf", "rep1", "SRR2", "SE", "H3K9me2", "WT_Input_rep1", "broad", "a # in a comment"),
  "",
  "# --- parked samples ---",
  "#" + row("WT_K9_rep2", "ChIP", "ColCEN", "genotype:WT,tissue:leaf", "rep2", "SRR3", "SE", "H3K9me2", "WT_Input_rep1", "broad", ""),
  "  ## " + row("old_K9_rep1", "ChIP", "ColCEN", "genotype:WT,tissue:leaf,time:0h", "rep1", "SRR4", "SE", "H3K9me2", "", "broad", "old design"),
  "#" + H,
  "# prose\twith a stray tab",
  row("WT_ATAC_rep1", "ATAC", "ColCEN", "genotype:WT,tissue:leaf", "rep1", "https://x.org/a.fq.gz#frag", "SE", "", "", "", ""),
  "# trailing note",
].join("\n") + "\n";

// --- import ---------------------------------------------------------------
let ctx = builderView(SHEET);
const ids = ctx.tableRows.map((r) => r.Sample_ID);
const commented = ctx.tableRows.map((r) => !!r._commented);
check("no import alerts: " + ctx.alerts.join("; "), ctx.alerts.length === 0);
check("five sample rows, got " + JSON.stringify(ids),
  JSON.stringify(ids) === JSON.stringify(["WT_Input_rep1", "WT_K9_rep1", "WT_K9_rep2", "old_K9_rep1", "WT_ATAC_rep1"]));
check("'#' stripped and rows flagged commented",
  JSON.stringify(commented) === JSON.stringify([false, false, true, true, false]));
check("factors come from active rows, not the parked 3-factor row",
  JSON.stringify(ctx.factors) === JSON.stringify(["genotype", "tissue"]));
check("'#' inside a field is kept", ctx.tableRows[4].Read_files === "https://x.org/a.fq.gz#frag");
check("preamble kept", JSON.stringify(ctx.sheetPreamble) === JSON.stringify(["# preamble note", ""]));
check("notes ride on the next row",
  ctx.tableRows[2]._notesBefore.length === 2 && ctx.tableRows[4]._notesBefore.length === 2);
check("trailing note kept", JSON.stringify(ctx.sheetTrailer) === JSON.stringify(["# trailing note"]));

// --- round trip -----------------------------------------------------------
ctx.exportTSV();
check("round trip is byte-identical", ctx.downloads[0] === SHEET);
if (ctx.downloads[0] !== SHEET) console.log(ctx.downloads[0]);

// --- toggling -------------------------------------------------------------
ctx = builderView(SHEET);
ctx.tableRows[0]._commented = true;       // comment out an active row
ctx.tableRows[2]._commented = false;      // uncomment a parked one
ctx.exportTSV();
const out = ctx.downloads[0].split("\n");
check("newly commented row gets a plain '#'", out.indexOf("#" + row("WT_Input_rep1", "ChIP", "ColCEN",
  "genotype:WT,tissue:leaf", "rep1", "SRR1", "SE", "Input", "", "", "")) !== -1);
check("uncommented row loses its '#'", out.indexOf(row("WT_K9_rep2", "ChIP", "ColCEN",
  "genotype:WT,tissue:leaf", "rep2", "SRR3", "SE", "H3K9me2", "WT_Input_rep1", "broad", "")) !== -1);

// Editing a parked row's factor cell rebuilds its Levels from the factor set.
ctx = builderView(SHEET);
ctx.tableRows[3]._factor_genotype = "mut";
ctx.exportTSV();
check("edited parked row's Levels rebuilt",
  ctx.downloads[0].indexOf("old_K9_rep1\tChIP\tColCEN\tgenotype:mut,tissue:leaf\t") !== -1);

// --- no header ------------------------------------------------------------
ctx = builderView("# only\n# comments\n");
check("all-comment file alerts instead of importing", ctx.alerts.length === 1 && ctx.tableRows.length === 0);

console.log(passed + " passed, " + failed + " failed");
process.exit(failed ? 1 : 0);
