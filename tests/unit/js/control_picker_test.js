// Behavioral tests for the Control pick list in tools/epicc-builder.html (#91).
//
// Hovering (or arrowing to) a Control candidate highlights the row it refers
// to. It must never scroll the table: the list editor is anchored to the cell
// being edited, and scrolling that cell away closed the editor mid-pick, so
// only the first, already-visible candidate could ever be chosen.
//
// Run via tests/unit/test_builder_sample_id.py, or directly:
//   node tests/unit/js/control_picker_test.js

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

function fakeElement() {
  const classes = new Set();
  return { classList: { add: (c) => classes.add(c), remove: (c) => classes.delete(c),
                        contains: (c) => classes.has(c) } };
}

const DATA = [
  { Sample_ID: "WT_ChIP_H3K9me2_rep1", Assay: "ChIP", IP_target: "H3K9me2", Genome: "G" },
  { Sample_ID: "WT_ChIP_Input_rep1", Assay: "ChIP", IP_target: "Input", Genome: "G" },
  { Sample_ID: "WT_ChIP_Input_rep2", Assay: "ChIP", IP_target: "Input", Genome: "G" },
];
const rows = DATA.map((d) => {
  const el = fakeElement();
  return { getData: () => d, getElement: () => el };
});
const ctx = {
  scrolls: 0,
  table: {
    getData: () => DATA,
    getRows: () => rows,
    scrollToRow() { ctx.scrolls++; return Promise.resolve(); },
  },
};
vm.createContext(ctx);
vm.runInContext([
  "var _ctrlHighlighted = null;",
  ...["escapeHtml", "controlItemLabel", "highlightControlRow"].map(extract),
].join("\n"), ctx);

let passed = 0, failed = 0;
function check(desc, cond) {
  if (cond) { passed++; } else { failed++; console.log("FAIL: " + desc); }
}

// Move through every candidate, as hovering or arrowing down the list does.
for (const d of DATA.slice(1)) ctx.highlightControlRow(d.Sample_ID);
check("highlighting never scrolls the table", ctx.scrolls === 0);
check("only the last hovered row stays highlighted",
  !rows[1].getElement().classList.contains("row-ref-highlight") &&
  rows[2].getElement().classList.contains("row-ref-highlight"));
ctx.highlightControlRow("");
check("clearing removes the highlight", !rows[2].getElement().classList.contains("row-ref-highlight"));
ctx.highlightControlRow("not_a_sample");
check("unknown ID highlights nothing", rows.every((r) => !r.getElement().classList.contains("row-ref-highlight")));

const label = ctx.controlItemLabel("WT_ChIP_Input_rep2");
check("label names the row number, got " + label, /row 3\b/.test(label));
check("label carries assay and target", /ChIP/.test(label) && /Input/.test(label));
check("label escapes HTML", ctx.controlItemLabel("<x>") === "&lt;x&gt;");

console.log(passed + " passed, " + failed + " failed");
process.exit(failed ? 1 : 0);
