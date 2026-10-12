// Behavioral tests for the builder's Levels pair check (levelsPairIssue) in
// tools/epicc-builder.html. The function is extracted from the shipped HTML.
//
// Run via tests/unit/test_builder_levels.py, or directly:
//   node tests/unit/js/levels_validation_test.js
// Exits non-zero on failure.

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

const ctx = {};
vm.createContext(ctx);
vm.runInContext(extract("levelsPairIssue"), ctx);
const levelsPairIssue = ctx.levelsPairIssue;

let passed = 0, failed = 0;
function check(desc, cond) {
  if (cond) { passed++; } else { failed++; console.log("FAIL: " + desc); }
}

const EQ = "factor names and levels must not contain '='";
const EMPTY = "empty factor name or level in Levels";

check("plain pair passes", levelsPairIssue("genotype", "WT") === "");
check("dash and dot pass", levelsPairIssue("genotype", "WT-1.2") === "");
check("'=' in level rejected", levelsPairIssue("genotype", "WT=1") === EQ);
check("'=' in factor rejected", levelsPairIssue("geno=type", "WT") === EQ);
check("empty level", levelsPairIssue("genotype", "") === EMPTY);
check("blank factor", levelsPairIssue("  ", "WT") === EMPTY);
check("empty wins over '='", levelsPairIssue("=", "") === EMPTY);

console.log(passed + " passed, " + failed + " failed");
process.exit(failed ? 1 : 0);
