"""The TMPDIR shell.prefix must not cost the workflow bash strict mode.

`set -euo pipefail` is Snakemake's own default for shell commands, but it
applies only while no shell.prefix is registered: `snakemake.shell.
_get_process_prefix` gates it on `_process_prefix is None`. Registering the
TMPDIR prefix therefore drops strict mode for *every* rule unless the prefix
re-asserts it, and a command failing anywhere but last in a rule body then
lets the job exit 0 with truncated output -- `cat missing | pigz > out` leaves
a 23-byte empty gzip that runs on to zero-coverage methylation calls.

These tests assert the *behaviour*, not the shape of the prefix string: the
separators and where the `set` lands are free to change as long as the rule
body still runs strict. `test_stripping_strict_mode_reintroduces_silent_
success` is the control that keeps the other assertions honest -- it shows
the same rule body exiting 0 once strict mode is removed, so a pass here
cannot come from strict mode arriving by some other route.

The prefix under test is read out of workflow/Snakefile via its AST, so it
cannot drift from what actually ships.
"""

import ast
import os
import re
import shutil
import subprocess
import textwrap

import pytest

_REPO_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", ".."))
_SNAKEFILE = os.path.join(_REPO_ROOT, "workflow", "Snakefile")

# Bash brace groups are doubled because Snakemake .format()s the shell string;
# `{output}` is the one real field. The trailing marker must never be reached:
# errexit has to stop the body at the failing command.
_PIPELINE_FAILURE = '''
rule all:
    output: "out.txt"
    shell:
        """
        {{
        cat /nonexistent/file | gzip > {output}
        printf "REACHED_END\\n"
        }} 2>&1 | tee -a log.txt
        """
'''

# Same, without a pipeline: covers errexit on its own rather than pipefail.
_PLAIN_FAILURE = '''
rule all:
    output: "out.txt"
    shell:
        """
        {{
        false
        printf "REACHED_END\\n" > {output}
        }} 2>&1 | tee -a log.txt
        """
'''

# Covers nounset.
_UNSET_VARIABLE = '''
rule all:
    output: "out.txt"
    shell:
        """
        {{
        printf "%s\\n" "$DEFINITELY_NOT_SET_BY_ANYTHING" > {output}
        printf "REACHED_END\\n"
        }} 2>&1 | tee -a log.txt
        """
'''

_HEALTHY = '''
rule all:
    output: "out.txt"
    shell:
        """
        {{
        printf "%s\\n" "$TMPDIR"
        test -d "$TMPDIR" && printf "TMPDIR_EXISTS\\n"
        }} > {output} 2>&1
        """
'''


def _snakefile_source():
    with open(_SNAKEFILE) as fh:
        return fh.read()


def _prefix_call_source():
    """The shell.prefix(...) call as written, located by balancing parens.

    The call's arguments hold no parentheses today. If that changes the slice
    stops parsing and the AST walk below fails loudly rather than silently
    testing the wrong string.
    """
    src = _snakefile_source()
    starts = [m.start() for m in re.finditer(r"\bshell\.prefix\s*\(", src)]
    assert len(starts) == 1, (
        f"expected exactly one shell.prefix(...) in {_SNAKEFILE}, "
        f"found {len(starts)}")
    start = starts[0]
    depth = 0
    for i in range(src.index("(", start), len(src)):
        if src[i] == "(":
            depth += 1
        elif src[i] == ")":
            depth -= 1
            if depth == 0:
                return src[start:i + 1]
    raise AssertionError("unbalanced shell.prefix(...) call")


def _shipped_prefix(tmpdir_base):
    """The prefix string as Snakemake will receive it.

    A Snakefile is not valid Python (rule/configfile/include are DSL), so only
    the extracted call is parsed, and no eval happens: the argument is a
    literal whose one interpolation is `_tmpdir_base`, substituted here.
    """
    tree = ast.parse(_prefix_call_source())
    call = next(
        node for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and node.func.attr == "prefix"
    )
    out = []

    def emit(node):
        if isinstance(node, ast.Constant):
            out.append(node.value)
        elif isinstance(node, ast.JoinedStr):
            for part in node.values:
                emit(part)
        elif isinstance(node, ast.FormattedValue):
            expr = ast.unparse(node.value)
            assert expr == "_tmpdir_base", (
                f"unexpected interpolation {expr!r} in the prefix; this test "
                "only knows how to substitute _tmpdir_base")
            out.append(str(tmpdir_base))
        else:  # pragma: no cover - guards against a non-literal argument
            raise AssertionError(f"prefix is not a literal: {ast.dump(node)}")

    emit(call.args[0])
    return "".join(out)


def _without_strict_mode(prefix):
    """The same prefix with the strict-mode assertion removed.

    Takes the following separator with it so the remainder is still valid
    bash whether the prefix chains with ';' or '&&'.
    """
    stripped, n = re.subn(r"set -euo pipefail\s*(?:;|&&)?\s*", "", prefix, count=1)
    assert n == 1, "no strict-mode assertion to strip"
    return stripped


class TestPrefixSource:
    def test_reasserts_strict_mode(self):
        assert "set -euo pipefail" in _prefix_call_source(), (
            "the prefix must re-assert strict mode; registering any prefix "
            "disables Snakemake's own")

    def test_strict_mode_is_not_disabled_afterwards(self):
        source = _prefix_call_source()
        if "set -euo pipefail" not in source:
            pytest.skip("no strict mode to re-disable; see "
                        "test_reasserts_strict_mode")
        tail = source.rsplit("set -euo pipefail", 1)[1]
        assert not re.search(r"set \+(?:[eux]|o\s+(?:errexit|nounset|pipefail))",
                             tail), (
            f"strict mode is turned back off later in the prefix: {tail!r}")


@pytest.mark.skipif(shutil.which("snakemake") is None,
                    reason="Snakemake not installed")
class TestPrefixBehaviour:
    def _run(self, tmp_path, rule_text, prefix=None):
        if prefix is None:
            prefix = _shipped_prefix(tmp_path / "results")
        snakefile = tmp_path / "Snakefile"
        # %r keeps the doubled braces intact: Snakemake .format()s the prefix
        # itself, exactly as it does the one in workflow/Snakefile.
        snakefile.write_text("shell.prefix(%r)\n" % prefix
                             + textwrap.dedent(rule_text))
        return subprocess.run(
            ["snakemake", "-j1", "--nolock", "-s", str(snakefile),
             "-d", str(tmp_path)],
            # cwd is the temp dir so the repo's profiles/default is not picked up.
            cwd=str(tmp_path),
            capture_output=True, text=True, timeout=300,
        )

    def _assert_body_stopped(self, tmp_path, result):
        assert result.returncode != 0, (
            "a command failing mid-rule must fail the job:\n" + result.stdout)
        assert not (tmp_path / "out.txt").exists(), \
            "Snakemake should have removed the incomplete output"
        log = tmp_path / "log.txt"
        if log.exists():
            assert "REACHED_END" not in log.read_text(), \
                "the body ran past the failing command"

    def test_pipeline_failure_fails_the_job(self, tmp_path):
        """pipefail: the failing command is the head of a pipe."""
        self._assert_body_stopped(
            tmp_path, self._run(tmp_path, _PIPELINE_FAILURE))

    def test_plain_failure_fails_the_job(self, tmp_path):
        """errexit on its own, with no pipeline involved."""
        self._assert_body_stopped(
            tmp_path, self._run(tmp_path, _PLAIN_FAILURE))

    def test_unset_variable_fails_the_job(self, tmp_path):
        """nounset."""
        self._assert_body_stopped(
            tmp_path, self._run(tmp_path, _UNSET_VARIABLE))

    def test_stripping_strict_mode_reintroduces_silent_success(self, tmp_path):
        """Control: without strict mode the same body exits 0 and ships junk.

        If this ever starts failing, the assertions above have stopped
        measuring the prefix and are passing for some other reason.
        """
        prefix = _without_strict_mode(_shipped_prefix(tmp_path / "results"))
        result = self._run(tmp_path, _PIPELINE_FAILURE, prefix=prefix)
        assert result.returncode == 0, (
            "expected the known silent-success failure mode without strict "
            "mode:\n" + result.stdout + result.stderr)
        assert (tmp_path / "out.txt").exists(), \
            "expected a truncated output file to survive"
        assert "REACHED_END" in (tmp_path / "log.txt").read_text(), \
            "expected the body to run past the failing command"

    def test_healthy_rule_runs_with_tmpdir_routed(self, tmp_path):
        result = self._run(tmp_path, _HEALTHY)
        assert result.returncode == 0, result.stdout + result.stderr
        written = (tmp_path / "out.txt").read_text()
        tmpdir = written.splitlines()[0].strip()
        assert "TMPDIR_EXISTS" in written, \
            f"TMPDIR was not created: {written!r}"
        assert tmpdir.startswith(str(tmp_path)), \
            f"TMPDIR should sit under the output dir: {tmpdir!r}"
        assert "/.tmp/" in tmpdir, f"TMPDIR not routed: {tmpdir!r}"
        # The doubled braces must have survived Snakemake's .format() and been
        # expanded by bash, not left as a literal.
        assert "{" not in tmpdir and "}" not in tmpdir, \
            f"brace escaping leaked into TMPDIR: {tmpdir!r}"
