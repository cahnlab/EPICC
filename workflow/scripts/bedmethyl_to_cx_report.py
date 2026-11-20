#!/usr/bin/env python3
"""Convert modkit bedMethyl format to Bismark CX_report format.

Determines methylation context (CG/CHG/CHH) from the reference genome.

Only 5mC calls (column 4 == 'm') are retained. Other modification types
(e.g. 5hmC 'h') are skipped. If the input was produced with modkit's
--combine-mods flag, all modifications are already merged into a single
row per position and column 4 is 'C' — these rows are kept as-is.

When the bedMethyl has combined-strand data (strand column == '.'), the
strand is inferred from the reference: if the base at the position is C
on the forward strand, the call is assigned to '+'; if it is G (i.e. a
C on the reverse strand), the call is assigned to '-'. Note that for CpG
sites produced with --combine-strands, each position appears only once
(on the reference C strand), so the CX_report will have half the CpG
entries compared to Bismark output.

bedMethyl columns (modkit output):
    0: chrom
    1: start (0-based)
    2: end
    3: mod_code ('m' for 5mC, 'h' for 5hmC, 'C' for combined)
    4: score
    5: strand (. = combined, + or -)
    6-8: thickStart, thickEnd, itemRgb
    9: coverage
    10: percent methylated
    11: N_mod (modified count)
    12: N_canonical (unmodified count)

CX_report format (Bismark):
    0: chrom
    1: position (1-based)
    2: strand (+/-)
    3: count methylated
    4: count unmethylated
    5: context (CG/CHG/CHH)
    6: trinucleotide sequence

Usage:
    python bedmethyl_to_cx_report.py input.bed.gz reference.fa output.CX_report.txt.gz [required_contexts sample_name]

    required_contexts ("CG|CHG|CHH", from methylation_contexts): each must have
    at least one covered call, or the script exits 1 naming the fix.
"""

import sys
import gzip
import pysam
from contextlib import ExitStack


def get_context_and_trinuc_from_seq(chrom_seq, chrom_len, pos, strand):
    """Determine methylation context and trinucleotide from cached chromosome sequence.

    Args:
        chrom_seq: uppercase chromosome sequence string
        chrom_len: length of chromosome
        pos: 0-based position of the C
        strand: '+' or '-'

    Returns:
        (context, trinucleotide) tuple
    """
    if strand == '+':
        # Forward strand: look at C and next 2 bases
        end = min(pos + 3, chrom_len)
        seq = chrom_seq[pos:end]
        if len(seq) < 2:
            return ("CHH", seq + "N" * (3 - len(seq)))

        trinuc = seq[:3] if len(seq) >= 3 else seq + "N" * (3 - len(seq))

        # Determine context
        if len(seq) >= 2 and seq[1] == 'G':
            context = "CG"
        elif len(seq) >= 3 and seq[2] == 'G':
            context = "CHG"
        else:
            context = "CHH"
    else:
        # Reverse strand: C is actually a G on forward strand
        # Look at 2 bases before and the G position
        start = max(0, pos - 2)
        seq = chrom_seq[start:pos + 1]

        # Reverse complement
        comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        seq_rc = ''.join(comp.get(b, 'N') for b in reversed(seq))

        if len(seq_rc) < 3:
            seq_rc = "N" * (3 - len(seq_rc)) + seq_rc
        trinuc = seq_rc[:3]

        # Determine context (same logic as forward, but on rev comp)
        if len(trinuc) >= 2 and trinuc[1] == 'G':
            context = "CG"
        elif len(trinuc) >= 3 and trinuc[2] == 'G':
            context = "CHG"
        else:
            context = "CHH"

    return (context, trinuc)


def infer_strand_from_reference(chrom_seq, pos):
    """Infer strand from the reference base at a 0-based position.

    Returns '+' if the base is C (cytosine on forward strand),
    '-' if the base is G (cytosine on reverse strand),
    or None if the base is neither (unexpected).
    """
    if pos < 0 or pos >= len(chrom_seq):
        return None
    base = chrom_seq[pos]
    if base == 'C':
        return '+'
    elif base == 'G':
        return '-'
    return None


def missing_contexts_error(sample, missing, counts):
    """Message for a dmC sample with no calls in a listed context."""
    have = ", ".join(f"{c} {n:,}" for c, n in counts.items() if n) or "none"
    return (
        f"Error: dmC sample {sample} has no covered "
        f"{'/'.join(missing)} position, but methylation_contexts lists "
        f"{'it' if len(missing) == 1 else 'them'} (calls per context: {have}).\n"
        "  Most ONT basecalling models call 5mC at CpG only (Dorado's "
        "5mCG_5hmCG models), so their bedMethyl has nothing in CHG or CHH, "
        "and every DMR, bigwig and PCA step in those contexts would have "
        "nothing to work on. Either:\n"
        "  - set methylation_contexts to the contexts the ONT/PacBio calls "
        "cover (e.g. [\"CG\"]); it applies to every mC sample in the run; or\n"
        "  - call the reads with an all-context model (Dorado's 5mC_5hmC "
        "models) and rerun modkit pileup.\n"
    )


def convert_bedmethyl_to_cx_report(input_file, reference_file, output_file,
                                   require_contexts=None, sample=None):
    """Convert bedMethyl to CX_report format with context from reference.

    Filters to 5mC calls only (mod_code 'm') or combined calls ('C').
    Infers strand from reference when bedMethyl has combined strands.

    ``require_contexts`` lists contexts (from ``methylation_contexts``) that
    must each have at least one covered call. If any has none, the error is
    printed and ``SystemExit(1)`` raised after the conversion. Returns the
    number of covered calls per context.
    """

    fasta = pysam.FastaFile(reference_file)

    # Cache for current chromosome sequence
    current_chrom = None
    chrom_seq = None
    chrom_len = 0

    n_kept = 0
    covered = {"CG": 0, "CHG": 0, "CHH": 0}
    n_skipped_mod = 0
    n_skipped_strand = 0
    n_skipped_no_ref = 0

    with ExitStack() as stack:
        # Handle gzipped or plain input
        if input_file.endswith('.gz'):
            infile = stack.enter_context(gzip.open(input_file, 'rt'))
        else:
            infile = stack.enter_context(open(input_file, 'r'))

        # Handle stdout or file output
        if output_file == '/dev/stdout':
            outfile = sys.stdout
        elif output_file.endswith('.gz'):
            outfile = stack.enter_context(gzip.open(output_file, 'wt'))
        else:
            outfile = stack.enter_context(open(output_file, 'w'))

        for line in infile:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 13:
                continue

            # Filter by modification type: keep only 5mC ('m') or
            # combined ('C', from --combine-mods)
            mod_code = fields[3]
            if mod_code not in ('m', 'C'):
                n_skipped_mod += 1
                continue

            chrom = fields[0]
            start = int(fields[1])
            strand = fields[5]
            n_mod = int(fields[11])
            n_canonical = int(fields[12])

            # Load chromosome sequence if we've moved to a new chromosome
            if chrom != current_chrom:
                try:
                    chrom_seq = fasta.fetch(chrom).upper()
                    chrom_len = len(chrom_seq)
                    current_chrom = chrom
                    print(f"Loaded chromosome {chrom} ({chrom_len:,} bp)", file=sys.stderr)
                except KeyError:
                    print(f"Warning: chromosome {chrom} not found in reference", file=sys.stderr)
                    chrom_seq = ""
                    chrom_len = 0
                    current_chrom = chrom

            # Skip positions on chromosomes absent from the reference: their
            # context cannot be determined, and emitting a placeholder "CNN"
            # leaks into the bigwig demux (mC.smk), where the catch-all branch
            # would miscount them as CG methylation.
            if chrom_len == 0:
                n_skipped_no_ref += 1
                continue

            # Infer strand from reference when bedMethyl has combined strands
            if strand == '.':
                strand = infer_strand_from_reference(chrom_seq, start)
                if strand is None:
                    n_skipped_strand += 1
                    continue

            # Get context and trinucleotide from cached sequence
            context, trinuc = get_context_and_trinuc_from_seq(chrom_seq, chrom_len, start, strand)

            # Convert to 1-based position for CX_report
            position = start + 1

            # Write CX_report format
            outfile.write(f"{chrom}\t{position}\t{strand}\t{n_mod}\t{n_canonical}\t{context}\t{trinuc}\n")
            n_kept += 1
            if n_mod + n_canonical > 0 and context in covered:
                covered[context] += 1

    fasta.close()

    print(f"Kept {n_kept:,} 5mC calls", file=sys.stderr)
    if n_skipped_mod:
        print(f"Skipped {n_skipped_mod:,} non-5mC rows (e.g. 5hmC)", file=sys.stderr)
    if n_skipped_strand:
        print(f"Skipped {n_skipped_strand:,} positions where reference base was neither C nor G", file=sys.stderr)
    if n_skipped_no_ref:
        print(f"Skipped {n_skipped_no_ref:,} positions on chromosomes absent from the reference", file=sys.stderr)
    print("Covered calls per context: "
          + ", ".join(f"{c} {n:,}" for c, n in covered.items()), file=sys.stderr)

    # Each context methylation_contexts lists must have data, or every
    # per-context step downstream would run on nothing (metilene segfaults on
    # an empty group). A CpG-only ONT model is the usual cause, so say so
    # here, at the first step that knows the contexts, instead of at DMRs.
    missing = [c for c in (require_contexts or []) if covered.get(c, 0) == 0]
    if missing:
        print(missing_contexts_error(sample or input_file, missing, covered),
              file=sys.stderr)
        raise SystemExit(1)
    return covered


def main():
    if len(sys.argv) not in (4, 6):
        print(f"Usage: {sys.argv[0]} input.bed.gz reference.fa output.CX_report.txt.gz "
              "[required_contexts sample_name]\n"
              "  required_contexts: '|'-separated, e.g. 'CG|CHG|CHH'")
        sys.exit(1)

    input_file = sys.argv[1]
    reference_file = sys.argv[2]
    output_file = sys.argv[3]
    require = sys.argv[4].split("|") if len(sys.argv) == 6 else None
    sample = sys.argv[5] if len(sys.argv) == 6 else None

    convert_bedmethyl_to_cx_report(input_file, reference_file, output_file,
                                   require_contexts=require, sample=sample)


if __name__ == "__main__":
    main()
