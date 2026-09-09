#!/usr/bin/env bash
#
# validate_pombe.sh — Orchestrate S. pombe pipeline validation stages.
#
# Usage:
#   bash scripts/validate_pombe.sh --dry     # Run dry-run integration tests
#   bash scripts/validate_pombe.sh --full    # Execute the full Snakemake pipeline
#   bash scripts/validate_pombe.sh --check   # Validate completed pipeline outputs
#   bash scripts/validate_pombe.sh --all     # All three in sequence
#
# On a cluster, --full needs a profile carrying your site's account/partition.
# Pass --profile NAME (or set EPICC_PROFILE); profiles/slurm is only a template.
#
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_DIR"

CONFIG_FILE="tests/integration/data/test_options_pombe.yaml"
# Execution profile for --full. Empty means "decide at run time"; see pick_profile.
PROFILE="${EPICC_PROFILE:-}"
SAMPLE_FILE="tests/integration/data/test_samples_pombe.tsv"

# Colors for output (if terminal supports it)
if [ -t 1 ]; then
    GREEN='\033[0;32m'
    RED='\033[0;31m'
    YELLOW='\033[0;33m'
    NC='\033[0m'
else
    GREEN='' RED='' YELLOW='' NC=''
fi

usage() {
    echo "Usage: $0 [--dry] [--full] [--check] [--all]"
    echo ""
    echo "  --dry    Run dry-run integration tests (pytest)"
    echo "  --full   Execute the full Snakemake pipeline"
    echo "  --check  Validate completed pipeline outputs (pytest)"
    echo "  --all    Run all three stages in sequence"
    echo "  --profile NAME  Snakemake profile for --full (default: \$EPICC_PROFILE,"
    echo "                  else \$SNAKEMAKE_PROFILE, else profiles/slurm)"
    echo ""
    echo "Stages abort on failure. --all stops before --full if --dry fails."
    exit 1
}

run_dry() {
    echo -e "${YELLOW}=== Stage: Dry-Run Tests ===${NC}"
    conda run -n epicc pytest tests/integration/test_pombe_dryrun.py -v
    echo -e "${GREEN}=== Dry-run tests PASSED ===${NC}"
}

run_full() {
    echo -e "${YELLOW}=== Stage: Full Pipeline Run ===${NC}"

    if ! [ -f "$CONFIG_FILE" ]; then
        echo -e "${RED}ERROR: Options file not found: $CONFIG_FILE${NC}" >&2
        exit 1
    fi
    if ! [ -f "$SAMPLE_FILE" ]; then
        echo -e "${RED}ERROR: Sample file not found: $SAMPLE_FILE${NC}" >&2
        exit 1
    fi

    if command -v sbatch &>/dev/null; then
        # profiles/slurm is a template with no account or partition, so picking
        # it just because sbatch exists sends jobs in without site settings.
        # Prefer an explicit choice; fall back with a warning that names the
        # site profiles actually present.
        local profile_args=()
        if [ -n "$PROFILE" ]; then
            echo -e "  SLURM detected — using profile: $PROFILE"
            profile_args=(--profile "$PROFILE")
        elif [ -n "${SNAKEMAKE_PROFILE:-}" ]; then
            echo -e "  SLURM detected — using \$SNAKEMAKE_PROFILE ($SNAKEMAKE_PROFILE)"
        else
            local others
            others="$(find profiles -mindepth 1 -maxdepth 1 -type d \
                        ! -name default ! -name slurm ! -name uge \
                        -printf '%f ' 2>/dev/null || true)"
            echo -e "${YELLOW}  SLURM detected — falling back to profiles/slurm," \
                    "a template with no account/partition.${NC}" >&2
            [ -n "$others" ] && echo -e "${YELLOW}  Site profiles available:" \
                    "${others}— pass --profile <name> to use one.${NC}" >&2
            profile_args=(--profile profiles/slurm)
        fi
        conda run -n epicc snakemake \
            "${profile_args[@]}" \
            --configfile "$CONFIG_FILE"
    else
        CORES=$(( $(nproc) / 2 ))
        [ "$CORES" -lt 1 ] && CORES=1
        echo -e "  No SLURM — running locally with $CORES cores"
        conda run -n epicc snakemake \
            --use-conda --conda-frontend conda \
            --cores "$CORES" \
            --configfile "$CONFIG_FILE"
    fi

    echo -e "${GREEN}=== Full pipeline run COMPLETED ===${NC}"
}

run_check() {
    echo -e "${YELLOW}=== Stage: Post-Run Validation ===${NC}"
    conda run -n epicc pytest tests/integration/test_pombe_postrun.py -v -m slow
    echo -e "${GREEN}=== Post-run validation PASSED ===${NC}"
}

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
if [ $# -eq 0 ]; then
    usage
fi

DO_DRY=false
DO_FULL=false
DO_CHECK=false

while [ $# -gt 0 ]; do
    case "$1" in
        --dry)   DO_DRY=true ;;
        --full)  DO_FULL=true ;;
        --check) DO_CHECK=true ;;
        --all)   DO_DRY=true; DO_FULL=true; DO_CHECK=true ;;
        --profile)
            [ $# -ge 2 ] || { echo -e "${RED}--profile needs a value${NC}" >&2; usage; }
            PROFILE="$2"; shift ;;
        --profile=*) PROFILE="${1#*=}" ;;
        -h|--help) usage ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}" >&2
            usage
            ;;
    esac
    shift
done

# ---------------------------------------------------------------------------
# Execute stages in order
# ---------------------------------------------------------------------------
if $DO_DRY; then
    run_dry
fi

if $DO_FULL; then
    run_full
fi

if $DO_CHECK; then
    run_check
fi

echo -e "${GREEN}=== All requested stages completed ===${NC}"
