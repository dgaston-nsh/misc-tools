#!/usr/bin/env bash
# =============================================================================
# run_salmon_quant.sh
# Run salmon quant on paired-end RNA-seq samples in a directory.
#
# Sample ID format: YY-XXXXAB R_
#   YY   = two-digit year (e.g. 26, 25, 19)
#   XXXX = 2–4 digit sample number, optionally zero-padded (e.g. 0042, 123)
#   AB   = two letters for first/last initial
#   R    = RNA (literal)
#
# Usage:
#   bash run_salmon_quant.sh <fastq_dir> <output_dir>
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# !! CONFIGURE THESE BEFORE RUNNING !!
# ---------------------------------------------------------------------------
SALMON_INDEX="/mnt/shared-data/Resources/RNA/Salmon/salmon_index"   # Path to your pre-built salmon index
THREADS=8                               # Number of threads to use
LIBTYPE="A"                             # Library type; "A" = auto-detect
# Optional extra salmon quant arguments (leave empty string to skip)
# ---------------------------------------------------------------------------

# --- Colour helpers ---------------------------------------------------------
RED='\033[0;31m'; YELLOW='\033[1;33m'; GREEN='\033[0;32m'; NC='\033[0m'
info()  { echo -e "${GREEN}[INFO]${NC}  $*"; }
warn()  { echo -e "${YELLOW}[WARN]${NC}  $*"; }
error() { echo -e "${RED}[ERROR]${NC} $*" >&2; }

# --- Argument handling ------------------------------------------------------
if [[ $# -ne 2 ]]; then
    echo "Usage: $(basename "$0") <fastq_dir> <output_dir>"
    exit 1
fi

FASTQ_DIR="${1%/}"   # strip trailing slash if present
OUT_DIR="${2%/}"

# --- Pre-flight checks ------------------------------------------------------
if [[ ! -d "$FASTQ_DIR" ]]; then
    error "FastQ directory not found: $FASTQ_DIR"
    exit 1
fi

if [[ "$SALMON_INDEX" == "/path/to/salmon_index" ]]; then
    error "SALMON_INDEX is still set to the placeholder value. Edit the script first."
    exit 1
fi

if [[ ! -d "$SALMON_INDEX" ]]; then
    error "Salmon index directory not found: $SALMON_INDEX"
    exit 1
fi

if ! command -v salmon &>/dev/null; then
    error "salmon is not found in PATH. Please install or activate the relevant environment."
    exit 1
fi

mkdir -p "$OUT_DIR"

# ---------------------------------------------------------------------------
# Regex for the sample ID prefix: YY-XXXXABR_
#   ^([0-9]{2}-[0-9]{2,4}[A-Za-z]{2}R)_
# Positive control HD789R is detected separately as a fixed special case.
# ---------------------------------------------------------------------------
ID_PATTERN='^([0-9]{2}-[0-9]{2,4}[A-Za-z]{2}R)_'
POSITIVE_CONTROL_ID="HD789R"

# ---------------------------------------------------------------------------
# Discover unique sample IDs from R1 files
# We look for files whose basename matches the pattern and contain _R1_ or _1.
# Adjust the glob/grep below if your R1 naming convention differs.
# ---------------------------------------------------------------------------
declare -A seen_ids   # associative array used as a set

info "Scanning $FASTQ_DIR for samples..."

while IFS= read -r -d '' fq; do
    basename_fq="$(basename "$fq")"
    if [[ "$basename_fq" =~ $ID_PATTERN ]]; then
        sample_id="${BASH_REMATCH[1]}"
        seen_ids["$sample_id"]=1
    elif [[ "$basename_fq" == "${POSITIVE_CONTROL_ID}_"* ]]; then
        seen_ids["$POSITIVE_CONTROL_ID"]=1
    fi
done < <(find "$FASTQ_DIR" -maxdepth 1 -type f \
            \( -name "*.fastq.gz" -o -name "*.fq.gz" \
               -o -name "*.fastq"  -o -name "*.fq" \) \
            -print0 | sort -z)

if [[ ${#seen_ids[@]} -eq 0 ]]; then
    error "No samples matching the expected ID pattern (or positive control $POSITIVE_CONTROL_ID) were found in $FASTQ_DIR"
    exit 1
fi

info "Found ${#seen_ids[@]} sample(s): ${!seen_ids[*]}"
echo

# ---------------------------------------------------------------------------
# Process one sample at a time
# ---------------------------------------------------------------------------
PASS=0
FAIL=0
SKIPPED=0

for sample_id in $(printf '%s\n' "${!seen_ids[@]}" | sort); do

    info "Processing sample: $sample_id"

    # -----------------------------------------------------------------------
    # Find the two FastQ files for this sample.
    # Supports common naming conventions:
    #   *_R1_001.fastq.gz / *_R2_001.fastq.gz
    #   *_R1.fastq.gz     / *_R2.fastq.gz
    #   *_1.fastq.gz      / *_2.fastq.gz
    # -----------------------------------------------------------------------
    # Patterns cover both:
    #   Direct:  <ID>_R1_001.fastq.gz  (no intervening fields)
    #   Illumina: <ID>_S##_L###_R1_001.fastq.gz  (sample + lane fields before R1/R2)
    mapfile -t r1_files < <(find "$FASTQ_DIR" -maxdepth 1 -type f \
        \( -name "${sample_id}_*_R1_*.fastq.gz" -o -name "${sample_id}_*_R1_*.fq.gz" \
           -o -name "${sample_id}_*_R1_*.fastq"  -o -name "${sample_id}_*_R1_*.fq"   \
           -o -name "${sample_id}_R1_*.fastq.gz" -o -name "${sample_id}_R1_*.fq.gz"  \
           -o -name "${sample_id}_R1_*.fastq"    -o -name "${sample_id}_R1_*.fq"     \
           -o -name "${sample_id}_R1.fastq.gz"   -o -name "${sample_id}_R1.fq.gz"    \
           -o -name "${sample_id}_R1.fastq"      -o -name "${sample_id}_R1.fq"       \
           -o -name "${sample_id}_1.fastq.gz"    -o -name "${sample_id}_1.fq.gz"     \
           -o -name "${sample_id}_1.fastq"       -o -name "${sample_id}_1.fq"        \
        \) | sort)

    mapfile -t r2_files < <(find "$FASTQ_DIR" -maxdepth 1 -type f \
        \( -name "${sample_id}_*_R2_*.fastq.gz" -o -name "${sample_id}_*_R2_*.fq.gz" \
           -o -name "${sample_id}_*_R2_*.fastq"  -o -name "${sample_id}_*_R2_*.fq"   \
           -o -name "${sample_id}_R2_*.fastq.gz" -o -name "${sample_id}_R2_*.fq.gz"  \
           -o -name "${sample_id}_R2_*.fastq"    -o -name "${sample_id}_R2_*.fq"     \
           -o -name "${sample_id}_R2.fastq.gz"   -o -name "${sample_id}_R2.fq.gz"    \
           -o -name "${sample_id}_R2.fastq"      -o -name "${sample_id}_R2.fq"       \
           -o -name "${sample_id}_2.fastq.gz"    -o -name "${sample_id}_2.fq.gz"     \
           -o -name "${sample_id}_2.fastq"       -o -name "${sample_id}_2.fq"        \
        \) | sort)

    # Validate we found exactly one R1 and one R2
    if [[ ${#r1_files[@]} -ne 1 || ${#r2_files[@]} -ne 1 ]]; then
        warn "Expected exactly 1 R1 and 1 R2 for $sample_id; found ${#r1_files[@]} R1 and ${#r2_files[@]} R2. Skipping."
        (( SKIPPED++ )) || true
        continue
    fi

    R1="${r1_files[0]}"
    R2="${r2_files[0]}"
    SAMPLE_OUT="$OUT_DIR/$sample_id"

    # Skip if output already exists (remove this block to force re-runs)
    if [[ -f "$SAMPLE_OUT/quant.sf" ]]; then
        warn "Output already exists for $sample_id ($SAMPLE_OUT/quant.sf). Skipping."
        (( SKIPPED++ )) || true
        continue
    fi

    info "  R1 : $R1"
    info "  R2 : $R2"
    info "  Out: $SAMPLE_OUT"

    # -----------------------------------------------------------------------
    # Run salmon quant
    # -----------------------------------------------------------------------
    # shellcheck disable=SC2086
    if salmon quant \
        --index  "$SALMON_INDEX" \
        --libType "$LIBTYPE" \
        --mates1  "$R1" \
        --mates2  "$R2" \
        --output  "$SAMPLE_OUT" \
        --threads "$THREADS" \
        --gcBias \
        2>&1 | tee "$SAMPLE_OUT/../${sample_id}_salmon.log"; then
        info "  ✓ Completed: $sample_id"
        (( PASS++ )) || true
    else
        error "  ✗ salmon quant failed for $sample_id (exit code $?)"
        (( FAIL++ )) || true
    fi

    echo  # blank line between samples
done

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo "========================================="
info "Done. Results:"
echo "  Completed : $PASS"
echo "  Failed    : $FAIL"
echo "  Skipped   : $SKIPPED"
echo "========================================="

[[ $FAIL -eq 0 ]]  # exit 1 if any sample failed
