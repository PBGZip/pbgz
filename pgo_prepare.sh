#!/bin/bash
#
# pgo_prepare.sh - build the PGO training set used by build-release.sh
#
# build-release.sh --pgo only *consumes* training data: it reads train.tsv and
# profiles the workloads listed there. This script produces that data:
#   1. the reference indexes (*.ni) via `pbgz makeref`, unless one already exists
#   2. the train.tsv manifest (one workload per line)
#
# Usage:
#   ./pgo_prepare.sh --out DIR --pair INPUT REF [--pair INPUT REF]...
#
#   INPUT  a SAM/BAM/FASTQ file (any input pbgz accepts)
#   REF    the reference FASTA/FASTA.gz for that workload, or '-' for none
#
# Options:
#   --out DIR         output directory for train.tsv (required)
#   --pbgz PATH       pbgz binary used for makeref
#                     (default: ./release-release/bin/pbgz, else pbgz from PATH)
#   --jobs N          parallel jobs for makeref (default: nproc)
#   --index-dir DIR   where to put/find the *.ni indexes
#                     (default: next to the reference, i.e. REF.ni)
#   --force           rebuild the index even if it already exists
#   --dry-run         print what would be done, write nothing
#   -h, --help        this help
#
# Example:
#   ./pgo_prepare.sh --out /data/pgo --index-dir /data/refs \
#       --pair /data/reads.sam /data/refs/hs37d5.fa.gz \
#       --pair /data/reads.fq  /data/refs/genome.fna.gz
#   bash build-release.sh --pgo --pgo-data=/data/pgo
#
set -e

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

OUT=""
PBGZ_BIN=""
JOBS=$(nproc 2>/dev/null || echo 2)
INDEX_DIR=""
FORCE=0
DRY_RUN=0
PAIR_INPUTS=()
PAIR_REFS=()

usage() {
    sed -n '3,32p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --out)       OUT="${2:-}"; shift 2 ;;
        --out=*)     OUT="${1#*=}"; shift ;;
        --pbgz)      PBGZ_BIN="${2:-}"; shift 2 ;;
        --pbgz=*)    PBGZ_BIN="${1#*=}"; shift ;;
        --jobs|-j)   JOBS="${2:-}"; shift 2 ;;
        --jobs=*)    JOBS="${1#*=}"; shift ;;
        --index-dir) INDEX_DIR="${2:-}"; shift 2 ;;
        --index-dir=*) INDEX_DIR="${1#*=}"; shift ;;
        --force)     FORCE=1; shift ;;
        --dry-run)   DRY_RUN=1; shift ;;
        -h|--help)   usage; exit 0 ;;
        --pair)
            if [[ $# -lt 3 ]]; then
                echo "Error: --pair needs two arguments: INPUT REF" >&2
                exit 1
            fi
            PAIR_INPUTS+=("$2")
            PAIR_REFS+=("$3")
            shift 3
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ -z "$OUT" ]]; then
    echo "Error: --out DIR is required." >&2
    exit 1
fi
if [[ ${#PAIR_INPUTS[@]} -eq 0 ]]; then
    echo "Error: at least one --pair INPUT REF is required." >&2
    exit 1
fi

# Resolve the pbgz binary.
if [[ -z "$PBGZ_BIN" ]]; then
    if [[ -x "$ROOT_DIR/release-release/bin/pbgz" ]]; then
        PBGZ_BIN="$ROOT_DIR/release-release/bin/pbgz"
    elif command -v pbgz >/dev/null 2>&1; then
        PBGZ_BIN=$(command -v pbgz)
    else
        echo "Error: no pbgz binary found; pass --pbgz PATH." >&2
        exit 1
    fi
fi
if [[ ! -x "$PBGZ_BIN" ]]; then
    echo "Error: pbgz binary is not executable: $PBGZ_BIN" >&2
    exit 1
fi
PBGZ_BIN=$(readlink -f "$PBGZ_BIN")

[[ -n "$INDEX_DIR" ]] && INDEX_DIR=$(readlink -f "$INDEX_DIR")

if [[ $DRY_RUN -eq 0 ]]; then
    mkdir -p "$OUT"
fi
OUT=$(readlink -f "$OUT" 2>/dev/null || echo "$OUT")

MANIFEST="$OUT/train.tsv"
LINES=()
declare -A INDEX_DONE=()
built=0
reused=0

for i in "${!PAIR_INPUTS[@]}"; do
    input="${PAIR_INPUTS[$i]}"
    ref="${PAIR_REFS[$i]}"

    if [[ ! -f "$input" ]]; then
        echo "Error: input not found: $input" >&2
        exit 1
    fi
    input_abs=$(readlink -f "$input")

    if [[ "$ref" == "-" || -z "$ref" ]]; then
        LINES+=("$input_abs"$'\t'"-"$'\t'"-")
        echo "[pair] $input_abs  (no reference)"
        continue
    fi

    if [[ ! -f "$ref" ]]; then
        echo "Error: reference not found: $ref" >&2
        exit 1
    fi
    ref_abs=$(readlink -f "$ref")
    if [[ -n "$INDEX_DIR" ]]; then
        ni="$INDEX_DIR/$(basename "$ref_abs").ni"
    else
        ni="$ref_abs.ni"
    fi

    if [[ -n "${INDEX_DONE[$ni]:-}" ]]; then
        : # already handled for an earlier pair
    elif [[ -f "$ni" && $FORCE -eq 0 ]]; then
        echo "[pair] $input_abs  (reuse index $ni)"
        INDEX_DONE[$ni]=1
        reused=$((reused + 1))
    else
        echo "[pair] $input_abs  -> makeref $ref_abs -o $ni"
        if [[ $DRY_RUN -eq 0 ]]; then
            mkdir -p "$(dirname "$ni")"
            "$PBGZ_BIN" makeref "$ref_abs" -o "$ni" -t "$JOBS" -f >/dev/null
        fi
        INDEX_DONE[$ni]=1
        built=$((built + 1))
    fi
    LINES+=("$input_abs"$'\t'"$ref_abs"$'\t'"$ni")
done

if [[ $DRY_RUN -eq 1 ]]; then
    echo
    echo "[dry-run] would write $MANIFEST with ${#LINES[@]} workload(s):"
    for line in "${LINES[@]}"; do
        printf '  %s\n' "${line//$'\t'/  }"
    done
    exit 0
fi

{
    echo "# PGO training set generated by pgo_prepare.sh"
    echo "# <input>	<reference>	<index>   ('-' means none; relative paths resolve against this file)"
    for line in "${LINES[@]}"; do
        printf '%s\n' "$line"
    done
} > "$MANIFEST"

echo
echo "Wrote $MANIFEST (${#LINES[@]} workload(s); $built index(es) built, $reused reused)."
echo "Now build with:"
echo "  bash build-release.sh --pgo --pgo-data=$OUT"
