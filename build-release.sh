#!/bin/bash
set -e

# Get script directory
ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# Default compiler and jobs
COMPILER="gcc"
JOBS=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 2)

# PGO is opt-in: --pgo (or PBGZ_PGO=1) plus --pgo-data (or PBGZ_PGO_DATA).
PGO_ENABLED=0
case "${PBGZ_PGO:-0}" in
  1|on|true|yes|ON|TRUE|YES) PGO_ENABLED=1 ;;
  0|off|false|no|"") PGO_ENABLED=0 ;;
  *)
    echo "Error: PBGZ_PGO must be 1/0 (or on/off); got '$PBGZ_PGO'"
    exit 1
    ;;
esac
PGO_DATA=${PBGZ_PGO_DATA:-}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --compiler=*)
      COMPILER="${1#*=}"
      shift
      ;;
    --jobs=*)
      JOBS="${1#*=}"
      shift
      ;;
    -j)
      if [[ $# -gt 1 ]]; then
        JOBS="$2"
        shift 2
      else
        echo "Error: -j requires a number"
        exit 1
      fi
      ;;
    --pgo)
      PGO_ENABLED=1
      shift
      ;;
    --pgo-data=*)
      PGO_DATA="${1#*=}"
      shift
      ;;
    --no-pgo)
      PGO_ENABLED=0
      shift
      ;;
    --help)
      echo "Usage: $0 [--compiler=gcc|clang|icc] [--jobs=N|-j N] [--pgo --pgo-data=DIR]"
      echo "Build Release version of PBGZ."
      echo ""
      echo "PGO is disabled by default. To build with PGO on ARM/GCC, provide a"
      echo "training data directory containing the reference and reads:"
      echo "  NA21087.chrom11_1079997.sam + hs37d5.fa.gz"
      echo "  and/or Arab_test.fq + GCA_028009825.1.fasta"
      echo ""
      echo "Options:"
      echo "  --compiler=NAME    Select compiler (gcc, clang, icc)"
      echo "  --jobs=N, -j N     Number of parallel build jobs (default: auto)"
      echo "  --pgo              Enable PGO (requires --pgo-data or PBGZ_PGO_DATA)"
      echo "  --pgo-data=DIR     Directory with the PGO training data"
      echo "  --no-pgo           Disable PGO (default)"
      echo ""
      echo "Environment:"
      echo "  PBGZ_PGO=1|0       Equivalent to --pgo / --no-pgo"
      echo "  PBGZ_PGO_DATA=DIR  Equivalent to --pgo-data=DIR"
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      echo "Use --help for usage information"
      exit 1
      ;;
  esac
done

# Determine architecture
ARCH=$(uname -m)
if [[ "$ARCH" == "x86_64" ]]; then
    ARCH="x86_64"
elif [[ "$ARCH" == "arm64" || "$ARCH" == "aarch64" ]]; then
    ARCH="arm64"
else
    echo "Unknown architecture: $ARCH"
    exit 1
fi

# Check dependencies
echo "Checking dependencies..."
source "$ROOT_DIR/scripts/check_dependencies.sh"

# Set compiler
case $COMPILER in
  gcc)
    export CC=gcc
    export CXX=g++
    ;;
  clang)
    export CC=clang
    export CXX=clang++
    ;;
  icc)
    export CC=icc
    export CXX=icpc
    ;;
  *)
    echo "Unsupported compiler: $COMPILER"
    echo "Supported compilers: gcc, clang, icc"
    exit 1
    ;;
esac

# Check if compiler exists
if ! command -v $CC &> /dev/null; then
    echo "Error: Compiler $CC not found"
    exit 1
fi

BUILD_DIR="$ROOT_DIR/build-$ARCH-$COMPILER-release"
PGO_BUILD_DIR="$ROOT_DIR/build-$ARCH-$COMPILER-pgo"
INSTALL_DIR="$ROOT_DIR/release-release"

PGO_SAM="$PGO_DATA/NA21087.chrom11_1079997.sam"
PGO_REF_GZ="$PGO_DATA/hs37d5.fa.gz"
PGO_FQ="$PGO_DATA/Arab_test.fq"
PGO_FA_REF="$PGO_DATA/GCA_028009825.1.fasta"

PGO_TMP=""
cleanup_pgo_tmp() {
    if [[ -n "$PGO_TMP" && -d "$PGO_TMP" ]]; then
        rm -rf "$PGO_TMP"
    fi
}
trap cleanup_pgo_tmp EXIT

build_lto() {
    echo "Creating LTO/normal build directory: $BUILD_DIR"
    mkdir -p "$BUILD_DIR"
    echo "Configuring CMake project (PBGZ_PGO=OFF)..."
    cmake -S "$ROOT_DIR" -B "$BUILD_DIR" \
          -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" \
          -DPBGZ_PGO=OFF
    echo "Building project with $JOBS parallel jobs..."
    cmake --build "$BUILD_DIR" -- VERBOSE=1 -j"$JOBS"
    echo "Installing to $INSTALL_DIR..."
    mkdir -p "$INSTALL_DIR"
    cmake --install "$BUILD_DIR"
    echo "Build complete. Executable located at: $INSTALL_DIR/bin/pbgz"
}

pgo_supported() {
    [[ "$ARCH" == "arm64" ]] || return 1
    [[ "$COMPILER" == "gcc" ]] || return 1
    local probe
    probe=$(mktemp /tmp/pbgz_pgo_probe.XXXXXX) || return 1
    if echo 'int main(){return 0;}' | "$CXX" -x c++ -O2 -fprofile-generate - -o "$probe" >/dev/null 2>&1; then
        rm -f "$probe" "$probe".gcda
        return 0
    fi
    rm -f "$probe" "$probe".gcda
    return 1
}

validate_pgo_data() {
    if [[ -z "$PGO_DATA" ]]; then
        echo "Error: --pgo requires --pgo-data=DIR (or PBGZ_PGO_DATA=DIR)." >&2
        return 1
    fi
    if [[ ! -d "$PGO_DATA" ]]; then
        echo "Error: PGO data directory does not exist: $PGO_DATA" >&2
        return 1
    fi
    PGO_DATA=$(cd "$PGO_DATA" && pwd)
    PGO_SAM="$PGO_DATA/NA21087.chrom11_1079997.sam"
    PGO_REF_GZ="$PGO_DATA/hs37d5.fa.gz"
    PGO_FQ="$PGO_DATA/Arab_test.fq"
    PGO_FA_REF="$PGO_DATA/GCA_028009825.1.fasta"

    local have_sam=0 have_fq=0
    [[ -f "$PGO_SAM" && -f "$PGO_REF_GZ" ]] && have_sam=1
    [[ -f "$PGO_FQ" && -f "$PGO_FA_REF" ]] && have_fq=1
    if [[ $have_sam -eq 0 && $have_fq -eq 0 ]]; then
        echo "Error: no usable training set in $PGO_DATA" >&2
        echo "Expected at least one pair:" >&2
        echo "  NA21087.chrom11_1079997.sam + hs37d5.fa.gz" >&2
        echo "  Arab_test.fq + GCA_028009825.1.fasta" >&2
        return 1
    fi
    return 0
}

run_pgo_training() {
    local bin="$INSTALL_DIR/bin/pbgz"
    local trained=0

    if [[ -f "$PGO_SAM" && -f "$PGO_REF_GZ" ]]; then
        echo "[pgo] makeref $PGO_REF_GZ"
        "$bin" makeref "$PGO_REF_GZ" -o "$PGO_TMP/train.ni" -t "$JOBS" -f || return 1
        echo "[pgo] SAM compress/decompress"
        "$bin" compress "$PGO_SAM" -o "$PGO_TMP/train.sam.pbgz" \
               -r "$PGO_REF_GZ" --ni "$PGO_TMP/train.ni" -t "$JOBS" -f || return 1
        "$bin" decompress "$PGO_TMP/train.sam.pbgz" -o "$PGO_TMP/train.sam" \
               -r "$PGO_REF_GZ" --ni "$PGO_TMP/train.ni" -t "$JOBS" -f || return 1
        rm -f "$PGO_TMP/train.sam.pbgz" "$PGO_TMP/train.sam"
        trained=1
    fi

    if [[ -f "$PGO_FQ" && -f "$PGO_FA_REF" ]]; then
        echo "[pgo] FASTQ compress/decompress"
        "$bin" compress "$PGO_FQ" -o "$PGO_TMP/train.fq.pbgz" \
               -r "$PGO_FA_REF" -t "$JOBS" -f || return 1
        "$bin" decompress "$PGO_TMP/train.fq.pbgz" -o "$PGO_TMP/train.fq" \
               -r "$PGO_FA_REF" -t "$JOBS" -f || return 1
        rm -f "$PGO_TMP/train.fq.pbgz" "$PGO_TMP/train.fq"
        trained=1
    fi

    [[ $trained -eq 1 ]]
}

build_pgo() {
    PGO_TMP=$(mktemp -d /tmp/pbgz_pgo.XXXXXX) || return 1
    rm -rf "$PGO_BUILD_DIR"
    echo "Configuring PGO instrumentation build in $PGO_BUILD_DIR..."
    cmake -S "$ROOT_DIR" -B "$PGO_BUILD_DIR" \
          -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" \
          -DPBGZ_PGO=GENERATE || return 1

    echo "Building instrumented pbgz..."
    cmake --build "$PGO_BUILD_DIR" --target pbgz -- -j"$JOBS" || return 1

    # The build tree is not executable on some hardened hosts; run the
    # instrumented binary from the install directory instead.
    mkdir -p "$INSTALL_DIR/bin"
    cp -f "$PGO_BUILD_DIR/src/pbgz" "$INSTALL_DIR/bin/pbgz" || return 1

    echo "Running PGO training..."
    run_pgo_training || return 1

    local profile_count
    profile_count=$(find "$PGO_BUILD_DIR" -name '*.gcda' | wc -l)
    echo "PGO profile files: $profile_count"
    [[ "$profile_count" -gt 0 ]] || return 1

    # Keep the profiles while forcing every object through -fprofile-use.
    find "$PGO_BUILD_DIR" -name '*.o' -delete
    echo "Configuring profile-use build..."
    cmake -S "$ROOT_DIR" -B "$PGO_BUILD_DIR" \
          -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" \
          -DPBGZ_PGO=USE || return 1
    echo "Building optimized PGO binary..."
    cmake --build "$PGO_BUILD_DIR" --target pbgz -- -j"$JOBS" || return 1
    cp -f "$PGO_BUILD_DIR/src/pbgz" "$INSTALL_DIR/bin/pbgz" || return 1
    return 0
}

if [[ "$PGO_ENABLED" == "1" ]]; then
    if ! validate_pgo_data; then
        exit 1
    fi
    if ! pgo_supported; then
        echo "Error: PGO requires ARM + GCC with -fprofile-generate support." >&2
        exit 1
    fi
    echo "PGO enabled (training data: $PGO_DATA)."
    if ! build_pgo; then
        echo "Error: PGO build failed." >&2
        exit 1
    fi
    echo "PGO build complete. Executable located at: $INSTALL_DIR/bin/pbgz"
else
    if [[ -n "$PGO_DATA" ]]; then
        echo "Note: PBGZ_PGO_DATA / --pgo-data is ignored without --pgo (or PBGZ_PGO=1)."
    fi
    echo "PGO disabled; using the normal build."
    build_lto
fi
