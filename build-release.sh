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
      echo "PGO is disabled by default and only supported on ARM64 Linux + GCC."
      echo "To build with PGO, point --pgo-data"
      echo "at a directory containing train.tsv, which you prepare beforehand (the script"
      echo "does not build reference indexes or any training data). One workload per line:"
      echo "  <input.sam|input.bam|input.fq>  <reference|->  <index|->"
      echo "Example:"
      echo "  reads.sam   ref.fa.gz   ref.fa.gz.ni"
      echo "  reads.fq    ref.fa      -"
      echo ""
      echo "Options:"
      echo "  --compiler=NAME    Select compiler (gcc, clang, icc)"
      echo "  --jobs=N, -j N     Number of parallel build jobs (default: auto)"
      echo "  --pgo              Enable PGO (requires --pgo-data or PBGZ_PGO_DATA)"
      echo "  --pgo-data=DIR     Directory containing train.tsv (prepared by the caller)"
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
CC=$(command -v $CC)
CXX=$(command -v $CXX)

# CMake comes from PATH: the project only needs cmake >= 3.10 (see CMakeLists.txt),
# so a system package (e.g. /usr/local/bin/cmake) is enough - no bundled toolchain.
CMAKE_BIN=${CMAKE_BIN:-cmake}
if ! command -v "$CMAKE_BIN" &> /dev/null; then
    echo "Error: cmake not found in PATH; cmake >= 3.10 is required."
    exit 1
fi
CMAKE_VERSION=$("$CMAKE_BIN" --version | head -1 | awk '{print $3}')
if [ "$(printf '%s\n' 3.10 "$CMAKE_VERSION" | sort -V | head -1)" != "3.10" ]; then
    echo "Error: cmake >= 3.10 required, found $CMAKE_VERSION"
    exit 1
fi
echo "Using cmake $CMAKE_VERSION ($(command -v "$CMAKE_BIN"))"
echo "Using compiler $($CXX --version | head -1) ($CXX)"

BUILD_DIR="$ROOT_DIR/build-$ARCH-$COMPILER-release"
PGO_BUILD_DIR="$ROOT_DIR/build-$ARCH-$COMPILER-pgo"
INSTALL_DIR="$ROOT_DIR/release-release"

# A build tree configured with a different compiler cannot be reused; recreate it.
prepare_build_dir() {
    local dir="$1"
    if [ -f "$dir/CMakeCache.txt" ]; then
        local cached
        cached=$(sed -n 's/^CMAKE_CXX_COMPILER:[A-Z]*=//p' "$dir/CMakeCache.txt" | head -1)
        if [ -n "$cached" ] && [ "$cached" != "$CXX" ]; then
            echo "Build tree was configured with $cached; recreating it for $CXX"
            rm -rf "$dir"
        fi
    fi
    mkdir -p "$dir"
}

PGO_TMP=""
PGO_MANIFEST_RESOLVED=""
cleanup_pgo_tmp() {
    if [[ -n "$PGO_TMP" && -d "$PGO_TMP" ]]; then
        rm -rf "$PGO_TMP"
    fi
    if [[ -n "$PGO_MANIFEST_RESOLVED" && -f "$PGO_MANIFEST_RESOLVED" ]]; then
        rm -f "$PGO_MANIFEST_RESOLVED"
    fi
}
trap cleanup_pgo_tmp EXIT

build_lto() {
    echo "Creating LTO/normal build directory: $BUILD_DIR"
    prepare_build_dir "$BUILD_DIR"
    echo "Configuring CMake project (PBGZ_PGO=OFF)..."
    "$CMAKE_BIN" -S "$ROOT_DIR" -B "$BUILD_DIR" \
          -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" \
          -DCMAKE_C_COMPILER="$CC" \
          -DCMAKE_CXX_COMPILER="$CXX" \
          -DPBGZ_PGO=OFF
    echo "Building project with $JOBS parallel jobs..."
    "$CMAKE_BIN" --build "$BUILD_DIR" -- VERBOSE=1 -j"$JOBS"
    echo "Installing to $INSTALL_DIR..."
    mkdir -p "$INSTALL_DIR"
    "$CMAKE_BIN" --install "$BUILD_DIR"
    echo "Build complete. Executable located at: $INSTALL_DIR/bin/pbgz"
}

pgo_compiler_supports_generate() {
    local probe
    probe=$(mktemp /tmp/pbgz_pgo_probe.XXXXXX) || return 1
    if echo 'int main(){return 0;}' | "$CXX" -x c++ -O2 -fprofile-generate - -o "$probe" >/dev/null 2>&1; then
        rm -f "$probe" "$probe".gcda
        return 0
    fi
    rm -f "$probe" "$probe".gcda
    return 1
}

# The training set is a manifest that the caller prepares; this script never creates
# reference indexes or training data. One workload per line:
#   <input.sam|input.bam|input.fq>  <reference|->  <index|->
# Relative paths are resolved against the manifest's directory. Both a compress and a
# decompress run are profiled for every line (decompress benefits most from PGO).
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
    local manifest="$PGO_DATA/train.tsv"
    if [[ ! -f "$manifest" ]]; then
        echo "Error: $manifest not found." >&2
        echo "PGO training is driven by a manifest you prepare beforehand; this script does" >&2
        echo "not build reference indexes or any other training data. One workload per line:" >&2
        echo "  <input.sam|input.bam|input.fq>  <reference|->  <index|->" >&2
        echo "Example:" >&2
        echo "  reads.sam   ref.fa.gz   ref.fa.gz.ni" >&2
        echo "  reads.fq    ref.fa      -" >&2
        echo "Generate it with: ./pgo_prepare.sh --out $PGO_DATA \\" >&2
        echo "    --pair reads.sam ref.fa.gz [--pair reads.fq ref.fa]" >&2
        return 1
    fi

    local resolved
    resolved=$(mktemp "${TMPDIR:-/tmp}/pbgz_pgo_set.XXXXXX") || return 1
    : > "$resolved"
    local input ref index n=0 ip rp xp
    while read -r input ref index _; do
        [[ -z "${input:-}" || "$input" == \#* ]] && continue
        ip="$input"
        rp="${ref:-}"
        xp="${index:-}"
        [[ "$ip" != /* ]] && ip="$PGO_DATA/$ip"
        [[ -n "$rp" && "$rp" != "-" && "$rp" != /* ]] && rp="$PGO_DATA/$rp"
        [[ -n "$xp" && "$xp" != "-" && "$xp" != /* ]] && xp="$PGO_DATA/$xp"
        if [[ ! -f "$ip" ]]; then
            echo "Error: training input not found: $ip" >&2
            rm -f "$resolved"
            return 1
        fi
        if [[ -n "$rp" && "$rp" != "-" && ! -f "$rp" ]]; then
            echo "Error: training reference not found: $rp" >&2
            rm -f "$resolved"
            return 1
        fi
        if [[ -n "$xp" && "$xp" != "-" && ! -f "$xp" ]]; then
            echo "Error: training reference index not found: $xp" >&2
            rm -f "$resolved"
            return 1
        fi
        printf '%s\t%s\t%s\n' "$ip" "${rp:--}" "${xp:--}" >> "$resolved"
        n=$((n + 1))
    done < "$manifest"

    if [[ $n -eq 0 ]]; then
        echo "Error: $manifest contains no workloads." >&2
        rm -f "$resolved"
        return 1
    fi
    PGO_MANIFEST_RESOLVED="$resolved"
    echo "PGO training set: $n workload(s) from $manifest"
    return 0
}

run_pgo_training() {
    local bin="$INSTALL_DIR/bin/pbgz"
    local input ref index n=0
    local ref_args ni_args
    while IFS=$'\t' read -r input ref index; do
        [[ -z "${input:-}" ]] && continue
        ref_args=()
        ni_args=()
        [[ -n "${ref:-}" && "$ref" != "-" ]] && ref_args=(-r "$ref")
        [[ -n "${index:-}" && "$index" != "-" ]] && ni_args=(--ni "$index")
        n=$((n + 1))
        echo "[pgo] ($n) compress/decompress $input"
        "$bin" compress "$input" -o "$PGO_TMP/train.$n.pbgz" \
               "${ref_args[@]}" "${ni_args[@]}" -t "$JOBS" -f >/dev/null || return 1
        "$bin" decompress "$PGO_TMP/train.$n.pbgz" -o "$PGO_TMP/train.$n.out" \
               "${ref_args[@]}" "${ni_args[@]}" -t "$JOBS" -f >/dev/null || return 1
        rm -f "$PGO_TMP/train.$n.pbgz" "$PGO_TMP/train.$n.out"
    done < "$PGO_MANIFEST_RESOLVED"
    [[ $n -gt 0 ]]
}

build_pgo() {
    PGO_TMP=$(mktemp -d /tmp/pbgz_pgo.XXXXXX) || return 1
    rm -rf "$PGO_BUILD_DIR"
    echo "Configuring PGO instrumentation build in $PGO_BUILD_DIR..."
    "$CMAKE_BIN" -S "$ROOT_DIR" -B "$PGO_BUILD_DIR" \
          -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" \
          -DCMAKE_C_COMPILER="$CC" \
          -DCMAKE_CXX_COMPILER="$CXX" \
          -DPBGZ_PGO=GENERATE || return 1

    echo "Building instrumented pbgz..."
    "$CMAKE_BIN" --build "$PGO_BUILD_DIR" --target pbgz -- -j"$JOBS" || return 1

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
    "$CMAKE_BIN" -S "$ROOT_DIR" -B "$PGO_BUILD_DIR" \
          -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" \
          -DCMAKE_C_COMPILER="$CC" \
          -DCMAKE_CXX_COMPILER="$CXX" \
          -DPBGZ_PGO=USE || return 1
    echo "Building optimized PGO binary..."
    "$CMAKE_BIN" --build "$PGO_BUILD_DIR" --target pbgz -- -j"$JOBS" || return 1
    cp -f "$PGO_BUILD_DIR/src/pbgz" "$INSTALL_DIR/bin/pbgz" || return 1
    return 0
}

if [[ "$PGO_ENABLED" == "1" ]]; then
    # PGO is implemented and verified for ARM64 Linux + GCC only; refuse to build a
    # silently unprofiled binary on anything else.
    if [[ "$ARCH" != "arm64" ]]; then
        echo "Error: PGO is only supported on ARM64 (current: $ARCH)." >&2
        exit 1
    fi
    if [[ "$COMPILER" != "gcc" ]]; then
        echo "Error: PGO requires GCC (selected: $COMPILER)." >&2
        exit 1
    fi
    if ! pgo_compiler_supports_generate; then
        echo "Error: $CXX does not support -fprofile-generate; PGO unavailable." >&2
        exit 1
    fi
    if ! validate_pgo_data; then
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
