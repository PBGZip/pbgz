# shellcheck shell=bash
# Ensure TMPDIR exists and is writable.  autom4te and config.guess create
# temporary files there; a stale TMPDIR such as /tmp/$USER can otherwise make
# automake report "The installed version of autoconf does not work" even when
# autoconf itself is fine.
_tmp="${TMPDIR:-/tmp}"
if [ "$_tmp" != "/tmp" ]; then
    mkdir -p "$_tmp" 2>/dev/null || true
fi
if [ ! -d "$_tmp" ] || [ ! -w "$_tmp" ]; then
    _tmp=/tmp
fi
export TMPDIR="$_tmp"
unset _tmp

# ---------------------------------------------------------------------------
# ARM64 ISA level for the third-party libraries.
#
# Keep ISA-L (NEON inflate) and htslib/htscodecs (NEON rANS) at the same -march
# as the pbgz binary itself (see cmake/BuildOptions.cmake): the release tree is
# always built natively on the deployment host, so probing the running CPU is
# correct.  Cross builds and non-Linux hosts keep the portable baseline.
#   DotProd + FP16 -> armv8.2-a+dotprod+fp16+nolse  (LSE disabled on purpose)
#   DotProd only   -> armv8.2-a+dotprod+nolse
#   CRC32 only     -> armv8-a+crc
#   otherwise      -> compiler default (armv8-a)
#
# Set PBGZ_ARM_MARCH in the environment to override the probe (use the same
# value as the CMake cache variable so the binary and these libraries agree),
# e.g.  PBGZ_ARM_MARCH=armv8-a+crc ./build.all.gcc.sh
# ---------------------------------------------------------------------------
if [ -z "${PBGZ_ARM_MARCH+set}" ]; then
    PBGZ_ARM_MARCH=""
    if [ "$(uname -s)" = "Linux" ] &&
       { [ "$(uname -m)" = "aarch64" ] || [ "$(uname -m)" = "arm64" ]; }; then
        _pbgz_feat=$(sed -n 's/^Features[[:space:]]*:[[:space:]]*//p' /proc/cpuinfo 2>/dev/null | head -1)
        case " $_pbgz_feat " in
            *" asimddp "*)
                case " $_pbgz_feat " in
                    *" fphp "*|*" asimdhp "*) PBGZ_ARM_MARCH="-march=armv8.2-a+dotprod+fp16+nolse" ;;
                    *)                        PBGZ_ARM_MARCH="-march=armv8.2-a+dotprod+nolse" ;;
                esac
                ;;
            *" crc32 "*)
                PBGZ_ARM_MARCH="-march=armv8-a+crc"
                ;;
        esac
        unset _pbgz_feat
    fi
elif [ -n "$PBGZ_ARM_MARCH" ]; then
    # Accept both "armv8-a+crc" and "-march=armv8-a+crc".
    case "$PBGZ_ARM_MARCH" in
        -*) : ;;
        *)  PBGZ_ARM_MARCH="-march=$PBGZ_ARM_MARCH" ;;
    esac
fi
export PBGZ_ARM_MARCH
