#!/bin/bash
# Echo an autoconf/config.sub-compatible build triplet.
#
# Some hosts (notably newer AMD / Hygon servers) make autoconf's config.guess
# return an empty value, so every bundled ./configure aborts with:
#   "cannot guess build type; you must specify one"
# Passing --build=<triplet> makes configure use that value directly and avoids
# config.guess entirely.
#
# The triplet must also be accepted by the *oldest* config.sub in this tree
# (m4-1.4.9), which predates aarch64 and rejects several modern spellings.
# Normalize the known aarch64 forms before returning them:
#   aarch64-linux-gnu        -> aarch64-unknown-linux-gnu
#   arm64-apple-darwin...    -> aarch64-unknown-darwin
# These canonical forms are accepted by m4 1.4.9, autoconf 2.69, automake 1.17,
# libtool 2.4.6 and yasm 1.3.0 on both Linux and Apple hosts.
#
# The script always exits 0; an empty line means the caller should omit --build.

triplet=$(${CC:-cc} -dumpmachine 2>/dev/null | tr -d '[:space:]')

# Apple clang reports arm64-...; GNU config.sub canonicalizes it to aarch64.
case "$triplet" in
    arm64-*) triplet="aarch64-${triplet#arm64-}" ;;
esac

case "$triplet" in
    aarch64-*linux*)  triplet="aarch64-unknown-linux-gnu" ;;
    aarch64-*darwin*) triplet="aarch64-unknown-darwin" ;;
esac

# Compiler-less fallback.  Keep the same normalized spellings so old config.sub
# versions accept them as well.
if [ -z "$triplet" ]; then
    case "$(uname -s)-$(uname -m)" in
        Linux-x86_64|Linux-amd64)   triplet=x86_64-unknown-linux-gnu ;;
        Linux-aarch64|Linux-arm64)  triplet=aarch64-unknown-linux-gnu ;;
        Darwin-x86_64)              triplet=x86_64-apple-darwin ;;
        Darwin-arm64)               triplet=aarch64-unknown-darwin ;;
    esac
fi

printf '%s\n' "$triplet"
