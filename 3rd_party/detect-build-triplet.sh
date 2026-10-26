#!/bin/bash
# Echo an autoconf/config.sub-compatible build triplet.
#
# Some hosts (notably newer AMD servers) make autoconf's config.guess return an
# empty value, so every bundled ./configure aborts with:
#   "cannot guess build type; you must specify one"
# The compiler already knows the triplet, so query it first and only fall back
# to the common uname mappings when the compiler cannot report it.
triplet=$(${CC:-cc} -dumpmachine 2>/dev/null || true)
if [ -z "$triplet" ]; then
    case "$(uname -s)-$(uname -m)" in
        Linux-x86_64|Linux-amd64)  triplet=x86_64-pc-linux-gnu ;;
        Linux-aarch64|Linux-arm64) triplet=aarch64-unknown-linux-gnu ;;
        Darwin-x86_64)             triplet=x86_64-apple-darwin ;;
        Darwin-arm64)              triplet=aarch64-apple-darwin ;;
    esac
fi
printf '%s\n' "$triplet"
