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
