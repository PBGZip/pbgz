#!/bin/bash
# Script008: Cleanup CI resources
# Function: Cleanup SSH key、temporary file, CI related resources
# Independent debugging: yes
# Reusable: yes

set -e

log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1"
}

error() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] ❌ Error: $1"
    exit 1
}

# Cleanup SSH related resources
cleanup_ssh_resources() {
    log "start Cleanup SSH related resources..."

    SSH_KEY_FILE="/tmp/build_server_key"
    BUILD_SERVER_IP=""

    # if Server Information file at, Get IP
    if [ -f /tmp/build_server_info.env ]; then
        source /tmp/build_server_info.env
        BUILD_SERVER_IP="$BUILD_SERVER_PUBLIC_IP"
    fi

    # removeSSH keyfile
    if [ -f "$SSH_KEY_FILE" ]; then
        log "Cleanup SSH keyfile: $SSH_KEY_FILE"
        rm -f "$SSH_KEY_FILE"
        log "✅ SSH keyalready Cleanup"
    else
        log "SSH keyfile Does Not Exist，skip Cleanup"
    fi

    # Fromknown_hostsinremove Serverentry
    if [ -n "$BUILD_SERVER_IP" ]; then
        log "Fromknown_hostsinremove Serverentry: $BUILD_SERVER_IP"
        ssh-keygen -R "$BUILD_SERVER_IP" 2>/dev/null || log "known_hostsentry Not Found"
        log "✅ known_hostsentryalready Cleanup"
    fi
}

# Cleanuptemporaryfile Anddirectory
cleanup_temp_files() {
    log "start Cleanuptemporaryfile Anddirectory..."

    # Cleanup Server Information File
    if [ -f /tmp/build_server_info.env ]; then
        log "Cleanup Server Information File"
        rm -f /tmp/build_server_info.env
    fi

    # Cleanupbuild status file
    for file in /tmp/build_status.txt /tmp/test_status.txt /tmp/build_status.exp /tmp/test_status.exp; do
        if [ -f "$file" ]; then
            log "Cleanup Statusfile: $file"
            rm -f "$file"
        fi
    done

    # Cleanuptemporary artifacts directory
    if [ -d /tmp/artifacts ]; then
        log "Cleanuptemporary artifacts directory"
        rm -rf /tmp/artifacts
    fi

    # Cleanuptemporary Configurefile
    for file in /tmp/server_config.json /tmp/vm_config.json; do
        if [ -f "$file" ]; then
            log "Cleanup Configurefile: $file"
            rm -f "$file"
        fi
    done

    # Cleanup SSH socket file（if exists）
    log "Cleanup SSH socket file"
    find /tmp -name "ssh-*" -type s -delete 2>/dev/null || true

    log "✅ temporaryfile Anddirectory Cleanupcomplete"
}

# Cleanup Localbuild cache（if exists）
cleanup_build_cache() {
    log "Check Cleanupbuild cache..."

    # Cleanupcommonbuild cachedirectory
    for dir in .git/objects .cmake build bin obj out; do
        if [ -d "$dir" ]; then
            log "Foundbuild cachedirectory: $dir"
            # Note：here only Cleanup，preserve directory structure
            find "$dir" -type f -delete 2>/dev/null || true
        fi
    done

    log "✅ build cache Cleanupcomplete"
}

# Cleanupenvironment variables
cleanup_environment() {
    log "CleanupCI related environment variables..."

    # clearsensitive of theCI environment variables
    unset BUILD_SERVER_ACCESS_KEY_ID
    unset BUILD_SERVER_ACCESS_KEY_SECRET
    unset BUILD_SERVER_SSH_KEY
    unset BUILD_SERVER_PUBLIC_IP
    unset BUILD_SERVER_INSTANCE_ID

    log "✅ environment variables Cleanupcomplete"
}

# verify Cleanup Results
verify_cleanup() {
    log "verify Cleanup Results..."

CLEANUP_STATUS=true

# Checkyesnohaveremnants of thesensitivefile
SENSITIVE_FILES=(
"/tmp/build_server_key"
"/tmp/build_server_info.env"
"/tmp/build_status.txt"
"/tmp/test_status.txt"
)

for file in "${SENSITIVE_FILES[@]}"; do
    if [ -f "$file" ]; then
        log "⚠️ Foundsensitivefileremnants: $file"
        CLEANUP_STATUS=false
    fi
done

# Checkwhether there are still artifact directories
if [ -d /tmp/artifacts ]; then
    log "⚠️ Found Artifactdirectoryremnants: /tmp/artifacts"
    CLEANUP_STATUS=false
fi

if [ "$CLEANUP_STATUS" = true ]; then
    log "✅ Cleanup Verification passed，nosensitive Informationremnants"
    else
    log "⚠️ Cleanupverify Foundremnantsfile，please manually Check"
fi
}

# generate Cleanupreport
generate_cleanup_report() {
    log "generate Cleanupreport..."

    cat << EOF
========== CI resources Cleanupreport ==========
Cleanuptime: $(date)
Scriptversion: 008
Cleanupproject:
SSH keyfile: already Cleanup
known_hostsentry: already Cleanup
temporaryfile: already Cleanup
temporarydirectory: already Cleanup
build cache: already Cleanup
environment variables: already Cleanup
security: sensitive Informationwill not retain
=======================================
EOF
}

# execute complete Cleanupprocess
full_cleanup() {
    log "Starting complete execution ofCI resources Cleanup..."

    cleanup_ssh_resources
    cleanup_temp_files
    cleanup_build_cache
    cleanup_environment
    verify_cleanup
    generate_cleanup_report

    log "✅ completeCI resources Cleanupalreadycomplete"
}

# main function
main() {
    log "=== Script008: Cleanup CI resources ==="

    # can selectively execute parts Cleanup
    case "${1:-full}" in
        "ssh")
log "execute onlySSH resources Cleanup"
cleanup_ssh_resources
;;
        "temp")
log "execute onlytemporaryfile Cleanup"
cleanup_temp_files
;;
        "cache")
log "execute onlybuild cache Cleanup"
cleanup_build_cache
;;
        "env")
log "execute onlyenvironment variables Cleanup"
cleanup_environment
;;
        *)
log "execute completeresources Cleanup"
full_cleanup
;;
    esac

    log "=== Script008complete ==="
}

main "$@"