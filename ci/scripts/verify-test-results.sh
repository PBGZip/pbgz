#!/bin/bash
# Script 005: Verify Gradle Test Results
# Function: Extract test information from Gradle build results
# Independent debugging: Yes
# Reusable: Yes

set -e

log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1"
}

error() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] ❌ Error: $1"
    exit 1
}

# Get Server Information
get_server_info() {
    if [ ! -f /tmp/build_server_info.env ]; then
        error "Server Information File Does Not Exist, Please Run First Script001"
    fi

    source /tmp/build_server_info.env
    log "✅ get_server_IP: $BUILD_SERVER_PUBLIC_IP"
}

# configure_SSH
configure_ssh() {
    SSH_KEY="${BUILD_SERVER_SSH_KEY:-}"
    SSH_KEY_FILE="/tmp/build_server_key"
    echo "$SSH_KEY" > "$SSH_KEY_FILE"
    chmod 600 "$SSH_KEY_FILE"
    log "✅ SSH key Configurecomplete"
}

# Check Gradle Build Status
check_gradle_status() {
    log "Check Gradle Build Status..."

    ssh -i "$SSH_KEY_FILE" -o Strict Host Key Checking=no root@"$BUILD_SERVER_PUBLIC_IP" bash -s << \'ENDSSH\'
set -e

log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1"
}

# Check Gradle Success
if [ -f "/tmp/gradle_success.txt" ]; then
    log "✅ Gradle Build Successcomplete"
    GRADLE_BUILD_STATUS="SUCCESS"

    # Display Success Information
    cat /tmp/gradle_success.txt

    else
    log "❌ Gradle Buildnot Successcomplete"
    GRADLE_BUILD_STATUS="FAILED"

    # Try To Find Errors Information
    if [ -f "/tmp/gradle_build_status.env" ]; then
        source /tmp/gradle_build_status.env
        log "Build Status: $GRADLE_BUILD_STATUS"
    fi

    # View Latest Gradle Log（ifhave）
    LOG_PATH="/workspace/logs"
    if [ -d "$LOG_PATH" ]; then
        LATEST_LOGS=$(ls -t "$LOG_PATH"/gradle_ci_*.log 2>/dev/null | head -n 1)
        if [ -n "$LATEST_LOGS" ] && [ -f "$LATEST_LOGS" ]; then
            log " of the Gradle Build Log:"
            tail -n 50 "$LATEST_LOGS"
        fi
    fi

    echo "TEST_STATUS=FAILED" > /tmp/test_status_check.env
    echo "TEST_ERROR=Gradle Build Failed" >> /tmp/test_status_check.env

    exit 1
fi

# Check Test Run（because Testalreadyto Gradlein）
log "Check Gradle Test Run..."

LOG_PATH="/workspace/logs"
if [ -d "$LOG_PATH" ]; then
    # Find Containing Test Information of the Log
    for log_file in "$LOG_PATH"/gradle_ci_*.log; do
        if [ -f "$log_file" ]; then
            # Extract Testrelated Information
            if grep -q "test.*passed\|test.*failed\|TEST.*PASSED\|TEST.*FAILED" "$log_file" 2>/dev/null; then
                log "Found Testexecute Information:"
                grep -E "test.*passed|test.*failed|TEST.*PASSED|TEST.*FAILED|✅|❌" "$log_file" | head -n 20
            fi

            # Find Test Taskexecute
            if grep -q "run Tests.*SKIPPED\|run Tests.*UP-TO-DATE\|run Tests.*EXECUTED" "$log_file" 2>/dev/null; then
                log "Test Task Status:"
                grep -E "run Tests.*SKIPPED|run Tests.*UP-TO-DATE|run Tests.*EXECUTED" "$log_file"
            fi
        fi
    done
fi

# Check Test Configure
TEST_ENABLED="${TEST_ENABLED:-false}"
if [ "$TEST_ENABLED" = "true" ]; then
    log "Testalready（TEST_ENABLED=true）"

    # Checkyesnohave Testskip of the Sign
    if [ -d "/workspace/logs" ]; then
        if grep -q "test.*SKIPPED\|skip Test" /workspace/logs/gradle_ci_*.log 2>/dev/null; then
            log "⚠️ to Testskip"
            echo "TEST_STATUS=SKIPPED" > /tmp/test_status_check.env
            else
            log "✅ Testshouldalreadyexecute"
            echo "TEST_STATUS=RUN" > /tmp/test_status_check.env
        fi
    fi
    else
    log "⚠️ Test Not Enabled（TEST_ENABLED=false Not Set）"
    echo "TEST_STATUS=DISABLED" > /tmp/test_status_check.env
fi

log "✅ Gradle Test Status Checkcomplete"
ENDSSH

CHECK_RESULT=$?

if [ $CHECK_RESULT -ne 0 ]; then
    log "⚠️ Gradle Test Status Check Found Problem"
    # notexit，because Mayyes Checkin Has Warnings
fi

log "✅ Test Status Checkcomplete"
}

# Collect Test Summary
collect_test_summary() {
    log "Collect Test Summary..."

ssh -i "$SSH_KEY_FILE" -o Strict Host Key Checking=no root@"$BUILD_SERVER_PUBLIC_IP" "cat /tmp/test_status_check.env" > /tmp/test_status_check.txt 2>/dev/null || log "⚠️ no Get Test Status Check Results"

if [ -f /tmp/test_status_check.txt ]; then
    log "Test Status Check Results:"
    cat /tmp/test_status_check.txt
fi

# Collect Gradle Login of the Test Information
ssh -i "$SSH_KEY_FILE" -o Strict Host Key Checking=no root@"$BUILD_SERVER_PUBLIC_IP" bash -s << 'ENDSSH'
LOG_PATH="/workspace/logs"

if [ -d "$LOG_PATH" ]; then
    # Collect All Testrelated Log
    echo "=== Gradle Testexecute Summary ==="

    for log_file in "$LOG_PATH"/gradle_ci_*.log; do
        if [ -f "$log_file" ]; then
            # Extract Key Information
            grep -A 5 -B 5 "(test|Task.*:test|Building|Complete)" "$log_file" | head -n 30
            break
        fi
    done
fi
ENDSSH

log "✅ Test Summary Collectcomplete"
}

# main function
main() {
    log "=== Script005: Check Gradle Test Results ==="
get_server_info
configure_ssh
check_gradle_status
collect_test_summary
log "=== Script005complete ==="
}

main "$@"