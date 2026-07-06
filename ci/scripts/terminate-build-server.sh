#!/bin/bash
# Script007: Terminate Server
# Function: Terminate Server Instance
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

# Get Server Information
get_server_info() {
    log "Get Server Information..."

    if [ ! -f /tmp/build_server_info.env ]; then
        error "Server Information File Does Not Exist, Please Run First Script001"
    fi

    source /tmp/build_server_info.env

    if [ -z "$BUILD_SERVER_INSTANCE_ID" ]; then
        error "server instance ID Not Set"
    fi

    log "instance Information:"
    log " ID: $BUILD_SERVER_INSTANCE_ID"
    log " IP: $BUILD_SERVER_PUBLIC_IP"
    log " Status: $BUILD_SERVER_STATUS"

    log "✅ server Information Getcomplete"
}

# Configure SSH key
configure_ssh() {
    log "Configure SSH key..."

    SSH_KEY="${BUILD_SERVER_SSH_KEY:-}"

    if [ -z "$SSH_KEY" ]; then
        error "BUILD_SERVER_SSH_KEY environment variable Not Set"
    fi

    # Create SSH key file
    SSH_KEY_FILE="/tmp/build_server_key"
    echo "$SSH_KEY" > "$SSH_KEY_FILE"
    chmod 600 "$SSH_KEY_FILE"

    log "✅ SSH key configuration completed"
}

# verify Server Status
verify_server_status() {
    log "verify Server Status..."

    ssh -i "$SSH_KEY_FILE" -o Strict Host Key Checking=no root@"$BUILD_SERVER_PUBLIC_IP" "echo 'Server Access'" 2>/dev/null

    if [ $? -eq 0 ]; then
        log "✅ Server Access，Run"
    else
        log "⚠️ Serveralreadyno Access，Mayalready Terminateenvironment Exception"
    fi
}

# Terminate Server Instance
terminate_server() {
    log "start Terminate Server Instance..."

    INSTANCE_ID="$BUILD_SERVER_INSTANCE_ID"

    log "Currently Terminate instance: $INSTANCE_ID"

    # Checkinstance Current Status
    log "Checkinstance Current Status..."
    CURRENT_STATUS=$(aliyun ecs Describe Instances --Instance Ids $INSTANCE_ID --output json | jq -r '.Instances[0].Status // empty' || echo "")

    if [ "$CURRENT_STATUS" = "Stopped" ]; then
        log "instancealreadyyes Stopped Status"
        return 0
    elif [ "$CURRENT_STATUS" = "Terminated" ]; then
        log "instancealreadyyes Terminated Status"
        return 0
    elif [ -z "$CURRENT_STATUS" ]; then
        error "no Getinstance Status，instance Mayalready Does Not Exist"
    else
        log "instance Current Status: $CURRENT_STATUS"
    fi

# execute Terminate
log "execute Terminate..."
TERMINATE_RESULT=$(aliyun ecs Stop Instance --Instance Id $INSTANCE_ID --Force Stop true --Confirm Terminate true --output json 2>&1)

if [ $? -eq 0 ]; then
    log "✅ Terminateexecute Success"
else
    log "⚠️ Terminateexecute Failure（May Needconfirm）"
    log "Error Information: $TERMINATE_RESULT"
fi

# verify Terminate Status
log "verify Terminate Status..."
MAX_WAIT=180
WAIT_INTERVAL=10
elapsed=0

while [ $elapsed -lt $MAX_WAIT ]; do
    STATUS=$(aliyun ecs Describe Instances --Instance Ids $INSTANCE_ID --output json | jq -r '.Instances[0].Status // empty' || echo "")

    case "$STATUS" in
    "Stopped")
        log "✅ instancealready Successfully Stopped"
        return 0
        ;;
    "Stopping")
        log "instance Currently Stoppingin... ($elapsed/$MAX_WAIT seconds)"
        ;;
    "Terminated")
        log "✅ instancealready Terminate"
        return 0
        ;;
    "")
        log "⚠️ no Getinstance Status"
        ;;
    *)
        log "instance Status: $STATUS"
        ;;
esac

    if [ $elapsed -lt $MAX_WAIT ]; then
        sleep $WAIT_INTERVAL
        elapsed=$((elapsed + WAIT_INTERVAL))
    fi
done

    if [ "$STATUS" = "Stopped" ] || [ "$STATUS" = "Terminated" ]; then
        log "✅ instancealready Success Terminate"
    else
        log "⚠️ Waiting Timeout，Terminate Mayat Line"
    fi
}

# ultimately Cleanup
final_cleanup() {
    log "executeultimately Cleanup..."

    # Cleanup Localtemporaryfile
    rm -f /tmp/build_server_info.env
    rm -f /tmp/build_status.txt
    rm -f /tmp/test_status.txt

    log "✅ ultimately Cleanupcomplete"
}

# Record Terminate Results
record_termination_result() {
    log "Record Terminate Results..."

    cat << EOF
========== server Terminate Results ==========
instance_ID: $BUILD_SERVER_INSTANCE_ID
Public NetworkIP: $BUILD_SERVER_PUBLIC_IP
Status: $BUILD_SERVER_STATUS
Terminatetime: $(date)
ultimately Status: $STATUS
=======================================
EOF
}

# main function
main() {
    log "=== Script007: Terminate Server ==="
get_server_info
configure_ssh
verify_server_status
terminate_server
final_cleanup
record_termination_result
log "=== Script007complete ==="
}

main "$@"