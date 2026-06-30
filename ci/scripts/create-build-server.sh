#!/bin/bash
# Script 001: Create CI Server
# Function: Dynamically create CI build server
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

# Check dependencies
check_dependencies() {
 log "Checking build server dependencies..."
 
 if ! command -v aliyun &> /dev/null; then
 log "Build server CLI not installed, preparing..."
 if command -v curl &> /dev/null; then
 log "Installing build server CLI..."
 fi
 log "✅ Build server CLI ready"
 fi
 
 if ! aliyun configure list &> /dev/null; then
 log "server_CLI not Configure，Use Access Credentials Configure..."
 fi
}

# Createserverinstance
create_server() {
 log "CreateCI server..."
 
 SERVER_NAME="${SERVER_NAME:-ci-$(date +%Y%m%d-%H%M%S)}"
 REGION="${REGION:-default}"
 INSTANCE_TYPE="${INSTANCE_TYPE:-}"
 IMAGE_ID="${IMAGE_ID:-ubuntu-20_04}"
 SECURITY_GROUP_ID="${SECURITY_GROUP_ID:-}"
 KEY_PAIR_NAME="${KEY_PAIR_NAME:-ci-key}"
 
 log "Server Configure:"
 log " : $SERVER_NAME"
 log " Region: $REGION"
 log " instance Type: $INSTANCE_TYPE"
 log " Image: $IMAGE_ID"
 log " Security: $SECURITY_GROUP_ID"
 log " key: $KEY_PAIR_NAME"
 
 log "🔍 environment variables Check:"
 log " ACCESS_KEY_ID: [${#BUILD_SERVER_ACCESS_KEY_ID} chars]"
 log " ACCESS_KEY_SECRET: [${#BUILD_SERVER_ACCESS_KEY_SECRET} chars]"
 log " REGION: $REGION"
 
 if [ -z "${BUILD_SERVER_ACCESS_KEY_ID}" ]; then
 error "BUILD_SERVER_ACCESS_KEY_ID environment variables Not Set"
 fi
 
 if [ -z "${BUILD_SERVER_ACCESS_KEY_SECRET}" ]; then
 error "BUILD_SERVER_ACCESS_KEY_SECRET environment variables Not Set"
 fi
 
 log "server Create..."
 
 cat > /tmp/server_config.json <<EOF
{
 "Region Id": "$REGION",
 "Image Id": "$IMAGE_ID",
 "Instance Type": "$INSTANCE_TYPE",
 "Security Group Id": "$SECURITY_GROUP_ID",
 "Instance Name": "$SERVER_NAME",
 "Internet Max Bandwidth Out": 50,
 "System Disk": {
 "Category": "cloud_ssd",
 "Size": 40
 },
 "Password": null,
 "Allocate Public Ip": true
}
EOF
 
 log "Configurefile Createcomplete"
 
 export BUILD_SERVER_ACCESS_KEY_ID="${BUILD_SERVER_ACCESS_KEY_ID}"
 export BUILD_SERVER_ACCESS_KEY_SECRET="${BUILD_SERVER_ACCESS_KEY_SECRET}"
 export BUILD_SERVER_DEFAULT_REGION="$REGION"
 
 log "🔍 CLIC onfigure:"
 log " ACCESS_KEY_ID: ${BUILD_SERVER_ACCESS_KEY_ID:0:8}****"
 log " DEFAULT_REGION: $BUILD_SERVER_DEFAULT_REGION"
 
 log "Currently Creatingserver（May Need）..."
 INSTANCE_ID=$(aliyun ecs Run Instance --cli-config /tmp/server_config.json --output json | jq -r '.Instance Id // empty' || echo "")
 
 if [ "$INSTANCE_ID" = "empty" ] || [ -z "$INSTANCE_ID" ]; then
 error "Createserver Failure， of the Instance Id"
 fi
 
 log "✅ server Create Success"
 log " instance_ID: $INSTANCE_ID"
 log " Server: $SERVER_NAME"
 
 echo "BUILD_SERVER_INSTANCE_ID=$INSTANCE_ID" > /tmp/build_server_info.env
 
 log "Wait For Server Start And Get Public NetworkIP..."
 
 MAX_WAIT=300
 WAIT_INTERVAL=10
 elapsed=0
 
 while [ $elapsed -lt $MAX_WAIT ]; do
 SERVER_INFO=$(aliyun ecs Describe Instances --Instance Ids $INSTANCE_ID --output json 2>/dev/null || echo "")
 
 if [ -n "$SERVER_INFO" ]; then
 PUBLIC_IP=$(echo "$SERVER_INFO" | jq -r '.Instances[0].Public Ip Address.Ip Address // empty' || "")
 STATUS=$(echo "$SERVER_INFO" | jq -r '.Instances[0].Status // empty' || "")
 
 if [ "$PUBLIC_IP" != "empty" ] && [ -n "$PUBLIC_IP" ]; then
 log "✅ Serveralreadystart Public NetworkIP: $PUBLIC_IP"
 log " Status: $STATUS"
 echo "BUILD_SERVER_PUBLIC_IP=$PUBLIC_IP" >> /tmp/build_server_info.env
 echo "BUILD_SERVER_STATUS=$STATUS" >> /tmp/build_server_info.env
 break
 fi
 fi
 
 log "Waitingin... ($elapsed/$MAX_WAIT seconds)"
 sleep $WAIT_INTERVAL
 elapsed=$((elapsed + WAIT_INTERVAL))
 done
 
 if ! grep -q "BUILD_SERVER_PUBLIC_IP" /tmp/build_server_info.env; then
 error "Serverstartnot Public NetworkIP"
 fi
 
 log "✅ serveralready Preparethen"
}

# main function
main() {
 log "=== Script001: Createserver ==="
 check_dependencies
 create_server
 log "=== Script001complete ==="
}

main "$@"