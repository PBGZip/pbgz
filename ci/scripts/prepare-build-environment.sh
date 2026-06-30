#!/bin/bash
# Script 002: Prepare CI Build Environment
# Function: Install necessary tools and dependencies
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

# get_server_IP
get_server_info() {
 log "Get Server Information..."
 
 if [ ! -f /tmp/build_server_info.env ]; then
 error "Server Information File Does Not Exist, Please Run First Script001"
 fi
 
 source /tmp/build_server_info.env
 
 if [ -z "$BUILD_SERVER_PUBLIC_IP" ]; then
 error "server_IP_N ot Set"
 fi
 
 log "✅ serverIP: $BUILD_SERVER_PUBLIC_IP"
 return "$BUILD_SERVER_PUBLIC_IP"
}

# configure_SSH connection
configure_ssh() {
 log "configure_SSH connection..."
 
 SSH_KEY="${BUILD_SERVER_SSH_KEY:-}"
 
 if [ -z "$SSH_KEY" ]; then
 error "BUILD_SERVER_SSH_KEY environment variables Not Set"
 fi
 
 # Create SSH key file
 SSH_KEY_FILE="/tmp/build_server_key"
 echo "$SSH_KEY" > "$SSH_KEY_FILE"
 chmod 600 "$SSH_KEY_FILE"
 
 log "✅ SSH key already Configure"
}

# Test SSH connection
test_ssh_connection() {
 log "Test SSH connection..."
 
 local max_attempts=10
 local attempt=1
 local wait_time=10
 
 while [ $attempt -le $max_attempts ]; do
 log "Aattempt_SSH connection ($attempt/$max_attempts)..."
 
 if ssh -i "$SSH_KEY_FILE" -o Strict Host Key Checking=no -o Connect Timeout=30 root@"$BUILD_SERVER_PUBLIC_IP" "echo 'SSH connection Success'" > /dev/null 2>&1; then
 log "✅ SSH connection Success"
 return 0
 fi
 
 if [ $attempt -lt $max_attempts ]; then
 log "connection Failure，Waiting${wait_time}Retry After Seconds..."
 sleep $wait_time
 fi
 
 attempt=$((attempt + 1))
 done
 
 error "SSH connection Failure，alreadyto Maximum Attempt Count"
}

# Install Docker
install_docker() {
 log "atserver Install Docker..."
 
 ssh -i "$SSH_KEY_FILE" -o Strict Host Key Checking=no root@"$BUILD_SERVER_PUBLIC_IP" bash -s << 'ENDSSH'
set -e
log() {
 echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1"
}

log "start Docker Install..."

# Check Dockeryesnoalready Install
if command -v docker &> /dev/null; then
 log "Dockeralready Install，updatein..."
 apt-get update
 apt-get install -y docker-ce docker-ce-cli containerd.io
else
 log "Install Docker..."
 # Install Docker
 curl -fsSL https://get.docker.com | bash
fi

log "✅ Docker Installcomplete"
docker version
ENDSSH

 log "✅ Docker Install Success"
}

# Install Other Build Dependencies
install_build_dependencies() {
 log "Install Build Dependencies..."
 
 ssh -i "$SSH_KEY_FILE" -o Strict Host Key Checking=no root@"$BUILD_SERVER_PUBLIC_IP" bash -s << 'ENDSSH'
set -e

# Install Build Tools
apt-get update
apt-get install -y \
 git \
 build-essential \
 cmake \
 ninja-build \
 gcc \
 g++ \
 make \
 zip \
 unzip \
 curl \
 wget \
 jq

# Configure Git Network Optimization
git config --global http.post Buffer 5242880000
git config --global http.low Speed Limit 0
git config --global http.timeout 600

# Create Workspace
mkdir -p /workspace
mkdir -p /workspace/logs

ENDSSH

 log "✅ Build Dependency Installationcomplete"
}

# main function
main() {
 log "=== Script002: Pprepare_CI_B uild Environment ==="
 get_server_info
 configure_ssh
 test_ssh_connection
 install_docker
 install_build_dependencies
 log "=== Script002complete ==="
}

main "$@"