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

# Get server IP
get_server_info() {
    log "Get Server Information..."

    if [ ! -f /tmp/build_server_info.env ]; then
        error "Server Information file does not exist, please run Script 001 first"
    fi

    source /tmp/build_server_info.env

    if [ -z "$BUILD_SERVER_PUBLIC_IP" ]; then
        error "Server IP not set"
    fi

    log "✅ Server IP: $BUILD_SERVER_PUBLIC_IP"
    return "$BUILD_SERVER_PUBLIC_IP"
}

# Configure SSH connection
configure_ssh() {
    log "Configure SSH connection..."

    SSH_KEY="${BUILD_SERVER_SSH_KEY:-}"

    if [ -z "$SSH_KEY" ]; then
        error "BUILD_SERVER_SSH_KEY environment variable not set"
    fi

    # Create SSH key file
    SSH_KEY_FILE="/tmp/build_server_key"
    echo "$SSH_KEY" > "$SSH_KEY_FILE"
    chmod 600 "$SSH_KEY_FILE"

    log "✅ SSH key configured"
}

# Test SSH connection
test_ssh_connection() {
    log "Test SSH connection..."

    local max_attempts=10
    local attempt=1
    local wait_time=10

    while [ $attempt -le $max_attempts ]; do
        log "Attempt SSH connection ($attempt/$max_attempts)..."

        if ssh -i "$SSH_KEY_FILE" -o Strict Host Key Checking=no -o ConnectTimeout=30 root@"$BUILD_SERVER_PUBLIC_IP" "echo 'SSH connection success'" > /dev/null 2>&1; then
            log "✅ SSH connection successful"
            return 0
        fi

        if [ $attempt -lt $max_attempts ]; then
            log "Connection failed, waiting ${wait_time}s before retry..."
            sleep $wait_time
        fi

        attempt=$((attempt + 1))
    done

    error "SSH connection failed, reached maximum attempt count"
}

# Install Docker
install_docker() {
    log "Install Docker on server..."

    ssh -i "$SSH_KEY_FILE" -o Strict Host Key Checking=no root@"$BUILD_SERVER_PUBLIC_IP" bash -s << 'ENDSSH'
set -e
log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1"
}

log "Start Docker install..."

# Check if Docker is already installed
if command -v docker &> /dev/null; then
    log "Docker already installed, updating..."
    apt-get update
    apt-get install -y docker-ce docker-ce-cli containerd.io
else
    log "Install Docker..."
    # Install Docker
    curl -fsSL https://get.docker.com | bash
fi

log "✅ Docker installation complete"
docker version
ENDSSH
    
    log "✅ Docker installation successful"
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
git config --global http.postBuffer 5242880000
git config --global http.lowSpeedLimit 0
git config --global http.timeout 600

# Create Workspace
mkdir -p /workspace
mkdir -p /workspace/logs

ENDSSH
    
    log "✅ Build dependency installation complete"
}

# Main function
main() {
    log "=== Script 002: Prepare CI Build Environment ==="
    get_server_info
    configure_ssh
    test_ssh_connection
    install_docker
    install_build_dependencies
    log "=== Script 002 complete ==="
}

main "$@"