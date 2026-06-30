#!/bin/bash
# Script006: Collect Build Products
# Function: Collect Build Log And Binary Artifacts
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

# Collect Build Log
collect_build_logs() {
 log "start Collect Build Log..."
 
 SOURCE_PATH="/workspace/source-code"
 LOG_PATH="/workspace/logs"
 LOCAL_ARTIFACTS_PATH="/tmp/artifacts"
 TIMESTAMP=$(date +%Y%m%d-%H%M%S)
 
 # Create Local Artifactdirectory
 mkdir -p "$LOCAL_ARTIFACTS_PATH"
 
 ssh -i "$SSH_KEY_FILE" -o Strict Host Key Checking=no root@"$BUILD_SERVER_PUBLIC_IP" bash -s << ENDSSH
SOURCE_PATH="$SOURCE_PATH"
LOG_PATH="$LOG_PATH"
TIMESTAMP="$TIMESTAMP"

log() {
 echo "[\$(date +'%Y-%m-%d %H:%M:%S')] \$1"
}

# Check Logdirectoryyesnoat
if [ ! -d "\$LOG_PATH" ]; then
 log "❌ Logdirectory Does Not Exist"
 exit 1
fi

# Create Packagedirectory
PACKAGING_DIR="/tmp/artifact_packaging"
mkdir -p "\$PACKAGING_DIR"

# Collect All Build Logs
log "Collect Build Log..."
if [ -f "\$LOG_PATH/build.log" ]; then
 cp \$LOG_PATH/build.log "\$PACKAGING_DIR/"
fi

# time Collect Log
find \$LOG_PATH -name "build-\$TIMESTAMP.log" -type f -exec cp {} \$PACKAGING_DIR/ \; 2>/dev/null
find \$LOG_PATH -name "build-\$TIMESTAMP-error.log" -type f -exec cp {} \$PACKAGING_DIR/ \; 2>/dev/null

# Collect Test Log
find \$LOG_PATH -name "test-\$TIMESTAMP.log" -type f -exec cp {} \$PACKAGING_DIR/ \; 2>/dev/null
find \$LOG_PATH -name "test-\$TIMESTAMP-error.log" -type f -exec cp {} \$PACKAGING_DIR/ \; 2>/dev/null

# Collect of the Logfile（ifhavetime of the）
if [ -z "\$(ls -A \$PACKAGING_DIR/*.log 2>/dev/null)" ]; then
 log "Not Foundtime of the Log，Collect Latest Log..."
 find \$LOG_PATH -name "*.log" -type f -mtime -1 -exec cp {} \$PACKAGING_DIR/ \; 2>/dev/null
fi

# Create Log Compression Package
log "Create Log Compression Package..."
BUILD_LOGS_ARCHIVE="ci_logs_\$TIMESTAMP.tar.gz"
tar -czf "\$BUILD_LOGS_ARCHIVE" -C \$PACKAGING_DIR . 2>/dev/null

log "✅ Log Collectcomplete: \$BUILD_LOGS_ARCHIVE"
ENDSSH

 log "✅ Build Log Collectcomplete"
}

# Collect Binary Artifacts
collect_binary_artifacts() {
 log "start Collect Binary Artifacts..."
 
 SOURCE_PATH="/workspace/source-code"
 TIMESTAMP=$(date +%Y%m%d-%H%M%S)
 
 ssh -i "$SSH_KEY_FILE" -o Strict Host Key Checking=no root@"$BUILD_SERVER_PUBLIC_IP" bash -s << ENDSSH
SOURCE_PATH="$SOURCE_PATH"
TIMESTAMP="$TIMESTAMP"

log() {
 echo "[\$(date +'%Y-%m-%d %H:%M:%S')] \$1"
}

# Create Packagedirectory
PACKAGING_DIR="/tmp/artifact_packaging"
mkdir -p "\$PACKAGING_DIR"

# Find Build Outputdirectory
BUILD_OUTPUT="build"
if [ -d "\$BUILD_OUTPUT" ]; then
 log "Found Build Outputdirectory: \$BUILD_OUTPUT"
 
 # Collect All Build Artifacts
 log "Collect Build Products..."
 find \$BUILD_OUTPUT -type f \
 \( -name "*.a" -o -name "*.so" -o -name "*.lib" -o -name "*.dll" \) \
 -exec cp {} \$PACKAGING_DIR/ \; 2>/dev/null
 
 # Collectexecutefile
 find \$BUILD_OUTPUT -type f -executable -exec cp {} \$PACKAGING_DIR/ \; 2>/dev/null
 
 # Collect Other Importantfile
 find \$BUILD_OUTPUT -type f \
 \( -name "*.h" -o -name "*.hpp" -o -name "*.c" -o -name "*.cpp" \) \
 -exec cp --parents {} \$PACKAGING_DIR/ \; 2>/dev/null
fi

# Checkyesnohave Othercommon Artifactdirectory
for DIR in "bin" "lib" "dist" "out" "release" "debug"; do
 if [ -d "\$DIR" ]; then
 log "Found Artifactdirectory: \$DIR"
 tar -cf - \$DIR 2>/dev/null | (cd \$PACKAGING_DIR && tar -xf -)
 fi
done

# generate Artifact Manifest
log "generate Artifact Manifest..."
if [ -n "\$(ls -A \$PACKAGING_DIR 2>/dev/null)" ]; then
 find \$PACKAGING_DIR -type f > \$PACKAGING_DIR/files_manifest.txt
 DU_OUTPUT=\$(du -sh \$PACKAGING_DIR | cut -f1)
 log "Artifact Total Size: \$DU_OUTPUT"
 
 # Create Artifact Compression Package
 log "Create Artifact Compression Package..."
 BUILD_ARTIFACTS_ARCHIVE="build_artifacts_\$TIMESTAMP.tar.gz"
 tar -czf "\$BUILD_ARTIFACTS_ARCHIVE" -C /tmp/artifact_packaging . 2>/dev/null
 
 log "✅ Artifact Collectcomplete: \$BUILD_ARTIFACTS_ARCHIVE"
else
 log "⚠️ Not Found Build Products"
 BUILD_ARTIFACTS_ARCHIVE=""
fi
ENDSSH

 log "✅ Binary Artifacts Collectcomplete"
}

# Download Artifactto Local
download_artifacts_to_local() {
 log "Download Artifactto Local..."
 
 LOCAL_ARTIFACTS_PATH="/tmp/artifacts"
 mkdir -p "$LOCAL_ARTIFACTS_PATH"
 
 ssh -i "$SSH_KEY_FILE" -o Strict Host Key Checking=no root@"$BUILD_SERVER_PUBLIC_IP" bash -s << ENDSSH
# List To Download of the file
if [ -f "/tmp/artifact_packaging/ci_logs_*.tar.gz" ]; then
 ls -lh /tmp/artifact_packaging/ci_logs_*.tar.gz
fi

if [ -f "/tmp/artifact_packaging/build_artifacts_*.tar.gz" ]; then
 ls -lh /tmp/artifact_packaging/build_artifacts_*.tar.gz
fi
ENDSSH
 
 # UseSCPD ownloadfile
 scp -i "$SSH_KEY_FILE" -o Strict Host Key Checking=no \
 root@"$BUILD_SERVER_PUBLIC_IP":/tmp/artifact_packaging/*.tar.gz \
 "$LOCAL_ARTIFACTS_PATH/" 2>/dev/null || log "⚠️ part Artifact Download Failure"
 
 # Download AnyZIP file（ifhave）
 scp -i "$SSH_KEY_FILE" -o Strict Host Key Checking=no \
 root@"$BUILD_SERVER_PUBLIC_IP":/tmp/artifact_packaging/*.zip \
 "$LOCAL_ARTIFACTS_PATH/" 2>/dev/null || true
 
 log "✅ Artifact Downloadcomplete"
 ls -lh "$LOCAL_ARTIFACTS_PATH/"
}

# Upload To Git Hub Actions（ifat Git Hub Actionsenvironmentin Run）
upload_to_github_actions() {
 log "Checkyesnoneed Upload To Git Hub Actions..."
 
 if [ -n "$GITHUB_JOB" ] && [ -n "$GITHUB_RUN_ID" ]; then
 log "at Git Hub Actionsenvironmentin，Upload Artifact..."
 
 LOCAL_ARTIFACTS_PATH="/tmp/artifacts"
 TIMESTAMP=$(date +%Y%m%d-%H%M%S)
 
 # UploadCIL og
 for file in "$LOCAL_ARTIFACTS_PATH"/ci_logs_*.tar.gz; do
 if [ -f "$file" ]; then
 echo "Upload Log Artifacts: $file"
 fi
 done
 
 # Upload Build Products
 for file in "$LOCAL_ARTIFACTS_PATH"/build_artifacts_*.tar.gz; do
 if [ -f "$file" ]; then
 echo "Upload Build Products: $file"
 fi
 done
 
 log "✅ Artifactalready Prepared Git Hub Actions Upload"
 else
 log "notat Git Hub Actionsenvironmentin，skip Upload"
 fi
}

# Analyze Artifacts Results
analyze_artifacts() {
 log "Analyze Artifacts Results..."
 
 LOCAL_ARTIFACTS_PATH="/tmp/artifacts"
 
 if [ -d "$LOCAL_ARTIFACTS_PATH" ] && [ -n "$(ls -A $LOCAL_ARTIFACTS_PATH)" ]; then
 log "📊 Artifact Analysis:"
 
 # Check Logfile
 for file in "$LOCAL_ARTIFACTS_PATH"/ci_logs_*.tar.gz; do
 if [ -f "$file" ]; then
 SIZE=$(du -h "$file" | cut -f1)
 log " CIL og: $(basename "$file") - $SIZE"
 fi
 done
 
 # Check Build Products
 for file in "$LOCAL_ARTIFACTS_PATH"/build_artifacts_*.tar.gz; do
 if [ -f "$file" ]; then
 SIZE=$(du -h "$file" | cut -f1)
 log " Build Products: $(basename "$file") - $SIZE"
 fi
 done
 
 TOTAL_SIZE=$(du -sh "$LOCAL_ARTIFACTS_PATH" | cut -f1)
 log " : $TOTAL_SIZE"
 
 log "✅ Artifact Analysiscomplete"
 else
 log "⚠️ Not Found Any Artifactfile"
 fi
}

# main function
main() {
 log "=== Script006: Collect Build Products ==="
 get_server_info
 configure_ssh
 collect_build_logs
 collect_binary_artifacts
 download_artifacts_to_local
 upload_to_github_actions
 analyze_artifacts
 log "=== Script006complete ==="
}

main "$@"