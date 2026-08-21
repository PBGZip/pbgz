#!/bin/bash
# Script 004: Run Gradle Build
# Function: Use Gradle to run complete build process
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

# Get build server information
get_server_info() {
    if [ ! -f /tmp/build_server_info.env ]; then
        error "Build server information file does not exist, please run script 001 first"
    fi

    source /tmp/build_server_info.env
    log "✅ Obtained build server IP: $BUILD_SERVER_PUBLIC_IP"
}

# Configure SSH
configure_ssh() {
    SSH_KEY="${BUILD_SERVER_SSH_KEY:-}"
    SSH_KEY_FILE="/tmp/build_server_key"
    echo "$SSH_KEY" > "$SSH_KEY_FILE"
    chmod 600 "$SSH_KEY_FILE"
    log "✅ SSH key configuration completed"
}

# Download CI scripts on build server
download_ci_scripts() {
    log "Downloading CI scripts on build server..."

    SOURCE_PATH="/workspace/ci-scripts"

    ssh -i "$SSH_KEY_FILE" -o Strict Host Key Checking=no root@"$BUILD_SERVER_PUBLIC_IP" bash -s << ENDSSH
set -e

# Create CI scripts directory
mkdir -p /workspace/ci-scripts
cd /workspace/ci-scripts

log() {
    echo "[\\$(date +'%Y-%m-%d %H:%M:%S')] \\$1"
}

# Clone CI script repository (using the original CI project repository)
GITHUB_CI_REPO="${GITHUB_CI_REPO:-https://github.com/yourusername/pbgz_ci.git}"
GITHUB_CI_BRANCH="${GITHUB_CI_BRANCH:-main}"

log "Downloading CI script repository: \$GITHUB_CI_REPO (branch: \$GITHUB_CI_BRANCH)"

# Clean up existing CI scripts directory
if [ -d ".git" ]; then
    log "Updating CI scripts..."
    git fetch origin --prune
    git reset --hard origin/\$GITHUB_CI_BRANCH
    else
    log "Cloning CI scripts..."
    git clone -b \$GITHUB_CI_BRANCH \$GITHUB_CI_REPO .
fi

log "✅ CI scripts download completed"
ENDSSH

log "✅ CI scripts download completed"
}

# Configure Gradle environment
configure_gradle() {
    log "Configuring Gradle environment..."

    ssh -i "$SSH_KEY_FILE" -o Strict Host Key Checking=no root@"$BUILD_SERVER_PUBLIC_IP" bash -s << 'ENDSSH'
set -e

CI_SCRIPTS_PATH="/workspace/ci-scripts"
WORKSPACE_PATH="/workspace"

cd "$CI_SCRIPTS_PATH"

log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1"
}

# Check Gradle wrapper
if [ ! -f "gradlew" ]; then
    log "Gradle wrapper does not exist, attempting to create..."

    # Try to create wrapper using system Gradle
    if command -v gradle &> /dev/null; then
        log "Creating wrapper using system Gradle..."
        gradle wrapper --gradle-version 7.6
        else
        log "Install Gradle..."
        # Download Install Gradle
        GRADLE_VERSION="7.6"
        curl -s -o gradle.zip "https://services.gradle.org/distributions/gradle-${GRADLE_VERSION}-bin.zip"
        unzip -o gradle.zip -d /opt/
        rm gradle.zip

        export PATH="/opt/gradle-${GRADLE_VERSION}/bin:$PATH"
        gradle --version

        # Createwrapper
        gradle wrapper --gradle-version $GRADLE_VERSION
    fi
fi

# Settingexecute Permissions
chmod +x gradlew

# verify Gradle
log "verify Gradle Install..."
./gradlew --version || error "Gradlenotyesexecute of the"

log "✅ Gradleenvironment Configurecomplete"
ENDSSH

log "✅ Gradleenvironment Configurecomplete"
}

# execute Gradle Build
run_gradle_build() {
    log "startexecute Gradle Build..."

    SOURCE_PATH="/workspace/ci-scripts"
    WORKSPACE_PATH="/workspace"
    LOG_PATH="/workspace/logs"
    BUILD_TIMESTAMP=$(date +%Y%m%d-%H%M%S)

    # Fromenvironment variables Get Gradle Configure
    GITHUB_REPO="${GITHUB_REPO:-https://github.com/PBGZip/pbgz.git}"
    GITHUB_BRANCH="${GITHUB_BRANCH:-pbgz_v2.2.0}"
    BUILD_SCRIPT_GCC="${BUILD_SCRIPT_GCC:-3rd_party/build.all.gcc.sh}"
    BUILD_SCRIPT_RELEASE="${BUILD_SCRIPT_RELEASE:-build-release.sh}"
    TEST_ENABLED="${TEST_ENABLED:-false}"
    TEST_SCRIPT="${TEST_SCRIPT:-ci/scripts/run_all_tests.sh}"

    log "Gradle Build Configure:"
    log " GitHub Repository: $GITHUB_REPO"
    log " Branch: $GITHUB_BRANCH"
    log " GCC Build Script: $BUILD_SCRIPT_GCC"
    log " Release Build Script: $BUILD_SCRIPT_RELEASE"
    log " Test Enabled: $TEST_ENABLED"
    log " CI Script Path: $SOURCE_PATH"
    log " Workspace Path: $WORKSPACE_PATH"
    log " Log Path: $LOG_PATH"

    ssh -i "$SSH_KEY_FILE" -o Strict Host Key Checking=no root@"$BUILD_SERVER_PUBLIC_IP" bash -s << ENDSSH
set -e

SOURCE_PATH="$SOURCE_PATH"
WORKSPACE_PATH="$WORKSPACE_PATH"
LOG_PATH="$LOG_PATH"
BUILD_TIMESTAMP="$BUILD_TIMESTAMP"
GITHUB_REPO="$GITHUB_REPO"
GITHUB_BRANCH="$GITHUB_BRANCH"
BUILD_SCRIPT_GCC="$BUILD_SCRIPT_GCC"
BUILD_SCRIPT_RELEASE="$BUILD_SCRIPT_RELEASE"
TEST_ENABLED="$TEST_ENABLED"

log() {
    echo "[\\$(date +'%Y-%m-%d %H:%M:%S')] \\$1"
}

# CIS criptdirectory
cd "\$SOURCE_PATH"

# Create Logdirectory
mkdir -p "\$LOG_PATH"

# Configure Gradle Parameters
GRADLE_ARGS=(
"-Pgithub.repo=\$GITHUB_REPO"
"-Pgithub.branch=\$GITHUB_BRANCH"
"-Pbuild.script.gcc=\$BUILD_SCRIPT_GCC"
"-Pbuild.script.release=\$BUILD_SCRIPT_RELEASE"
"-Ptest.enabled=\$TEST_ENABLED"
)

# execute Gradle Build Task
log "execute Gradle Build Task: ci"

# RuncompleteCI process（Download、Build、Test）
if ./gradlew ci "\${GRADLE_ARGS[@]}" --info --stacktrace 2>&1 | tee "\$LOG_PATH/gradle_ci_\$BUILD_TIMESTAMP.log"; then
    log "✅ Gradle Build Successcomplete"
    GRADLE_BUILD_STATUS="SUCCESS"
    else
    log "❌ Gradle Build Failed"
    GRADLE_BUILD_STATUS="FAILED"

    # Display Build Failed of therelated Information
    log "Check Build Log..."
    if [ -f "\$LOG_PATH/gradle_ci_\$BUILD_TIMESTAMP.log" ]; then
        log "Build Failure Login of thefinally50Line:"
        tail -n 50 "\$LOG_PATH/gradle_ci_\$BUILD_TIMESTAMP.log"
    fi

    exit 1
fi

# Build Status
echo "GRADLE_BUILD_STATUS=\$GRADLE_BUILD_STATUS" >> /tmp/gradle_build_status.env

# Check Build Products
BUILD_WORKSPACE="\$WORKSPACE_PATH/workspace"
if [ -d "\$BUILD_WORKSPACE" ]; then
    log "Found Build Workspace: \$BUILD_WORKSPACE"

    artifact_count=\$(find "\$BUILD_WORKSPACE" -type f | wc -l)
    log "Build Productsfile: \$artifact_count"

    if [ \$artifact_count -gt 0 ]; then
        artifact_size=\$(du -sh "\$BUILD_WORKSPACE" | cut -f1)
        log "Build Artifact Total Size: \$artifact_size"
    fi
fi

log "✅ Gradle Build Taskexecutecomplete"
echo "✅ Gradle CI process Successcomplete" > /tmp/gradle_success.txt

# Run tests after build if enabled
if [ "\$TEST_ENABLED" = "true" ]; then
    log "Test enabled, running unit tests and black-box test cases..."

    # run_all_tests.sh lives in the CI repository under ci/scripts/
    TEST_SCRIPT_PATH="\$SOURCE_PATH/ci/scripts/run_all_tests.sh"
    if [ -f "\$TEST_SCRIPT_PATH" ]; then
        # Determine actual build test binary path (Gradle downloads source to workspace/pbgz)
        # If test binary not found under default path, probe common locations
        PBGZ_TEST_BIN="\$WORKSPACE_PATH/workspace/pbgz/build-x86_64-gcc-release/test/pbgz_test"
        if [ ! -f "\$PBGZ_TEST_BIN" ]; then
            PBGZ_TEST_BIN="\$(find \$WORKSPACE_PATH/workspace -name 'pbgz_test' -type f 2>/dev/null | head -n 1)"
        fi

        if [ -n "\$PBGZ_TEST_BIN" ] && [ -f "\$PBGZ_TEST_BIN" ]; then
            log "Running tests with C++ test binary: \$PBGZ_TEST_BIN"
            bash "\$TEST_SCRIPT_PATH" --pbgz-test-path="\$PBGZ_TEST_BIN" \
                --python-test-path="\$SOURCE_PATH/ci/run_all_testcase.py" \
                > "\$LOG_PATH/test_all_\$BUILD_TIMESTAMP.log" 2>&1
            TEST_RESULT=\$?
            if [ \$TEST_RESULT -ne 0 ]; then
                log "❌ Tests FAILED (exit code \$TEST_RESULT)"
                tail -n 50 "\$LOG_PATH/test_all_\$BUILD_TIMESTAMP.log"
                exit 1
            else
                log "✅ Tests PASSED"
            fi
        else
            log "❌ Test binary not found under workspace, cannot run tests"
            exit 1
        fi
    else
        log "❌ Test script not found: \$TEST_SCRIPT_PATH"
        exit 1
    fi
else
    log "Tests disabled (TEST_ENABLED=false), skipping tests"
fi
ENDSSH

GRADLE_RESULT=$?

if [ $GRADLE_RESULT -ne 0 ]; then
    error "Gradle Build Taskexecute Failure"
fi

log "✅ Gradle Build Taskexecutecomplete"
}

# Collect Gradle Build Status
collect_gradle_status() {
    log "Collect Gradle Build Status..."

ssh -i "$SSH_KEY_FILE" -o Strict Host Key Checking=no root@"$BUILD_SERVER_PUBLIC_IP" "cat /tmp/gradle_build_status.env" > /tmp/gradle_build_status.txt 2>/dev/null || log "⚠️ no Get Gradle Build Status"

if [ -f /tmp/gradle_build_status.txt ]; then
    log "Gradle Build Status:"
    cat /tmp/gradle_build_status.txt
fi
}

# verify Gradle Install
verify_gradle_installation() {
    log "verify Gradle Install..."

    ssh -i "$SSH_KEY_FILE" -o Strict Host Key Checking=no root@"$BUILD_SERVER_PUBLIC_IP" bash -s << 'ENDSSH'
# Preventive Check Gradleyesno
cd /workspace/ci-scripts 2>/dev/null || (mkdir -p /workspace/ci-scripts && cd /workspace/ci-scripts)

if [ -f "./gradlew" ]; then
    echo "✅ Gradle wrapper already Install"
    ./gradlew --version || echo "⚠️ Gradle wrapper Test Failure"
else
    echo "⚠️ Gradle wrapper Not Found"
fi
ENDSSH

    log "✅ Gradle Installverifycomplete"
}

# main function
main() {
    log "=== Script004: Run Gradle Build ==="
get_server_info
configure_ssh
verify_gradle_installation
download_ci_scripts
configure_gradle
run_gradle_build
collect_gradle_status
log "=== Script004complete ==="
}

main "$@"