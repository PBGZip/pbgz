#!/bin/bash
# MinIO Storage Integration Script for CI/CD
# Provides functionality for uploading test results, backing up build artifacts, etc.

set -e

# Load environment variables
sourceci_env() {
      if [ -f "ci/ci-config.properties" ]; then
          while IFS='=' read -r key value; do
              [[ $key =~ ^#.*$ ]] && continue
              [[ -z $key ]] && continue
              value=$(echo "$value" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
              export "$key"="$value"
          done < ci/ci-config.properties
      fi
}

# Initialize MinIO adapter
init_minio() {
          echo "🗄️ Initializing MinIO storage..."

          # Check required sensitive configuration
          if [ -z "$MINIO_ACCESS_KEY" ] || [ -z "$MINIO_SECRET_KEY" ]; then
                  echo "❌ Error: MinIO credentials not found"
                  echo "Required: MINIO_ACCESS_KEY and MINIO_SECRET_KEY from GitHub Secrets"
                  return 1
          fi

          # Check loaded configuration
          echo "✅ MinIO Configuration:"
          echo "   Endpoint: ${MINIO_ENDPOINT:-not set}"
          echo "   Secure: ${MINIO_SECURE:-false}"
          echo "   Region: ${MINIO_REGION:-not set}"

          return 0
}

# Upload test results to MinIO
upload_test_results() {
      local test_file="$1"
      local object_name="$2"
      
      if [ ! -f "$test_file" ]; then
          echo "❌ Error: Test file not found: $test_file"
          return 1
      fi
      
      echo "📤 Uploading test results to MinIO..."
      
      # Use Python MinIO tool for upload
      if python ci/testcase/pbgz_minio_tools.py compress-upload \
          "$test_file" \
          --bucket "${MINIO_RESULTS_BUCKET:-pbgz-results}" \
          --object "$object_name" \
          --meta "{\"workflow\":\"${WORKFLOW_TYPE:-unknown}\",\"run_id\":\"${GITHUB_RUN_ID:-$(date +%s)}\"}"; then
          echo "✅ Test results uploaded successfully: $object_name"
          return 0
      else
          echo "⚠️ Warning: Failed to upload test results"
          return 1
      fi
}

# Upload build artifacts to MinIO
upload_build_artifacts() {
      local artifact_pattern="$1"
      
      echo "📦 Uploading build artifacts to MinIO..."
      
      for artifact in $artifact_pattern; do
          if [ -f "$artifact" ]; then
              local artifact_name=$(basename "$artifact")
              local timestamp=$(date +%Y%m%d-%H%M%S)
              local object_name="artifacts/${WORKFLOW_TYPE:-unknown}/${timestamp}/${artifact_name}"
              
              if python ci/testcase/pbgz_minio_tools.py compress-upload \
                  "$artifact" \
                  --bucket "${MINIO_COMPRESSED_BUCKET:-pbgz-compressed}" \
                  --object "$object_name"; then
                  echo "✅ Artifact uploaded: $artifact_name"
              else
                  echo "⚠️ Warning: Failed to upload $artifact_name"
              fi
          fi
      done
      
      return 0
}

# Download test data from MinIO
download_test_data() {
      local object_name="$1"
      local local_file="$2"
      
      if [ -z "$object_name" ] || [ -z "$local_file" ]; then
          echo "❌ Error: Missing arguments"
          echo "Usage: download_test_data \<object_name\> \<local_file\>"
          return 1
      fi
      
      echo "📥 Downloading test data from MinIO..."
      
      python ci/testcase/pbgz_minio_tools.py download-decompress \
          "$object_name" \
          --bucket "${MINIO_SOURCE_BUCKET:-pbgz-source}" \
          --output "$local_file" || true
      
      if [ -f "$local_file" ]; then
          echo "✅ Test data downloaded successfully: $local_file"
          return 0
      else
          echo "❌ Error: Failed to download test data"
          return 1
      fi
}

# Backup important files to MinIO
backup_to_minio() {
      local local_file="$1"
      local backup_name="$2"
      
      if [ ! -f "$local_file" ]; then
          echo "❌ Error: File not found: $local_file"
          return 1
      fi
      
      echo "💾 Backing up file to MinIO..."
      
      local timestamp=$(date +%Y%m%d-%H%M%S)
      local object_name="backups/${WORKFLOW_TYPE:-unknown}/${backup_name}_${timestamp}"
      
      python ci/testcase/pbgz_minio_tools.py compress-upload \
          "$local_file" \
          --bucket "${MINIO_BACKUP_BUCKET:-pbgz-backups}" \
          --object "$object_name" || true
      
      if [ $? -eq 0 ]; then
          echo "✅ Backup completed: $object_name"
          return 0
      else
          echo "⚠️ Warning: Backup failed"
          return 1
      fi
}

# Clean up old files from MinIO
cleanup_old_files() {
      local bucket="$1"
      local keep_days="${2:-7}"
      
      echo "🧹 Cleaning up old files from MinIO bucket: $bucket (keeping ${keep_days} days)"
      
      # Cleanup logic can be added here
      # Currently just an example, actual implementation needs to be based on specific requirements
      
      echo "⚠️ Cleanup not implemented yet"
      return 0
}

# List files in specific bucket
list_minio_files() {
      local bucket="$1"
      local prefix="$2"
      
      echo "📋 Listing files in MinIO bucket: $bucket"
      echo "   Prefix: $prefix"
      
      # Use Python SDK to list files
      local bucket_name="${bucket:-pbgz-results}"
      local prefix_name="${prefix:-}"
      BUCKET="$bucket_name" PREFIX="$prefix_name" python3 << 'PYEOF'
import sys
import os
sys.path.insert(0, 'ci/testcase')
try:
    from minio_storage_adapter import create_minio_adapter_from_env

    adapter = create_minio_adapter_from_env()
    if adapter:
        files = adapter.list_files(os.environ.get('BUCKET', 'pbgz-results'), os.environ.get('PREFIX', ''))
        print(f"Found {len(files)} files:")
        for file in files:
            print(f"  - {file}")
    else:
        print("MinIO adapter not configured")
except Exception as e:
    print(f"Error: {e}")
    sys.exit(1)
PYEOF
}

# Main function
main() {
      sourceci_env
      
      local command="$1"
      shift
      
      case "$command" in
          init)
              init_minio
              ;;
          upload-test)
              upload_test_results "$@"
              ;;
          upload-artifacts)
              upload_build_artifacts "$@"
              ;;
          download-data)
              download_test_data "$@"
              ;;
          backup)
              backup_to_minio "$@"
              ;;
          cleanup)
              cleanup_old_files "$@"
              ;;
          list)
              list_minio_files "$@"
              ;;
          *)
              echo "MinIO Storage Integration Script"
              echo "Usage: $0 {init|upload-test|upload-artifacts|download-data|backup|cleanup|list} [args...]"
              echo ""
              echo "Commands:"
              echo "  init                    - Initialize MinIO connection"
              echo "  upload-test <file> <name>  - Upload test results"
              echo "  upload-artifacts <pattern>  - Upload build artifacts"
              echo "  download-data <object> <local> - Download test data"
              echo "  backup <file> <name>     - Backup file to MinIO"
              echo "  cleanup <bucket> <days>  - Clean up old files"
              echo "  list <bucket> <prefix>   - List files in bucket"
              exit 1
              ;;
      esac
}

# Execute main function
main "$@"