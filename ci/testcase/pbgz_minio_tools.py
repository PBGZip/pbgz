#!/usr/bin/env python3
"""
PBGZ + MinIO Integration Example Script
Demonstrates how to combine PBGZ compression with MinIO object storage
"""
import os
import sys
import argparse
import subprocess
from pathlib import Path

# Add test directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'ci', 'testcase'))

try:
    from minio_storage_adapter import MinIOStorageAdapter, create_minio_adapter_from_env
except ImportError:
    print("Error: minio package not installed. Install with: pip install minio")
    sys.exit(1)


def setup_minio_client(config=None):
    """
    Set up MinIO client
    
    Args:
        config: Configuration dictionary, if None then read from environment variables
    
    Returns:
        MinIOStorageAdapter instance
    
    Configuration Priority:
    1. Explicitly passed config parameter (highest)
    2. Environment variables (medium) - sensitive configuration read from GitHub Secrets
    3. ci-config.properties (lowest) - non-sensitive configuration
    
    Sensitive Configuration (must be read from GitHub Secrets):
    - MINIO_ACCESS_KEY: MinIO access key
    - MINIO_SECRET_KEY: MinIO secret key
    
    Non-sensitive Configuration (can be configured in ci-config.properties):
    - MINIO_ENDPOINT: MinIO service endpoint
    - MINIO_SECURE: Whether to use HTTPS
    - MINIO_REGION: Region setting
    - MINIO_VERIFY_SSL: Whether to verify SSL certificate
    """
    # If explicit configuration is provided, use it directly
    if config and all(key in config for key in ['endpoint', 'access_key', 'secret_key']):
        return MinIOStorageAdapter(
            endpoint=config['endpoint'],
            access_key=config['access_key'],
            secret_key=config['secret_key'],
            secure=config.get('secure', True),
            verify_ssl=config.get('verify_ssl', True)
        )
    
    # Otherwise read configuration from environment variables
    # Sensitive configuration must be read from GitHub Secrets
    access_key = os.getenv('MINIO_ACCESS_KEY')
    secret_key = os.getenv('MINIO_SECRET_KEY')
    
    if not access_key or not secret_key:
        print("Warning: MinIO sensitive credentials not found")
        print("Required environment variables:")
        print("  - MINIO_ACCESS_KEY (from GitHub Secrets)")
        print("  - MINIO_SECRET_KEY (from GitHub Secrets)")
        return None
    
    # Non-sensitive configuration read from ci-config.properties (already loaded into environment variables by workflow)
    endpoint = os.getenv('MINIO_ENDPOINT', 'play.min.io:9000')
    secure = os.getenv('MINIO_SECURE', 'true').lower() == 'true'
    region = os.getenv('MINIO_REGION', 'us-east-1')
    verify_ssl = os.getenv('MINIO_VERIFY_SSL', 'true').lower() == 'true'
    
    return MinIOStorageAdapter(
        endpoint=endpoint,
        access_key=access_key,
        secret_key=secret_key,
        secure=secure,
        region=region,
        verify_ssl=verify_ssl
    )


def compress_and_upload(
    source_file, 
    minio_adapter, 
    bucket_name, 
    object_name=None,
    use_reference=None,
    metadata=None
):
    """
    Compress file and upload to MinIO
    
    Args:
        source_file: Source file path
        minio_adapter: MinIO adapter
        bucket_name: Bucket name
        object_name: Object name (defaults to compressed file name)
        use_reference: Reference genome file path
        metadata: File metadata
    
    Returns:
        Compressed file path
    """
    if not os.path.exists(source_file):
        print(f"Error: Source file not found: {source_file}")
        return None
    
    # Build PBGZ compression command
    compressed_file = f"{source_file}.pbgz"
    command = ["./release-release/bin/pbgz", "compress", source_file, "-o", compressed_file]
    
    # Optional: Use reference genome
    if use_reference:
        command.extend(["-r", use_reference])
    
    print(f"Compressing {source_file}...")
    result = subprocess.run(command, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Compression failed: {result.stderr}")
        return None
    
    print(f"Compression successful: {compressed_file}")
    
    # Upload to MinIO
    if minio_adapter:
        if object_name is None:
            object_name = os.path.basename(compressed_file)
        
        success = minio_adapter.upload_file(
            compressed_file,
            bucket_name,
            object_name,
            metadata=metadata or {}
        )
        
        if success:
            print(f"Successfully uploaded to MinIO: {bucket_name}/{object_name}")
        else:
            print("Failed to upload to MinIO")
    
    return compressed_file


def download_and_decompress(
    minio_adapter,
    bucket_name,
    object_name,
    local_compressed_path,
    local_decompressed_path=None,
    use_reference=None
):
    """
    Download file from MinIO and decompress
    
    Args:
        minio_adapter: MinIO adapter
        bucket_name: Bucket name
        object_name: Object name
        local_compressed_path: Local compressed file path
        local_decompressed_path: Local decompressed file path
        use_reference: Reference genome file path
    
    Returns:
        Decompressed file path
    """
    if not minio_adapter:
        print("Error: MinIO adapter not configured")
        return None
    
    # Download compressed file from MinIO
    print(f"Downloading from MinIO: {bucket_name}/{object_name}")
    downloaded_path = minio_adapter.download_file(
        bucket_name,
        object_name,
        local_compressed_path,
        overwrite=True
    )
    
    if not downloaded_path:
        print("Failed to download from MinIO")
        return None
    
    # Determine decompressed file path
    if local_decompressed_path is None:
        local_decompressed_path = downloaded_path.replace('.pbgz', '')
    
    # Build PBGZ decompression command
    command = ["./release-release/bin/pbgz", "decompress", downloaded_path, "-o", local_decompressed_path]
    
    # Optional: Use reference genome
    if use_reference:
        command.extend(["-r", use_reference])
    
    print(f"Decompressing {downloaded_path}...")
    result = subprocess.run(command, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Decompression failed: {result.stderr}")
        return None
    
    print(f"Decompression successful: {local_decompressed_path}")
    return local_decompressed_path


def main():
    parser = argparse.ArgumentParser(description='PBGZ + MinIO Integration Tool')
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Compress upload subcommand
    compress_parser = subparsers.add_parser('compress-upload', help='Compress file and upload to MinIO')
    compress_parser.add_argument('source_file', help='Source file path')
    compress_parser.add_argument('--bucket', default='pbgz-compressed', help='MinIO bucket name')
    compress_parser.add_argument('--object', help='MinIO object name')
    compress_parser.add_argument('--reference', help='Reference genome file path')
    compress_parser.add_argument('--meta', help='File metadata (JSON format)')
    
    # Download decompress subcommand
    download_parser = subparsers.add_parser('download-decompress', help='Download from MinIO and decompress file')
    download_parser.add_argument('object_name', help='MinIO object name')
    download_parser.add_argument('--bucket', default='pbgz-compressed', help='MinIO bucket name')
    download_parser.add_argument('--output', help='Local output file path')
    download_parser.add_argument('--reference', help='Reference genome file path')
    
    # Configuration options
    parser.add_argument('--minio-endpoint', help='MinIO endpoint')
    parser.add_argument('--minio-access-key', help='MinIO access key')
    parser.add_argument('--minio-secret-key', help='MinIO secret key')
    parser.add_argument('--minio-secure', action='store_true', default=True, help='Use HTTPS')
    parser.add_argument('--minio-no-verify-ssl', action='store_true', help='Disable SSL verification')
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        return 1
    
    # Set up MinIO client
    config = None
    if args.minio_endpoint or os.getenv('MINIO_ENDPOINT'):
        config = {
            'endpoint': args.minio_endpoint or os.getenv('MINIO_ENDPOINT'),
            'access_key': args.minio_access_key or os.getenv('MINIO_ACCESS_KEY'),
            'secret_key': args.minio_secret_key or os.getenv('MINIO_SECRET_KEY'),
            'secure': args.minio_secure,
            'verify_ssl': not args.minio_no_verify_ssl
        }
    
    minio_adapter = setup_minio_client(config)
    
    if not minio_adapter and not args.minio_endpoint:
        print("Warning: MinIO not configured, will operate locally only")
    
    # Execute command
    if args.command == 'compress-upload':
        metadata = {}
        if args.meta:
            try:
                import json
                metadata = json.loads(args.meta)
            except json.JSONDecodeError:
                print("Warning: Invalid metadata JSON, using empty metadata")
        
        result = compress_and_upload(
            args.source_file,
            minio_adapter,
            args.bucket,
            args.object,
            args.reference,
            metadata
        )
        
        return 0 if result else 1
    
    elif args.command == 'download-decompress':
        result = download_and_decompress(
            minio_adapter,
            args.bucket,
            args.object_name,
            f"downloaded_{os.path.basename(args.object_name)}",
            args.output,
            args.reference
        )
        
        return 0 if result else 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())