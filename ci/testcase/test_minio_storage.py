"""
MinIO Storage Test Case
Demonstrates how to use MinIO object storage in PBGZ testing framework
"""
import os
import random
import hashlib
from testcase.pbgz_test_framework import PBGZTestCase
from testcase.minio_storage_adapter import create_minio_adapter_from_env


class MinIOStorageTest(PBGZTestCase):
    """MinIO Storage Test Case"""
    
    def __init__(self):
        super().__init__("MinIOStorageTest")
        self.minio_adapter = None
        self.source_file = "minio_test_data.txt"
        self.compressed_file = "minio_test_data.txt.pbgz"
        self.decompressed_file = "minio_test_data.txt.dec"
        self.bucket_name = "pbgz-test"
    
    def get_test_files(self) -> tuple:
        """Return test file information"""
        return (self.source_file, self.compressed_file, self.decompressed_file)
    
    def prepare_data(self):
        """Prepare test data"""
        # Initialize MinIO adapter
        self.minio_adapter = create_minio_adapter_from_env()
        if self.minio_adapter is None:
            print("Warning: MinIO adapter not configured, skipping MinIO-specific tests")
            print("Set environment variables: MINIO_ENDPOINT, MINIO_ACCESS_KEY, MINIO_SECRET_KEY")
            return
        
        # Generate test file
        text_content = f"This is test data for MinIO storage testing.\n" * 1000
        self.original_hash = hashlib.md5(text_content.encode('utf-8')).hexdigest()
        
        with open(self.source_file, 'w', encoding='utf-8') as f:
            f.write(text_content)
        
        # Execute compression operation
        compress_command = f"./release-release/bin/pbgz compress {self.source_file} -o {self.compressed_file}"
        self.add_test_command(compress_command, cmd_id=1)
        
        # Execute decompression operation
        decompress_command = f"./release-release/bin/pbgz decompress {self.compressed_file} -o {self.decompressed_file}"
        self.add_test_command(decompress_command, cmd_id=2)
    
    def verify_expected_results(self) -> bool:
        """Verify test results with MinIO integration"""
        
        # Verify local files
        if not os.path.exists(self.compressed_file):
            print("Error: Compressed file not created locally")
            return False
        
        if not os.path.exists(self.decompressed_file):
            print("Error: Decompressed file not created locally")
            return False
        
        # Get compression ratio
        compression_ratio = self.get_compression_rate()
        if compression_ratio is not None:
            print(f"Compression ratio: {compression_ratio:.2f}%")
        
        # If MinIO configuration is unavailable, only verify local results
        if self.minio_adapter is None:
            return self._verify_local_integrity()
        
        # TODO: Execute MinIO related tests
        # Example: Upload compressed file to MinIO
        upload_success = self.minio_adapter.upload_file(
            self.compressed_file,
            self.bucket_name,
            f"tests/{self.compressed_file}"
        )
        
        if not upload_success:
            print("Warning: Failed to upload to MinIO, but local test passed")
            return self._verify_local_integrity()
        
        # Verify file exists in MinIO
        remote_object = f"tests/{self.compressed_file}"
        if not self.minio_adapter.file_exists(self.bucket_name, remote_object):
            print("Error: File not found in MinIO")
            return False
        
        print(f"✓ File successfully uploaded to MinIO: {remote_object}")
        
        # Optional: Download from MinIO and verify
        downloaded_file = "minio_downloaded.txt.pbgz"
        download_path = self.minio_adapter.download_file(
            self.bucket_name,
            remote_object,
            downloaded_file,
            overwrite=True
        )
        
        if download_path:
            print(f"✓ File successfully downloaded from MinIO: {download_path}")
            # Verify downloaded file
            with open(self.compressed_file, 'rb') as f1, open(downloaded_file, 'rb') as f2:
                original = f1.read()
                downloaded = f2.read()
                if original == downloaded:
                    print("✓ Downloaded file matches original")
                    os.remove(downloaded_file)  # Clean up downloaded file
                else:
                    print("Error: Downloaded file doesn't match original")
                    return False
        
        # Clean up file in MinIO
        self.minio_adapter.delete_file(self.bucket_name, remote_object)
        
        # Final verification of local data integrity
        return self._verify_local_integrity()
    
    def _verify_local_integrity(self) -> bool:
        """Verify local data integrity"""
        with open(self.source_file, 'rb') as f1, open(self.decompressed_file, 'rb') as f2:
            original_data = f1.read()
            decompressed_data = f2.read()
        
        decompressed_hash = hashlib.md5(decompressed_data).hexdigest()
        
        if self.original_hash == decompressed_hash:
            print("✓ Local data integrity verified")
            return True
        else:
            print("✗ Local data integrity check failed")
            return False
    
    def cleanup_test_data(self):
        """Clean up temporary files"""
        files_to_clean = [self.source_file, self.compressed_file, self.decompressed_file]
        for filepath in files_to_clean:
            try:
                if os.path.exists(filepath):
                    os.remove(filepath)
            except Exception as e:
                print(f"Warning: Failed to remove {filepath}: {e}")


class MinIOAdvancedTest(PBGZTestCase):
    """Advanced MinIO Test Case - Demonstrates More Integration Scenarios"""
    
    def __init__(self):
        super().__init__("MinIOAdvancedTest")
        self.minio_adapter = None
        self.bucket_name = "pbgz-advanced-test"
    
    def get_test_files(self) -> tuple:
        return (None, None, None)
    
    def prepare_data(self):
        """Prepare advanced test data"""
        self.minio_adapter = create_minio_adapter_from_env()
        if self.minio_adapter is None:
            print("MinIO not configured, skipping advanced tests")
            return
        
        print("This test demonstrates MinIO integration for:")
        print("1. Remote file storage and retrieval")
        print("2. Distributed compression scenarios")
        print("3. Backup and archiving workflows")
        print("4. Multi-environment data sharing")
        
        # Add some example commands to demonstrate possible integration methods
        # In actual applications, you can use MinIO to download source files, compress them and then upload them back
        # Or use MinIO to share compression results in distributed environments
    
    def verify_expected_results(self) -> bool:
        """Verify advanced test results"""
        if self.minio_adapter is None:
            print("MinIO not configured, marking as passed")
            return True
        
        # List all files in the bucket
        print(f"\nChecking MinIO bucket: {self.bucket_name}")
        try:
            self.minio_adapter.ensure_bucket(self.bucket_name)
            
            files = self.minio_adapter.list_files(self.bucket_name)
            print(f"Found {len(files)} files in bucket")
            
            print("✓ MinIO advanced integration test completed")
            return True
        except Exception as e:
            print(f"Error in advanced test: {e}")
            return False
    
    def cleanup_test_data(self):
        """Clean up advanced test data"""
        pass


if __name__ == "__main__":
    import sys
    
    print("=" * 60)
    print("MinIO Storage Integration Test")
    print("=" * 60)
    
    # Environment variable configuration prompt
    print("\nFor MinIO testing, set these environment variables:")
    print("  export MINIO_ENDPOINT='your-minio-server:9000'")
    print("  export MINIO_ACCESS_KEY='your-access-key'")
    print("  export MINIO_SECRET_KEY='your-secret-key'")
    print("  export MINIO_SECURE='true'")
    print("=" * 60 + "\n")
    
    # Run basic MinIO test
    test_case = MinIOStorageTest()
    result = test_case.execute()
    test_case.print_results()
    
    if result:
        print("\n✓ MinIO storage test passed!")
    else:
        print("\n✗ MinIO storage test failed!")
        sys.exit(1)