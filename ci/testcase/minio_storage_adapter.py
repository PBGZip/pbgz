"""
MinIO Storage Adapter for PBGZ
Provides file upload and download functionality for MinIO object storage
"""
import os
import ssl
from typing import Optional
from urllib.parse import urlparse
from minio import Minio
from minio.error import S3Error


class MinIOStorageAdapter:
    """MinIO Storage Adapter"""
    
    def __init__(
        self,
        endpoint: str,
        access_key: str,
        secret_key: str,
        secure: bool = True,
        region: Optional[str] = None,
        verify_ssl: bool = True
    ):
        """
        Initialize MinIO client
        
        Args:
            endpoint: MinIO service address, e.g., 'play.min.io:9000'
            access_key: Access key
            secret_key: Secret key
            secure: Whether to use HTTPS
            region: Region setting
            verify_ssl: Whether to verify SSL certificate
        """
        self.endpoint = endpoint
        self.access_key = access_key
        self.secret_key = secret_key
        
        # Handle SSL verification
        http_client = None
        if not verify_ssl:
            http_client = self._create_insecure_http_client()
        
        self.client = Minio(
            endpoint,
            access_key=access_key,
            secret_key=secret_key,
            secure=secure,
            region=region,
            http_client=http_client
        )
    
    def _create_insecure_http_client(self):
        """Create HTTP client that does not verify SSL"""
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
        import urllib3
        return urllib3.PoolManager(
            cert_reqs='CERT_NONE',
            assert_hostname=False,
            ssl_context=ctx
        )
    
    def ensure_bucket(self, bucket_name: str):
        """
        Ensure bucket exists
        
        Args:
            bucket_name: Bucket name
        """
        try:
            if not self.client.bucket_exists(bucket_name):
                self.client.make_bucket(bucket_name)
                print(f"Created bucket: {bucket_name}")
            else:
                print(f"Bucket exists: {bucket_name}")
        except S3Error as e:
            print(f"Error ensuring bucket: {e}")
            raise
    
    def upload_file(
        self,
        local_file_path: str,
        bucket_name: str,
        object_name: Optional[str] = None,
        metadata: Optional[dict] = None
    ) -> bool:
        """
        Upload file to MinIO
        
        Args:
            local_file_path: Local file path
            bucket_name: Bucket name
            object_name: Object name (uses filename if not provided)
            metadata: File metadata
        
        Returns:
            Whether the upload was successful
        """
        if not os.path.exists(local_file_path):
            print(f"Local file not found: {local_file_path}")
            return False
        
        if object_name is None:
            object_name = os.path.basename(local_file_path)
        
        try:
            self.ensure_bucket(bucket_name)
            self.client.fput_object(
                bucket_name,
                object_name,
                local_file_path,
                metadata=metadata
            )
            print(f"Successfully uploaded: {local_file_path} -> {bucket_name}/{object_name}")
            return True
        except S3Error as e:
            print(f"Error uploading file: {e}")
            return False
    
    def download_file(
        self,
        bucket_name: str,
        object_name: str,
        local_file_path: Optional[str] = None,
        overwrite: bool = False
    ) -> Optional[str]:
        """
        Download file from MinIO
        
        Args:
            bucket_name: Bucket name
            object_name: Object name
            local_file_path: Local file path (uses object name if not provided)
            overwrite: Whether to overwrite existing file
        
        Returns:
            Downloaded file path, returns None on failure
        """
        if local_file_path is None:
            local_file_path = os.path.basename(object_name)
        
        if os.path.exists(local_file_path) and not overwrite:
            print(f"File already exists: {local_file_path}")
            return local_file_path
        
        try:
            directory = os.path.dirname(local_file_path)
            if directory and not os.path.exists(directory):
                os.makedirs(directory)
            
            self.client.fget_object(bucket_name, object_name, local_file_path)
            print(f"Successfully downloaded: {bucket_name}/{object_name} -> {local_file_path}")
            return local_file_path
        except S3Error as e:
            print(f"Error downloading file: {e}")
            return None
    
    def delete_file(self, bucket_name: str, object_name: str) -> bool:
        """
        Delete file from MinIO
        
        Args:
            bucket_name: Bucket name
            object_name: Object name
        
        Returns:
            Whether the deletion was successful
        """
        try:
            self.client.remove_object(bucket_name, object_name)
            print(f"Successfully deleted: {bucket_name}/{object_name}")
            return True
        except S3Error as e:
            print(f"Error deleting file: {e}")
            return False
    
    def file_exists(self, bucket_name: str, object_name: str) -> bool:
        """
        Check if file exists
        
        Args:
            bucket_name: Bucket name
            object_name: Object name
        
        Returns:
            Whether the file exists
        """
        try:
            self.client.stat_object(bucket_name, object_name)
            return True
        except S3Error:
            return False
    
    def list_files(self, bucket_name: str, prefix: str = "") -> list:
        """
        List files in bucket
        
        Args:
            bucket_name: Bucket name
            prefix: File prefix (for filtering)
        
        Returns:
            List of files
        """
        try:
            objects = self.client.list_objects(bucket_name, prefix=prefix)
            return [obj.object_name for obj in objects]
        except S3Error as e:
            print(f"Error listing files: {e}")
            return []


def create_minio_adapter_from_env() -> Optional[MinIOStorageAdapter]:
    """
    Create MinIO adapter from environment variables
    
    Environment Variables:
        MINIO_ENDPOINT: MinIO service endpoint
        MINIO_ACCESS_KEY: Access key
        MINIO_SECRET_KEY: Secret key
        MINIO_SECURE: Whether to use HTTPS (default True)
        MINIO_REGION: Region setting
        MINIO_VERIFY_SSL: Whether to verify SSL (default True)
    
    Returns:
        MinIOStorageAdapter instance, returns None when configuration is missing
    """
    endpoint = os.getenv('MINIO_ENDPOINT')
    access_key = os.getenv('MINIO_ACCESS_KEY')
    secret_key = os.getenv('MINIO_SECRET_KEY')
    
    if not all([endpoint, access_key, secret_key]):
        print("MinIO configuration not found in environment variables")
        return None
    
    return MinIOStorageAdapter(
        endpoint=endpoint,
        access_key=access_key,
        secret_key=secret_key,
        secure=os.getenv('MINIO_SECURE', 'true').lower() == 'true',
        region=os.getenv('MINIO_REGION'),
        verify_ssl=os.getenv('MINIO_VERIFY_SSL', 'true').lower() == 'true'
    )


def create_minio_adapter_from_url(minio_url: str) -> Optional[MinIOStorageAdapter]:
    """
    Create adapter from MinIO URL
    
    Args:
        minio_url: MinIO URL, format:
                   minio://access_key:secret_key@endpoint:port/bucket
                   or minio://access_key:secret_key@endpoint/bucket
    
    Returns:
        MinIOStorageAdapter instance
    """
    try:
        parsed = urlparse(minio_url)
        if parsed.scheme != 'minio':
            raise ValueError("URL scheme must be 'minio'")
        
        # Extract authentication information
        if '@' not in parsed.netloc:
            raise ValueError("Invalid MinIO URL format")
        
        auth, endpoint = parsed.netloc.split('@')
        access_key, secret_key = auth.split(':')
        
        return MinIOStorageAdapter(
            endpoint=endpoint,
            access_key=access_key,
            secret_key=secret_key,
            secure=parsed.scheme == 'https'
        )
    except Exception as e:
        print(f"Error creating MinIO adapter from URL: {e}")
        return None