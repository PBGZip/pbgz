#!/usr/bin/env python3
"""
MinIO Configuration Verification Script
Verify all required configurations are set correctly
"""
import os
import sys
import json

def load_ci_config():
    """Load CI configuration file"""
    config = {}
    config_file = os.path.join(os.path.dirname(__file__), '..', '..', 'ci', 'ci-config.properties')
    
    try:
        with open(config_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                if '=' in line:
                    key, value = line.split('=', 1)
                    key = key.strip()
                    value = value.strip()
                    config[key] = value
    
        return config
    except Exception as e:
        print(f"Warning: Could not load CI config: {e}")
        return {}

def minio_config_check():
    """Check MinIO configuration"""
    print("=" * 60)
    print("🗄️ MinIO Configuration Check")
    print("=" * 60)
    
    # Load non-sensitive configuration
    ci_config = load_ci_config()
    
    # Check sensitive configuration (environment variables)
    sensitive_configs = [
        ('MINIO_ACCESS_KEY', 'MinIO Access Key'),
        ('MINIO_SECRET_KEY', 'MinIO Secret Key')
    ]
    
    # Check non-sensitive configuration
    non_sensitive_configs = [
        ('MINIO_ENDPOINT', 'MinIO Service Endpoint'),
        ('MINIO_SECURE', 'Use HTTPS'),
        ('MINIO_REGION', 'Storage Region'),
        ('MINIO_VERIFY_SSL', 'Verify SSL Certificate'),
        ('MINIO_COMPRESSED_BUCKET', 'Compressed Files Bucket'),
        ('MINIO_SOURCE_BUCKET', 'Source Files Bucket'),
        ('MINIO_RESULTS_BUCKET', 'Results Files Bucket'),
        ('MINIO_BACKUP_BUCKET', 'Backup Files Bucket')
    ]
    
    issues = []
    warnings = []
    
    print("\n🔐 Sensitive Configuration Check (should be read from GitHub Secrets):")
    print("-" * 60)
    
    for env_key, description in sensitive_configs:
        value = os.getenv(env_key)
        if value:
            masked = value[:4] + '*' * (len(value) - 4) if len(value) > 4 else '****'
            print(f"✅ {description:30} : {masked}")
        else:
            print(f"❌ {description:30} : Not Set")
            issues.append(f"Missing Environment Variable: {env_key}")
    
    print("\n📝 Non-sensitive Configuration Check (read from ci-config.properties):")
    print("-" * 60)
    
    for config_key, description in non_sensitive_configs:
        value = os.getenv(config_key)
        if value:
            print(f"✅ {description:30} : {value}")
        else:
            print(f"⚠️  {description:30} : Not Set (using default value)")
            warnings.append(f"Optional Configuration: {config_key}")
    
    # Test importing MinIO module
    print("\n🔧 Module Import Test:")
    print("-" * 60)
    
    try:
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
        from minio_storage_adapter import MinIOStorageAdapter, create_minio_adapter_from_env
        print("✅ MinIO module import successful")
        
        # Test creating adapter
        adapter = create_minio_adapter_from_env()
        if adapter:
            print("✅ MinIO adapter created successfully")
        else:
            print("⚠️  MinIO adapter not created (incomplete configuration)")
            warnings.append("MinIO adapter not created")
            
    except ImportError as e:
        print(f"❌ MinIO module import failed: {e}")
        issues.append("Installation Missing: pip install minio")
    except Exception as e:
        print(f"❌ MinIO module error: {e}")
        issues.append(f"Module Error: {str(e)}")
    
    # Output result summary
    print("\n" + "=" * 60)
    print("📊 Check Results Summary")
    print("=" * 60)
    
    if not issues and not warnings:
        print("✅ All configuration checks passed! MinIO integration can be used normally.")
        return 0
    elif issues and not warnings:
        print(f"❌ Found {len(issues)} issue(s) that need to be fixed:")
        for issue in issues:
            print(f"   - {issue}")
        return 1
    elif not issues and warnings:
        print(f"⚠️  Found {len(warnings)} warning(s), but basic functionality is not affected:")
        for warning in warnings:
            print(f"   - {warning}")
        return 0
    else:
        print(f"❌ Found {len(issues)} issue(s) and {len(warnings)} warning(s):")
        print("Issues:")
        for issue in issues:
            print(f"   - {issue}")
        print("Warnings:")
        for warning in warnings:
            print(f"   - {warning}")
        return 1

def main():
    """Main function"""
    print("🔍 PBGZ MinIO Integration Configuration Verification")
    print()
    
    result = minio_config_check()
    
    if result == 0:
        print("\n🎉 Configuration verification passed! You can start using MinIO integration.")
        print("📚 Detailed Usage Instructions: ci/testcase/MINIO_INTEGRATION.md")
        print("🔐 GitHub Secrets Configuration: ci/GITHUB_SECRETS.md")
    else:
        print("\n❌ Configuration verification failed, please check the issues above.")
        print("💡 Configuration Help: ci/GITHUB_SECRETS.md")
    
    return result

if __name__ == "__main__":
    sys.exit(main())