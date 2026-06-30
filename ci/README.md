# PBGZ CI System - Gradle Integration Version

Feature-rich continuous integration system supporting multiple trigger modes, modular script design, complete build verification process, deep Gradle integration, and configuration file management for non-sensitive parameters.

## Features

- 🔧 **Modular Design**: 7 independent, debuggable, reusable functional scripts
- 🚀 **Multiple Trigger Modes**: Manual build, scheduled build, PR build  
- ☁️ **Dynamic Cloud Server**: Automatic creation and release of build servers
- 🎯 **Gradle Integration**: Use Gradle to drive complete CI process
- 📦 **Unified Management**: Download, build, test uniformly scheduled through Gradle
- ⚙️ **Configuration File Management**: Non-sensitive parameters managed through configuration files
- 🔄 **Smart Download**: Incremental updates, 3 retries, integrity checking
- 🧪 **Complete Testing**: Automated test suite integration
- 📊 **Artifact Collection**: Logs and binary artifacts automatically collected
- 🛡️ **Security Cleanup**: SSH keys and temporary resources automatically cleaned up

## Quick Start

### 1. Configure GitHub Secrets (Required)

Configure the following Secrets in GitHub repository settings (Sensitive configuration):

| Secret Name | Description | Example |
|------------|------------|--------| 
| `BUILD_SERVER_ACCESS_KEY_ID` | Cloud server access key ID | `LTAI5t1234567890abcdef` |
| `BUILD_SERVER_ACCESS_KEY_SECRET` | Cloud server access key secret | `abc123def456...` |
| `BUILD_SERVER_SSH_KEY` | SSH private key complete content | `-----BEGIN RSA PRIVATE KEY-----...` |

### 2. Configure Parameter File (Required)

Edit the `ci-config.properties` file to configure non-sensitive parameters:

```properties
# Cloud server configuration
BUILD_SERVER_REGION=cn-hangzhou
INSTANCE_TYPE=ecs.c6.xlarge

# Target project configuration
GITHUB_REPO=https://github.com/PBGZip/pbgz.git
GITHUB_BRANCH=pbgz_v2.2.0

# Build script configuration
BUILD_SCRIPT_GCC=3rd_party/build.all.gcc.sh
BUILD_SCRIPT_RELEASE=build-release.sh

# CI script repository configuration
GITHUB_CI_REPO=https://github.com/your-username/pbgz_ci.git
GITHUB_CI_BRANCH=main

# Test configuration
TEST_ENABLED=false
```

**Note**: The `ci-config.properties` file contains detailed configuration instructions and comments. Modify configuration values as needed.

### 3. Choose Trigger Mode

#### Manual Build
- **Use Case**: Manual verification, debugging, testing
- **Trigger Method**: GitHub Actions UI → manual-build → Run workflow

#### Scheduled Build
- **Use Case**: Daily automated build verification
- **Trigger Method**: Automatically executes daily at UTC 02:00 (Beijing time 10:00)

#### PR Build
- **Use Case**: PR code verification, quality assurance
- **trigger method**：CreateupdatePRAutomaticexecute

## Configuration Details

### Division Between GitHub Secrets and Configuration File

**GitHub Secrets** (Sensitive Information):
- `BUILD_SERVER_ACCESS_KEY_ID` - Access credential ID
- `BUILD_SERVER_ACCESS_KEY_SECRET` - Access credential secret  
- `BUILD_SERVER_SSH_KEY` - SSH private key

**ci-config.properties** (Non-Sensitive Configuration):
- `BUILD_SERVER_REGION` - Server region
- `INSTANCE_TYPE` - Instance type
- `GITHUB_REPO` - Target project repository URL
- `GITHUB_BRANCH` - Target branch
- `BUILD_SCRIPT_GCC` - GCC build script path
- `BUILD_SCRIPT_RELEASE` - Release build script path
- `GITHUB_CI_REPO` - CI script repository URL
- `GITHUB_CI_BRANCH` - CI script branch
- `TEST_ENABLED` - Whether to enable testing

### Configuration File Priority

1. **ci-config.properties**: Configuration file in project root directory
2. **Default Values**: Use system defaults if configuration file does not exist
3. **Environment Variables**: Can be overridden in workflows through environment variables

## Project Structure

```
pbgz_ci/
├── .github/workflows/           # GitHub Actions configuration
│   ├── manual-build.yml        # Manual trigger build
│   ├── scheduled-build.yml     # Scheduled trigger build
│   └── pr-build.yml            # PR trigger build
├── scripts/                     # CI script library (functional naming)
│   ├── create-build-server.sh   # Create build server
│   ├── prepare-build-environment.sh  # Prepare build environment
│   ├── run-gradle-build.sh     # Run Gradle build (integrated version)
│   ├── verify-test-results.sh   # verifyTest Results
│   ├── collect-artifacts.sh     # Collect Build Artifacts
│   ├── terminate-build-server.sh # terminateBuildServer
│   └── cleanup-resources.sh     # cleanupresources
├── ci-config.properties         # CIConfigurationfile (alreadyContains Detailed Explanation)
├── build.gradle.kts            # GradleBuild Configuration
├── gradle.properties           # GradlePropertiesConfiguration
├── gradlew                     # GradleScript
├── test-gradle-integration.sh   # Integration Test Script
├── README.md                   # Projectdocumentation
└── INTEGRATION_COMPLETE_REPORT.md # Integrationcompletereport
```

## Usage

### Local TestingGradle CI

```bash
# Run complete CI Process
./gradlew ci

# View Available Tasks
./gradlew tasks --group=ci
```

### executeScript

```bash
bash scripts/create-build-server.sh
bash scripts/prepare-build-environment.sh
bash scripts/run-gradle-build.sh
```

### RunIntegration Test

```bash
./test-gradle-integration.sh
```

## Configuration

### BuildParameters

 `ci-config.properties` fileBuildParameters：

```properties
# Change Cloud Server Configuration
BUILD_SERVER_REGION=cn-beijing
INSTANCE_TYPE=ecs.c6.2xlarge

# Change Target Project
GITHUB_REPO=https://github.com/your-username/your-project.git
GITHUB_BRANCH=develop

# Build Script
BUILD_SCRIPT_GCC=scripts/custom-gcc-build.sh
BUILD_SCRIPT_RELEASE=scripts/custom-release-build.sh

# Test
TEST_ENABLED=true
```

**Configuration**：`ci-config.properties` fileinAllConfiguration of theDetailed Explanation、exampleAnd，file of the。

### environment variables

ifneedtemporaryConfigurationfileIn TheParameters，canatWorkflowinSettingenvironment variables：

```yaml
env:
  BUILD_SERVER_REGION: cn-shanghai
  INSTANCE_TYPE: ecs.c8.xlarge
  TEST_ENABLED: "true"
```

## 

### Common Problems

1. **Configurationfilenot**
   - Ensure `ci-config.properties` fileatProjectdirectory
   - checkfileformatwhether（key=valueformat）
   - View Workflow LogsIn TheConfigurationInformation

2. **GradleBuild Failed**
   - verifyGradleversionAndservice
   - checkGradleConfigurationwhether
   - View Detailed of theGradleBuild Log

3. **CIScriptDownloadFailure**
   - verifyGITHUB_CI_REPOConfiguration
   - checkRepositoryAccessPermissions
   - confirmBranch Namewhether

4. **SecretsConfiguration Error**
   - confirm of theSecretsalreadyConfiguration
   - checkkeyformatwhether
   - verifykeyPermissionswhetherenough

## SecurityBest Practices

✅ **Recommended Practices**:
- UseConfigurationfileSensitiveParameters
- SensitiveInformationstorageatGitHub Secretsin
- checkAndupdateConfiguration
- notenvironmentUsenot of theConfigurationfile

❌ **Avoidance Practices**:
- atConfigurationfileinstorageSensitiveInformation
- atinencodingkey
- ConfigurationfileCommittoversionSensitiveInformation

## guide

1. ForkProject
2. CreateBranch
3. Commit
4. Local Run`./gradlew ci`Test
5. toBranch
6. CreatePull Request

## version

### v3.1.1 - Configurationfileversion
- ✅ SensitiveParametersConfigurationfile
- ✅ Createci-config.propertiesConfigurationfile（Contains Detailed Explanation）
- ✅ updateWorkflowSupportConfigurationfile
- ✅ Configurationsecurity
- ✅ Simplify Configuration Management Process

### v3.1.0 - GradleIntegration Version
- ✅ CIProcessWithGradleDeep Integration
- ✅ UseGradleBuild、Test Process
- ✅ Supports Local TestingCIProcess
- ✅ optimizeWorkflowStep，More Concise

### v3.0.0 - modularversion
- ✅ 7FunctionalScriptmodule
- ✅ no，
- ✅ trigger method（Manual、Scheduled、PR）

## 

Project of theProject。

## Contact Information

If There Are Problems，Please Contact Via The Following Methods：
- CommitIssue
- Create
- Pull Request

---

**Currentversion**: v3.1.1（Configurationfileversion）

**Configuration**: 
- SensitiveInformationUseGitHub Secrets (3)
- SensitiveParametersUseci-config.properties (alreadyContains Detailed Explanation)
- ConfigurationfilealreadyatWorkflowinIntegration
- ci-config.propertiesConfiguration