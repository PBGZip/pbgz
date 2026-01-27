# Security Policy

## Supported Versions

| Version | Supported          |
|---------|--------------------|
| 2.1     | ✅                 |
| 2.0     | ✅                 |
| < 2.0   | ❌                 |

## Security Contact Information

To report a security vulnerability, please use the following channels:

### Primary Contact
- **Email**: wenjinyang2729@phytiumn.com.cn
- **PGP Key**: 642C 2089 D438 0E27 E045 59D2 6068 D10B AF91 9AA7

### Alternative Contacts
- **GitHub Issues**: Use the private vulnerability reporting feature on GitHub

## Vulnerability Disclosure Process

### 1. Discovery and Reporting
- **Initial Contact**: Report discovered vulnerabilities to the primary email
- **PGP Encryption**: Please encrypt sensitive vulnerability details using the PGP key
- **Information Required**: Include detailed vulnerability descriptions and proof-of-concept
- **Response Time**: Initial response within 48-72 hours

### 2. Assessment and Classification
- **Triage Process**: Security team assesses impact and severity within 5 business days
- **Severity Levels**: 
  - **Critical**: Direct code execution, privilege escalation, data compromise
  - **High**: Significant security impact with limited exploitation potential
  - **Medium**: Moderate security impact with user interaction required
  - **Low**: Minimal security impact

### 3. Resolution Timeline
| Severity | Fix Time | Public Disclosure |
|----------|-----------|-------------------|
| Critical | 2-4 weeks | At fix release |
| High | 4-8 weeks | At fix release |
| Medium | 8-12 weeks | At fix release |
| Low | Next release | At fix release |

### 4. Disclosure Policy
- **Coordinated Disclosure**: We follow coordinated vulnerability disclosure
- **Public Disclosure**: Details disclosed after fix is available
- **Credit**: Security researchers credited in advisories and acknowledgments

## Security Scope

### In Scope
- PBGZ core compression/decompression algorithms
- File handling and parsing components
- Network communication protocols
- Authentication and authorization mechanisms
- Data encryption and integrity verification
- Input validation and sanitization
- Memory management and buffer handling
- Configuration and parameter processing

### Out of Scope
- Dependencies and third-party libraries (handled separately)
- Physical security of hardware components
- Social engineering attacks
- DoS attacks against network infrastructure
- Vulnerabilities in outdated versions

## Security Assurance

### Threat Modeling
- **Data Centered Approach**: Focus on genomic data protection throughout lifecycle
- **Attack Surfaces**: Identify and mitigate potential attack vectors
- **Data Flows**: Secure data handling across compression pipeline
- **Trust Boundaries**: Clear separation between components and privilege levels

### Secure Development
- **Input Validation**: Comprehensive validation of FASTQ/FASTA inputs
- **Memory Safety**: Secure memory management to prevent data leaks
- **Configuration Security**: Secure handling of compression parameters and keys
- **Error Handling**: Secure error reporting without information disclosure

### Testing and Verification
- **Fuzz Testing**: Comprehensive fuzz testing of input parsing
- **Penetration Testing**: Regular security assessments
- **Code Audits**: Static and dynamic security analysis
- **Regression Testing**: Security regression test suite

## Data Security Classification

### Genomic Data Sensitivity
- **PHI (Protected Health Information)**: Personal identifiers
- **Genetic Information**: Sensitive genomic variants, disease markers
- **Research Data**: Proprietary or confidential research results
- **Clinical Data**: Patient diagnostic and treatment information

### Protection Requirements
- **Encryption**: AES encryption for data at rest and in transit
- **Access Control**: Role-based access to compressed data
- **Audit Trail**: Complete audit logging of data operations
- **Integrity**: MD5 integrity verification of compressed files
- **Retention**: Secure data handling and deletion procedures

## Security Best Practices

### For Users
- **Access Control**: Restrict access to compressed files using file permissions
- **Network Security**: Use secure channels for file transfer
- **Backup Security**: Secure backup of compression configuration files
- **Monitoring**: Monitor compression/decompression operations
- **Updates**: Keep software updated to latest secure versions

### For Administrators
- **System Hardening**: Apply OS and filesystem security configurations
- **Resource Limits**: Set appropriate CPU, memory, and disk limits
- **Logging**: Enable comprehensive logging and monitoring
- **Vulnerability Management**: Regular security scanning and updates
- **Incident Response**: Establish security incident procedures

## Security Related Documents

For detailed technical security information, please refer to:
- **[Security Analysis Report](security_analysis_report.md)**: Comprehensive security analysis
- **[Threat Assessment](docs/threat_model.md)**: Threat model and analysis
- **[Security Best Practices](docs/security_guide.md)**: Security configuration guide
- **[Vulnerability Management](docs/vulnerability_process.md)**: Vulnerability handling process

## Incidence Response

### Security Incident Categories
- **Data Compromise**: Unauthorized access to sensitive genomic data
- **Service Disruption**: Attacks affecting compression service availability
- **_privilege Escalation**: Unauthorized elevation of system privileges
- **Data Integrity**: Unauthorized modification of compressed data
- **System Compromise**: Full system takeover or malware infection

### Response Process
1. **Detection**: Security monitoring and alerting
2. **Assessment**: Incident classification and impact analysis
3. **Containment**: Immediate isolation of affected systems
4. **Eradication**: Remove threat and security vulnerabilities
5. **Recovery**: Restore secure operations
6. **Lessons Learned**: Post-incident review and improvement

### Reporting Contacts
| Time Zone | Contact Information | Response Time |
|-----------|-------------------|---------------|
| Weekdays (9:00-18:00 UTC+8) | Security Team | < 2 hours |
| Weekends/Holidays | Security Team | < 8 hours |
| Critical Incidents | Security Team | < 1 hour |

## Security Tools

### Recommended Security Tools
- **Static Analysis**: CodeQL, SonarQube
- **Dynamic Analysis**: Valgrind, AddressSanitizer
- **Fuzz Testing**: AFL++, libFuzzer
- **Dependency Scanning**: OWASP Dependency-Check
- **Container Security**: Docker Security Scanning
- **File Integrity**: AIDE, Tripwire

### Security Configuration
- **SELinux/AppArmor**: Mandatory access control
- **Firewall**: Network traffic filtering
- **Intrusion Detection**: OSSEC, Snort
- **Log Management**: ELK Stack, Splunk
- **Security Monitoring**: Wazuh, OpenVAS


### Dependency Security Updates

We regularly audit and update dependencies:

- 📅 **Regular Checks**: Monthly security advisory checks
- 🔄 **Automatic Updates**: Critical security patches updated automatically
- 🧪 **Test Verification**: Comprehensive testing after updates

## Threat Model

### 🎯 Potential Threats

#### Data Integrity Threats

- 🗂️ **File Tampering**: Malicious modification of compressed files
- 🔍 **Index Corruption**: Index file corruption or tampering
- 📊 **Data Leakage**: Sensitive data leakage through compressed files

#### Availability Threats

- 💥 **Denial of Service**: Program crashes through malicious input
- 🔄 **Resource Exhaustion**: Resource exhaustion from large file processing
- ⏰ **Performance Degradation**: Performance degradation through specially crafted input

#### Configuration Threats

- ⚙️ **Configuration Injection**: Malicious configuration files
- 🏗️ **Build Attacks**: Security vulnerabilities in build process
- 🔌 **Plugin Security**: Security risks from third-party plugins

### 🛡️ Mitigation Measures

#### Input Validation

```cpp
// Strict validation of all inputs
class InputValidator {
public:
    static bool validate_file_size(size_t size) {
        return size <= MAX_FILE_SIZE;
    }
    
    static bool validate_compression_level(int level) {
        return level >= MIN_COMPRESSION_LEVEL && 
               level <= MAX_COMPRESSION_LEVEL;
    }
    
    static bool validate_thread_count(int count) {
        return count > 0 && count <= MAX_THREADS;
    }
};
```

#### Error Handling

```cpp
// Secure error handling
try {
    // Dangerous operation
    process_file(file_path);
} catch (const std::exception& e) {
    // Log errors without exposing sensitive information
    log_error("Processing failed: " + std::string(e.what()));
    // Secure recovery
    cleanup();
    return false;
}
```

## Security Testing

### 🧪 Testing Strategy

#### Unit Testing

- Boundary condition testing
- Exceptional input testing
- Error path testing

#### Integration Testing

- File processing testing
- Network security testing
- Concurrency security testing

#### Penetration Testing

- Fuzz testing
- Boundary scanning
- Code audit

### 📊 Security Metrics

- 🐛 Vulnerability density
- ⏱️ Mean time to fix
- 📈 Security coverage rate
- ✅ Compliance checks

## Compliance

### 📋 Regulatory Compliance

PBGZ follows relevant security regulations and standards:

- **Data Protection**: Follow data protection best practices
- **Export Control**: Comply with relevant export control regulations
- **License**: See [LICENSE](LICENSE) file for license information

### 🔍 Audit and Certification

- 📊 **Code Audits**: Regular code security audits
- 🏆 **Security Certification**: Plan to obtain relevant security certifications
- 📝 **Compliance Reports**: Regular publication of compliance reports

## Community Security

### 👥 Security Team

Our security team is responsible for:

- 🔍 Vulnerability assessment and remediation
- 📋 Security policy development
- 🛠️ Security tool development
- 📚 Security awareness training

### 🤝 Community Participation

We encourage community participation in security work:

- 🐛 Report potential security issues
- 💡 Propose security improvements
- 🔬 Participate in security testing
- 📖 Contribute security documentation

## Contact Information

### 🆘 Emergency Contact

If you need to report urgent security issues:

- 📧 **Security Email**: wenjinyang2729@phytiumn.com.cn
- 🔑 **PGP Key**: 642C 2089 D438 0E27 E045 59D2 6068 D10B AF91 9AA7

### 📞 General Inquiries

For general security-related questions:

- 💬 **GitHub Discussions**: Discuss in security channel
- 📧 **Project Email**: wenjinyang2729@phytiumn.com.cn

## Acknowledgments

We thank the following individuals and organizations for their contributions to PBGZ security:

- 🌟 **Security Researchers**: Researchers who reported vulnerabilities
- 🏢 **Organizations**: Supporting institutions
- 👥 **Community Members**: Community members who participated in security testing

## Changelog

### Security Updates

| Date | Version | Update Type | Description |
|------|---------|-------------|-------------|
| 2026-1-26 | v2.1.0 | 🛡️ Security Enhancement | Initial security framework |
| | | | |

---

**Security is everyone's responsibility**. If you discover any security issues or have suggestions for improvement, please let us know. Together, let's maintain the security of the PBGZ project! 🔒

---

*Last updated: January 2026*
