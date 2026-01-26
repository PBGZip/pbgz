# Contributing Guidelines

Thank you for your interest in the PBGZ project! We welcome various forms of contributions, including but not limited to:

- 🐛 Bug reports
- 💡 New feature suggestions
- 📝 Documentation improvements
- 🔧 Code contributions
- 🧪 Test cases
- 🌐 Translation work

## Development Environment Setup

### System Requirements

- **Operating System**: Linux (Ubuntu, Kylin), macOS
- **Architecture**: x86_64, ARMv8
- **Compiler**: GCC 7+ or Clang 8+
- **CMake**: 3.10+
- **Git**: 2.0+

### Dependencies

PBGZ depends on the following third-party libraries:

- **zlib** (1.2.11) - Compression library
- **bzip2** - Compression library
- **zstd** (1.5.0) - Zstandard compression library
- **libdeflate** (1.7) - High-performance compression library
- **htslib** (1.12) - Genomic data processing library
- **jsoncpp** (1.9.4) - JSON processing library

All dependency source code is included in the `3rd_party/` directory.

### Build Steps

1. **Clone Repository**
   ```bash
   git clone https://github.com/PBGZip/pbgz.git
   cd pbgz
   ```

2. **Build Dependencies**
   ```bash
   cd 3rd_party
   ./build.all.gcc.sh
   cd ..
   ```

3. **Build Project**
   
   **Debug Version**:
   ```bash
   ./build-debug.sh
   ```
   
   **Release Version**:
   ```bash
   ./build-release.sh
   ```
   
   **Performance Optimized Version**:
   ```bash
   ./build-release-profile.sh
   ```

## Contribution Process

### 1. Fork and Clone

1. Fork the project to your GitHub account
2. Clone your fork:
   ```bash
   git clone https://github.com/YOUR_USERNAME/pbgz.git
   cd pbgz
   git remote add upstream https://github.com/PBGZip/pbgz.git
   ```

### 2. Create Branch

```bash
git checkout -b feature/your-feature-name
# or
git checkout -b fix/your-bug-fix
```

### 3. Development and Testing

- Please follow the project's [coding standards](#coding-standards) when writing code
- Add necessary test cases
- Ensure all existing tests pass

```bash
# Run tests
cd build
make test
# or
ctest
```

### 4. Commit Code

```bash
git add .
git commit -m "feat: add new feature description"
# or
git commit -m "fix: fix issue description"
```

#### Commit Message Format

Use [Conventional Commits](https://www.conventionalcommits.org/) format:

- `feat:` New feature
- `fix:` Bug fix
- `docs:` Documentation update
- `style:` Code formatting adjustments
- `refactor:` Code refactoring
- `test:` Test related
- `perf:` Performance optimization
- `ci:` CI/CD related

### 5. Push and Pull Request

```bash
git push origin feature/your-feature-name
```

Then create a Pull Request on GitHub.

## Coding Standards

### C++ Code Style

- **Indentation**: Use 4 spaces, no tabs
- **Naming**:
  - Class names: `PascalCase` (e.g., `PbgzEngine`)
  - Function names: `camelCase` (e.g., `compressFile`)
  - Variable names: `camelCase` (e.g., `filePath`)
  - Constants: `UPPER_SNAKE_CASE` (e.g., `MAX_BUFFER_SIZE`)
  - Macros: `UPPER_SNAKE_CASE` (e.g., `MAX_BUFFER_SIZE`)

- **Comments**:
  - Classes and functions in header files must have detailed comments
  - Complex algorithms need inline comments
  - Comments can be in Chinese or English, maintain consistency

- **Braces**: K&R style
  ```cpp
  if (condition) {
      // code
  } else {
      // code
  }
  ```

### File Organization

- **Header files**: Only contain declarations, avoid `using namespace`
- **Source files**: Contain implementations, grouped logically
- **Unit tests**: Place in `test/testcase/` directory

### Performance Considerations

PBGZ is a performance-sensitive project, please pay special attention to:

- Avoid unnecessary memory allocations
- Use move semantics to reduce copies
- Avoid virtual function calls on critical paths
- Consider cache-friendly data structures

## Testing Guidelines

### Test Types

1. **Unit Tests**: Test individual functions/classes
2. **Integration Tests**: Test component interactions
3. **Performance Tests**: Verify compression ratio and speed
4. **Regression Tests**: Ensure fixes don't introduce new issues

### Running Tests

```bash
# Build tests
cd build
make test

# Run specific tests
./test/pbgz_test --gtest_filter=YourTestCase
```

### Adding Tests

- New features must include test cases
- Bug fixes need regression tests
- Performance-critical code needs benchmark tests

## Bug Reports

Please include when reporting bugs:

1. **Environment Information**:
   - Operating system and version
   - CPU architecture
   - Compiler version
   - PBGZ version

2. **Reproduction Steps**:
   - Detailed operation steps
   - Input file information
   - Expected behavior vs actual behavior

3. **Error Information**:
   - Complete error logs
   - Stack traces (if any)

4. **Additional Information**:
   - Minimal reproduction example
   - Related configuration files

## Feature Requests

When requesting new features, please describe:

1. **Use Case**: Why do you need this feature
2. **Expected Behavior**: How should the feature work
3. **Alternatives**: Have you considered other solutions
4. **Additional Information**: Related links, reference materials, etc.

## Documentation Contributions

Documentation is an important part of the project, we welcome:

- API documentation improvements
- Usage examples additions
- Installation guide updates
- Performance tuning suggestions
- Translation work

## Release Process

The project's release process:

1. **Version Planning**: Discuss in GitHub Issues
2. **Development Phase**: Feature development and testing
3. **Code Review**: All PRs need at least one maintainer review
4. **Integration Testing**: Comprehensive testing on all platforms
5. **Release Preparation**: Update version numbers and CHANGELOG
6. **Official Release**: Create Git tags and GitHub Release

## Community Guidelines

- 🤝 **Respect**: Respect all participants, maintain friendliness and professionalism
- 🧠 **Learn**: Continuous learning, willing to share knowledge
- 🎯 **Focus**: Maintain focus on project goals
- 📈 **Improve**: Continuously improve code and processes
- 🤲 **Collaborate**: Encourage collaboration, avoid duplicate work

## Getting Help

If you need help or have questions:

1. Check [documentation](../docs/)
2. Search [Issues](https://github.com/PBGZip/pbgz/issues)
3. Create a new Issue
4. Discuss in GitHub Discussions

## Acknowledgments

Thank you to all developers who have contributed to the PBGZ project! Your contributions make this project better.

---

Thank you again for your contribution! 🎉
