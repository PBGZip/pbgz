# 贡献指南

感谢您对 PBGZ 项目的关注！我们欢迎各种形式的贡献，包括但不限于：

- 🐛 Bug 报告
- 💡 新功能建议
- 📝 文档改进
- 🔧 代码贡献
- 🧪 测试用例
- 🌐 翻译工作

## 开发环境设置

### 系统要求

- **操作系统**: Linux (Ubuntu, 麒麟), macOS
- **架构**: x86_64, ARMv8
- **编译器**: GCC 7+ 或 Clang 8+
- **CMake**: 3.10+
- **Git**: 2.0+

### 依赖项

PBGZ 依赖以下第三方库：

- **zlib** (1.2.11) - 压缩库
- **bzip2** - 压缩库
- **zstd** (1.5.0) - Zstandard 压缩库
- **libdeflate** (1.7) - 高性能压缩库
- **htslib** (1.12) - 基因组数据处理库
- **jsoncpp** (1.9.4) - JSON 处理库

所有依赖项的源码已包含在 `3rd_party/` 目录中。

### 构建步骤

1. **克隆仓库**
   ```bash
   git clone https://github.com/PBGZip/pbgz.git
   cd pbgz
   ```

2. **构建依赖项**
   ```bash
   cd 3rd_party
   ./build.all.gcc.sh
   cd ..
   ```

3. **构建项目**
   
   **调试版本**:
   ```bash
   ./build-debug.sh
   ```
   
   **发布版本**:
   ```bash
   ./build-release.sh
   ```
   
   **性能优化版本**:
   ```bash
   ./build-release-profile.sh
   ```

## 贡献流程

### 1. Fork 和 Clone

1. Fork 项目到您的 GitHub 账户
2. Clone 您的 fork:
   ```bash
   git clone https://github.com/YOUR_USERNAME/pbgz.git
   cd pbgz
   git remote add upstream https://github.com/PBGZip/pbgz.git
   ```

### 2. 创建分支

```bash
git checkout -b feature/your-feature-name
# 或
git checkout -b fix/your-bug-fix
```

### 3. 开发和测试

- 编写代码时请遵循项目的 [代码规范](#代码规范)
- 添加必要的测试用例
- 确保所有现有测试通过

```bash
# 运行测试
cd build
make test
# 或
ctest
```

### 4. 提交代码

```bash
git add .
git commit -m "feat: 添加新功能描述"
# 或
git commit -m "fix: 修复问题描述"
```

#### 提交信息格式

使用 [Conventional Commits](https://www.conventionalcommits.org/) 格式：

- `feat:` 新功能
- `fix:` Bug 修复
- `docs:` 文档更新
- `style:` 代码格式调整
- `refactor:` 代码重构
- `test:` 测试相关
- `perf:` 性能优化
- `ci:` CI/CD 相关

### 5. 推送和 Pull Request

```bash
git push origin feature/your-feature-name
```

然后在 GitHub 上创建 Pull Request。

## 代码规范

### C++ 代码风格

- **缩进**: 使用 4 个空格，不使用 Tab
- **命名**:
  - 类名: `PascalCase` (如 `PbgzEngine`)
  - 函数名: `snake_case` (如 `compress_file`)
  - 变量名: `snake_case` (如 `file_path`)
  - 常量: `kPascalCase` (如 `kMaxBufferSize`)
  - 宏定义: `UPPER_SNAKE_CASE` (如 `MAX_BUFFER_SIZE`)

- **注释**:
  - 头文件中的类和函数必须有详细注释
  - 复杂算法需要添加行内注释
  - 注释使用中文或英文，保持一致性

- **括号**: K&R 风格
  ```cpp
  if (condition) {
      // code
  } else {
      // code
  }
  ```

### 文件组织

- **头文件**: 只包含声明，避免 `using namespace`
- **源文件**: 包含实现，按逻辑分组
- **单元测试**: 放在 `test/testcase/` 目录

### 性能考虑

PBGZ 是性能敏感的项目，请特别注意：

- 避免不必要的内存分配
- 使用移动语义减少拷贝
- 关键路径避免虚函数调用
- 考虑缓存友好的数据结构

## 测试指南

### 测试类型

1. **单元测试**: 测试单个函数/类
2. **集成测试**: 测试组件间交互
3. **性能测试**: 验证压缩率和速度
4. **回归测试**: 确保修复不引入新问题

### 运行测试

```bash
# 构建测试
cd build
make test

# 运行特定测试
./test/pbgz_test_main --gtest_filter=YourTestCase

# 性能测试
./benchmark/pbgz_benchmark
```

### 添加测试

- 新功能必须包含测试用例
- Bug 修复需要添加回归测试
- 性能关键代码需要基准测试

## Bug 报告

报告 Bug 时请包含：

1. **环境信息**:
   - 操作系统和版本
   - CPU 架构
   - 编译器版本
   - PBGZ 版本

2. **复现步骤**:
   - 详细的操作步骤
   - 输入文件信息
   - 预期行为 vs 实际行为

3. **错误信息**:
   - 完整的错误日志
   - 栈跟踪（如果有）

4. **附加信息**:
   - 最小复现示例
   - 相关配置文件

## 功能请求

提出新功能时请描述：

1. **使用场景**: 为什么需要这个功能
2. **预期行为**: 功能应该如何工作
3. **替代方案**: 考虑过其他解决方案吗
4. **附加信息**: 相关链接、参考资料等

## 文档贡献

文档是项目的重要组成部分，我们欢迎：

- API 文档改进
- 使用示例补充
- 安装指南更新
- 性能调优建议
- 翻译工作

## 发布流程

项目的发布流程：

1. **版本规划**: 在 GitHub Issues 中讨论
2. **开发阶段**: 功能开发和测试
3. **代码审查**: 所有 PR 需要至少一个维护者审查
4. **集成测试**: 全面测试所有平台
5. **发布准备**: 更新版本号和 CHANGELOG
6. **正式发布**: 创建 Git 标签和 GitHub Release

## 社区准则

- 🤝 **尊重**: 尊重所有参与者，保持友好和专业
- 🧠 **学习**: 持续学习，乐于分享知识
- 🎯 **专注**: 保持对项目目标的专注
- 📈 **改进**: 不断改进代码和流程
- 🤲 **合作**: 鼓励协作，避免重复工作

## 获得帮助

如果您需要帮助或有疑问：

1. 查看 [文档](../docs/)
2. 搜索 [Issues](https://github.com/PBGZip/pbgz/issues)
3. 创建新的 Issue
4. 在 GitHub Discussions 中讨论

## 致谢

感谢所有为 PBGZ 项目做出贡献的开发者！您的贡献让这个项目变得更好。

---

再次感谢您的贡献！🎉
