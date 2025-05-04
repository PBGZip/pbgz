# PPMD编码器组件

# 查找所有C源文件
file(GLOB PPMD_SOURCES "Ppmd*.c")

# 使用C编译器
set_source_files_properties(${PPMD_SOURCES} PROPERTIES LANGUAGE C)

# 创建PPMD静态库
add_library(pbgz_coder_ppmd STATIC ${PPMD_SOURCES})

# 设置C编译选项
target_compile_options(pbgz_coder_ppmd PRIVATE -fPIC)

# 暴露头文件目录
target_include_directories(pbgz_coder_ppmd PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)