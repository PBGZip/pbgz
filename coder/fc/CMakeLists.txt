# FC编码器组件

# 源文件
set(FC_SOURCES
    fc_model.cpp
    fc.cpp
    transform.cpp
    libsais.c
)

# 设置libsais.c为C语言
set_source_files_properties(libsais.c PROPERTIES LANGUAGE C)

# 创建FC静态库
add_library(pbgz_coder_fc STATIC ${FC_SOURCES})

# 设置编译选项
target_compile_options(pbgz_coder_fc PRIVATE -fPIC)

# 暴露头文件目录
target_include_directories(pbgz_coder_fc PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)