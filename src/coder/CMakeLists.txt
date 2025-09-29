# 编码器组件CMake文件

# 源文件
set(CODER_SOURCES
    coder_json.cpp
    bcm_libsais.cpp
    coder.cpp
)

# 添加子目录
add_subdirectory(ppmd)
add_subdirectory(fc)

# 创建编码器库
add_library(pbgz_coder STATIC ${CODER_SOURCES})

# 添加头文件目录
target_include_directories(pbgz_coder PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/ppmd
    ${CMAKE_CURRENT_SOURCE_DIR}/fc
)

# 链接到子组件
target_link_libraries(pbgz_coder
    pbgz_coder_ppmd
    pbgz_coder_fc
)