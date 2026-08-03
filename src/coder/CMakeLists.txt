# Coder component CMake file

# Source files
set(CODER_SOURCES
    coder_json.cpp
    bcm_libsais.cpp
    coder.cpp
    coder_fcv2.cpp
)

# Add subdirectories
add_subdirectory(ppmd)
add_subdirectory(fc)

# Create coder library
add_library(pbgz_coder STATIC ${CODER_SOURCES})

# Add header file directories
target_include_directories(pbgz_coder PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/ppmd
    ${CMAKE_CURRENT_SOURCE_DIR}/fc
)

# Link to subcomponents
target_link_libraries(pbgz_coder
    pbgz_coder_ppmd
    pbgz_coder_fc
)