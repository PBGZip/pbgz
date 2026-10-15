# FC encoder component

# Source files
set(FC_SOURCES
    fc_model.cpp
    fc.cpp
    transform.cpp
    bsc_libsais_compat.cpp
)

# fc_model.cpp / bsc_libsais_compat.cpp are vendored/bridge code around
# third-party (libbsc / libsais); suppress warnings we don't need to fix.
set_source_files_properties(fc_model.cpp PROPERTIES COMPILE_FLAGS "-Wno-class-memaccess")
set_source_files_properties(bsc_libsais_compat.cpp PROPERTIES COMPILE_FLAGS "-Wno-unused-parameter")

# Create FC static library
add_library(pbgz_coder_fc STATIC ${FC_SOURCES})

# Set compile options
target_compile_options(pbgz_coder_fc PRIVATE -fPIC)

# The shared libsais copy (src/coder/libsais) is also linked by the BWT-CM
# coder in the parent directory.
target_link_libraries(pbgz_coder_fc PRIVATE pbgz_libsais)

# Expose header file directories
target_include_directories(pbgz_coder_fc PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)
