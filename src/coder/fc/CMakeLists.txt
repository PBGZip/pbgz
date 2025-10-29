# FC encoder component

# Source files
set(FC_SOURCES
    fc_model.cpp
    fc.cpp
    transform.cpp
    libsais.c
)

# Set libsais.c as C language
set_source_files_properties(libsais.c PROPERTIES LANGUAGE C)

# Create FC static library
add_library(pbgz_coder_fc STATIC ${FC_SOURCES})

# Set compile options
target_compile_options(pbgz_coder_fc PRIVATE -fPIC)

# Expose header file directories
target_include_directories(pbgz_coder_fc PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)