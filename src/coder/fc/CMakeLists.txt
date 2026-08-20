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

# libsais.c / fc_model.cpp are vendored third-party (libbsc); suppress warnings
# that come from upstream code we don't maintain.
set_source_files_properties(libsais.c PROPERTIES COMPILE_FLAGS "-Wno-unused-function")
set_source_files_properties(fc_model.cpp PROPERTIES COMPILE_FLAGS "-Wno-class-memaccess")

# Create FC static library
add_library(pbgz_coder_fc STATIC ${FC_SOURCES})

# Set compile options
target_compile_options(pbgz_coder_fc PRIVATE -fPIC)

# Expose header file directories
target_include_directories(pbgz_coder_fc PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)