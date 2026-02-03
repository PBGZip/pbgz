# PPMD encoder component

# Find all C source files
file(GLOB PPMD_SOURCES "Ppmd*.c")

# Use C compiler
set_source_files_properties(${PPMD_SOURCES} PROPERTIES LANGUAGE C)

# Create PPMD static library
add_library(pbgz_coder_ppmd STATIC ${PPMD_SOURCES})

# Set C compile options
target_compile_options(pbgz_coder_ppmd PRIVATE -fPIC)

# Expose header file directories
target_include_directories(pbgz_coder_ppmd PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)