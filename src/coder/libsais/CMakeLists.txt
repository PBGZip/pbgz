# Shared vendored upstream libsais 2.8.7 (Apache-2.0).
# See NOTICE.md and the file header in libsais.h for license details.
#
# This is the project's single libsais copy. It is shared by the FC coder
# (SEQ field / match stream, src/coder/fc) and by the BWT-CM coder
# (src/coder/coder_bwt_cm.h), which is why it lives in this common directory
# rather than under either coder.

set_source_files_properties(libsais.c PROPERTIES LANGUAGE C)

# Vendored third-party (libsais); suppress warnings from upstream code we don't
# maintain.
set_source_files_properties(libsais.c PROPERTIES COMPILE_FLAGS "-Wno-unused-function")

add_library(pbgz_libsais STATIC libsais.c)

target_compile_options(pbgz_libsais PRIVATE -fPIC)

target_include_directories(pbgz_libsais PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)
