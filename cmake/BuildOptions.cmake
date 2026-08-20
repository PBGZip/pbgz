# Build configuration options - Debug/Release differentiation with multi-compiler support

# Detect operating system
if(APPLE)
    # macOS specific definitions
    add_definitions(-DMACOSX)

    # Exclude sys/io.h which is Linux-specific
    add_definitions(-DNO_SYS_IO_H)
elseif(UNIX AND NOT APPLE)
    # Linux specific definitions
    add_definitions(-DLINUX)
endif()

# Detect compiler type
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(PBGZ_COMPILER "gcc")
    message(STATUS "Detected GNU compiler")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(PBGZ_COMPILER "clang")
    message(STATUS "Detected Clang compiler")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    set(PBGZ_COMPILER "icc")
    message(STATUS "Detected Intel compiler")
else()
    set(PBGZ_COMPILER "unknown")
    message(STATUS "Using unknown compiler")
endif()


# Detect processor architecture
if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64")
    set(PBGZ_ARCH "arm64")
    message(STATUS "Detected ARM64 architecture")
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|amd64|AMD64")
    set(PBGZ_ARCH "x86_64")
    message(STATUS "Detected x86_64 architecture")
else()
    set(PBGZ_ARCH "unknown")
    message(STATUS "Detected unknown architecture")
endif()

# Platform detection
if(APPLE)
    set(PBGZ_MACOS TRUE)
    message(STATUS "Building for macOS platform")

    # Platform-specific checks for macOS
    if(EXISTS "/usr/include/sys/io.h")
        add_definitions(-DHAVE_SYS_IO_H)
    endif()
elseif(UNIX AND NOT APPLE)
    set(PBGZ_LINUX TRUE)
    message(STATUS "Building for Linux platform")
    add_definitions(-DHAVE_SYS_IO_H) # Linux typically has this header
elseif(WIN32)
    set(PBGZ_WINDOWS TRUE)
    message(STATUS "Building for Windows platform")
endif()

# General compilation flags (common to all compilers)
set(PBGZ_COMMON_FLAGS "-fvisibility=hidden")

# Architecture-specific flags
if(PBGZ_ARCH STREQUAL "x86_64")
    set(PBGZ_ARCH_FLAGS "-mbmi2 -msse4.2 -mmovbe -mavx2")
    add_definitions(-D__BMI2__ -D__SSE4_2__)
endif()

# Compiler-specific flags
if(PBGZ_COMPILER STREQUAL "gcc")
    # GCC specific flags
    set(COMPILER_FLAGS "-Wall -Wextra -fpermissive")
    set(COMPILER_FLAGS_C "-Wall -Wextra")
    set(OPTIMIZATION_DEBUG "-O0")
    set(OPTIMIZATION_RELEASE "-O3")

elseif(PBGZ_COMPILER STREQUAL "clang")
    # Clang specific flags
    set(COMPILER_FLAGS "-Wall -Wextra -fpermissive -Wno-deprecated-register")
    set(COMPILER_FLAGS_C "-Wall -Wextra -Wno-deprecated-register")
    set(OPTIMIZATION_DEBUG "-O0")
    set(OPTIMIZATION_RELEASE "-O3")

elseif(PBGZ_COMPILER STREQUAL "icc")
    # Intel compiler specific flags
    set(COMPILER_FLAGS "-Wall -fpermissive")
    set(COMPILER_FLAGS_C "-Wall")
    set(OPTIMIZATION_DEBUG "-O0")
    set(OPTIMIZATION_RELEASE "-O3 -ipo") # Intel-specific optimization
endif()

# Debug/Release specific flags
set(CMAKE_CXX_FLAGS_DEBUG "${OPTIMIZATION_DEBUG} -g -DDEBUG -DFARMHASH_DEBUG=0")
set(CMAKE_C_FLAGS_DEBUG "${OPTIMIZATION_DEBUG} -g -DDEBUG -DFARMHASH_DEBUG=0")

set(CMAKE_CXX_FLAGS_RELEASE "${OPTIMIZATION_RELEASE} -DNDEBUG")
set(CMAKE_C_FLAGS_RELEASE "${OPTIMIZATION_RELEASE} -DNDEBUG")

# Release-Profile build type: optimized with debug symbols
set(CMAKE_CXX_FLAGS_RELEASEPROFILE "${OPTIMIZATION_RELEASE} -g -DNDEBUG")
set(CMAKE_C_FLAGS_RELEASEPROFILE "${OPTIMIZATION_RELEASE} -g -DNDEBUG")

# Apply common flags, architecture flags, and compiler flags.
# -fpermissive is C++-only; C files must not receive it.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PBGZ_COMMON_FLAGS} ${PBGZ_ARCH_FLAGS} ${COMPILER_FLAGS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${PBGZ_COMMON_FLAGS} ${PBGZ_ARCH_FLAGS} ${COMPILER_FLAGS_C}")

# Address Sanitizer option
option(ENABLE_SANITIZER "Enable Address Sanitizer" OFF)

if(ENABLE_SANITIZER)
    if(NOT PBGZ_COMPILER STREQUAL "icc") # Intel compiler might have different sanitizer syntax
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
    else()
        message(WARNING "Address Sanitizer may not be fully supported with Intel compiler")
    endif()
endif()

# Security hardening (similar to what's in the configure.ac files)
if(UNIX)
    set(SECURITY_FLAGS "-fstack-protector -D_FORTIFY_SOURCE=2")

    if(NOT APPLE)
        # Linux-specific security flags
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-z,noexecstack -Wl,-z,relro -Wl,-z,now")
    endif()

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SECURITY_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${SECURITY_FLAGS}")
endif()

# Display current configuration summary
message(STATUS "============ Build Information ============")
message(STATUS "Architecture: ${PBGZ_ARCH}")
message(STATUS "Compiler: ${PBGZ_COMPILER}")
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Address Sanitizer: ${ENABLE_SANITIZER}")
message(STATUS "Third-party include paths added to compilation")
message(STATUS "C flags: ${CMAKE_C_FLAGS}")
message(STATUS "C++ flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "C++ Debug flags: ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "C++ Release flags: ${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "C++ Release-Profile flags: ${CMAKE_CXX_FLAGS_RELEASEPROFILE}")
message(STATUS "Linker flags: ${CMAKE_EXE_LINKER_FLAGS}")
message(STATUS "=========================================")
