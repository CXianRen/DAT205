cmake_minimum_required ( VERSION 3.0.2 )

project ( OpenGL_SmokeSimulator )

# For windows we use our bundled binaries.
if(WIN32)
    # info
    message(STATUS "It's not supported to build this project on Windows.")
    return()
endif(WIN32)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(CMAKE_CONFIGURATION_TYPES "Debug;Release;RelWithDebInfo" CACHE STRING "" FORCE)

if(UNIX)
    set(OpenGL_GL_PREFERENCE "GLVND")
    list(APPEND CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/external/embree2")
    list(APPEND CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/external_src/glm")
    list(APPEND CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/external_linux/sdl2")
    list(APPEND CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/external_linux/glew")
    list(APPEND CMAKE_PREFIX_PATH "/chalmers/sw/sup64/cuda_toolkit-12.1.1")

endif(UNIX)
# EIGEN
set(EIGEN_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/external_src/Eigen")

# openMP
find_package(OpenMP REQUIRED)

# CUDA
find_package(CUDA REQUIRED)
set(CMAKE_CUDA_ARCHITECTURES 61 75)
enable_language(CUDA)
message(STATUS "CUDA_INCLUDE_DIRS: ${CUDA_INCLUDE_DIRS}")
message(STATUS "CUDA_LIBRARIES: ${CUDA_LIBRARIES}")
message(STATUS "CUDA_cusolver_LIBRARY: ${CUDA_cusolver_LIBRARY}")
message(STATUS "CUDA_cusparse_LIBRARY: ${CUDA_cusparse_LIBRARY}")
message(STATUS "CUDA_cublas_LIBRARY: ${CUDA_cublas_LIBRARY}")

# OpenGL
add_definitions(-DGLM_ENABLE_EXPERIMENTAL)
add_definitions(-D_CRT_SECURE_NO_WARNINGS)

set ( CMAKE_EXPORT_COMPILE_COMMANDS ON )
set ( CMAKE_CXX_STANDARD 14 )
if( NOT CMAKE_BUILD_TYPE )
  set ( CMAKE_BUILD_TYPE Release )
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3")

add_subdirectory ( labhelper )
add_subdirectory( particle_system )