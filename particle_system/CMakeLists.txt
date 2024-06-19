cmake_minimum_required ( VERSION 3.0.2 )

set(CUDA_TOOLKIT_ROOT_DIR "/chalmers/sw/sup64/cuda_toolkit-12.1.1")
find_package(CUDA REQUIRED)

project ( particle_system LANGUAGES CXX CUDA )

set(CMAKE_CUDA_ARCHITECTURES 61 75)

find_package(OpenMP REQUIRED)

file(GLOB_RECURSE CPP_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/common/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/particle/*.cpp"
)

file(GLOB_RECURSE CU_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/particle/*.cu"
)

# Build and link executable.
add_executable ( ${PROJECT_NAME}
    main.cpp
    ${CPP_FILES}
    ${CU_FILES}
)

target_include_directories( ${PROJECT_NAME} 
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/common
    ${CMAKE_CURRENT_SOURCE_DIR}/particle
    ${CMAKE_CURRENT_SOURCE_DIR}/../external_src/Eigen
    ${CUDA_INCLUDE_DIRS}
)

message(STATUS "CUDA_INCLUDE_DIRS: ${CUDA_INCLUDE_DIRS}")
message(STATUS "CUDA_LIBRARIES: ${CUDA_LIBRARIES}")
message(STATUS "CUDA_cusolver_LIBRARY: ${CUDA_cusolver_LIBRARY}")
message(STATUS "CUDA_cusparse_LIBRARY: ${CUDA_cusparse_LIBRARY}")
message(STATUS "CUDA_cublas_LIBRARY: ${CUDA_cublas_LIBRARY}")

target_link_libraries ( 
    ${PROJECT_NAME} 
    labhelper 
    ${CUDA_LIBRARIES}
    ${CUDA_cusolver_LIBRARY} 
    ${CUDA_cusparse_LIBRARY}
    ${CUDA_cublas_LIBRARY}
    OpenMP::OpenMP_CXX
)
config_build_output()

# add_subdirectory(test)