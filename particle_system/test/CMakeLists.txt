cmake_minimum_required ( VERSION 3.0.2 )

set(CUDA_TOOLKIT_ROOT_DIR "/chalmers/sw/sup64/cuda_toolkit-12.1.1")
find_package(CUDA REQUIRED)

project (simulator_test  CXX CUDA )

set(CMAKE_CUDA_ARCHITECTURES 61 75)

file(GLOB_RECURSE SRC_CPP
    "${CMAKE_CURRENT_SOURCE_DIR}/../common/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/../particle/*.cpp"
)

file(GLOB_RECURSE TEST_CPP
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
)

file(GLOB_RECURSE CU_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/../particle/*.cu"
)

# Build and link executable.
add_executable ( ${PROJECT_NAME}
    test.cpp
    ${SRC_CPP}
    ${TEST_CPP}
    ${CU_FILES}
    
)

target_include_directories( ${PROJECT_NAME} 
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../
    ${CMAKE_CURRENT_SOURCE_DIR}/../common
    ${CMAKE_CURRENT_SOURCE_DIR}/../particle
    ${CUDA_INCLUDE_DIRS}
)

message(STATUS "test CUDA_INCLUDE_DIRS: ${CUDA_INCLUDE_DIRS}")

target_link_libraries ( 
    ${PROJECT_NAME} 
    labhelper
    ${CUDA_LIBRARIES}
    ${CUDA_cusolver_LIBRARY} 
    ${CUDA_cusparse_LIBRARY}
    ${CUDA_cublas_LIBRARY}
    )
config_build_output()
