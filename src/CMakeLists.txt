cmake_minimum_required ( VERSION 3.0.2 )

project ( particle_system LANGUAGES CXX CUDA )

file(GLOB_RECURSE CPP_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/common/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/particle/*.cpp"
)

file(GLOB_RECURSE CU_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/particle/*.cu"
)

include_directories(${EIGEN_INCLUDE_DIR})
include_directories(${CUDA_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/common)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/particle)

link_libraries(
    ${CUDA_LIBRARIES}
    labhelper 
    ${CUDA_LIBRARIES}
    ${CUDA_cusolver_LIBRARY} 
    ${CUDA_cusparse_LIBRARY}
    ${CUDA_cublas_LIBRARY}
    OpenMP::OpenMP_CXX
)

# Build and link executable.
add_executable ( ${PROJECT_NAME}
    main.cpp
    ${CPP_FILES}
    ${CU_FILES}
)

target_include_directories( 
    ${PROJECT_NAME} 
    PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}
)


# target_link_libraries ( 
#     ${PROJECT_NAME} 
# )

add_subdirectory(test)