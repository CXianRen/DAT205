cmake_minimum_required ( VERSION 3.8 ) # 3.8 is required for CUDA support

project ( particle_system LANGUAGES CXX CUDA ) # Enable CUDA as a language

set(CMAKE_CUDA_ARCHITECTURES 60 70 75 80 86)

# Find *all* shaders.
file(GLOB_RECURSE SHADERS
    "${CMAKE_CURRENT_SOURCE_DIR}/shaders/*.vert"
    "${CMAKE_CURRENT_SOURCE_DIR}/shaders/*.frag"
)

file(GLOB_RECURSE CPP_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/common/*.cpp"
)

file(GLOB_RECURSE CU_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/particle/*.cu"
)

# Separate filter for shaders.
source_group("Shaders" FILES ${SHADERS})

# Build and link executable.
add_executable ( ${PROJECT_NAME}
    main.cpp
    ${SHADERS}
    ${CPP_FILES}
    ${CU_FILES}
)

target_include_directories( ${PROJECT_NAME} 
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/common
    ${CMAKE_CURRENT_SOURCE_DIR}/particle
    ${CMAKE_CURRENT_SOURCE_DIR}/../external_src/Eigen
)

target_link_libraries ( ${PROJECT_NAME} labhelper )
config_build_output()
