cmake_minimum_required ( VERSION 3.0.2 )

project ( particle_system )

# Find *all* shaders.
file(GLOB_RECURSE SHADERS
    "${CMAKE_CURRENT_SOURCE_DIR}/shaders/*.vert"
    "${CMAKE_CURRENT_SOURCE_DIR}/shaders/*.frag"
)

file(GLOB_RECURSE MCPP
    "${CMAKE_CURRENT_SOURCE_DIR}/common/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/particle/*.cpp"
)

# Build and link executable.
add_executable ( ${PROJECT_NAME}
    main.cpp
    ${MCPP}
)

target_include_directories( ${PROJECT_NAME} 
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/common
    ${CMAKE_CURRENT_SOURCE_DIR}/particle
    ${CMAKE_CURRENT_SOURCE_DIR}/../external_src/Eigen
)

target_link_libraries ( ${PROJECT_NAME} labhelper )
config_build_output()

# add_subdirectory(test)