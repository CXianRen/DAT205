cmake_minimum_required ( VERSION 3.0.2 )

project ( simulator_test)


file(GLOB_RECURSE TCPP
    "${CMAKE_CURRENT_SOURCE_DIR}/../common/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/../particle/*.cpp"
)

# Build and link executable.
add_executable ( ${PROJECT_NAME}
    test.cpp
    ${TCPP}
)

target_include_directories( ${PROJECT_NAME} 
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../
)

target_link_libraries ( ${PROJECT_NAME} labhelper )
config_build_output()
