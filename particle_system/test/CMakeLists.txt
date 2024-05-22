cmake_minimum_required ( VERSION 3.0.2 )

project ( simulator_test)


file(GLOB_RECURSE SRC_CPP
    "${CMAKE_CURRENT_SOURCE_DIR}/../common/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/../particle/*.cpp"
)

file(GLOB_RECURSE TEST_CPP
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
)

# Build and link executable.
add_executable ( ${PROJECT_NAME}
    test.cpp
    ${SRC_CPP}
    ${TEST_CPP}
)

target_include_directories( ${PROJECT_NAME} 
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../
    ${CMAKE_CURRENT_SOURCE_DIR}/../common
    ${CMAKE_CURRENT_SOURCE_DIR}/../particle
)

target_link_libraries ( ${PROJECT_NAME} labhelper )
config_build_output()
