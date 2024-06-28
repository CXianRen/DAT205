cmake_minimum_required ( VERSION 3.0.2 )

project (simulator_test  CXX CUDA )

file(GLOB_RECURSE TEST_CPP
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
)

# Build and link executable.
add_executable ( ${PROJECT_NAME}
    test.cpp
    ${SRC_CPP}
    ${CU_FILES} 
    ${TEST_CPP}
)

target_include_directories( ${PROJECT_NAME} 
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
)