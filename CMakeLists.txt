cmake_minimum_required(VERSION 3.24)

add_subdirectory("External/macros")

file(GLOB_RECURSE HEADERS
     RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
     "include/*.h"
)

file(GLOB_RECURSE SOURCES
     RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
     "src/*.cpp"
)

project(mattresses CXX)

add_library(${PROJECT_NAME}
            ${SOURCES}
            ${HEADERS}
)

target_include_directories(${PROJECT_NAME} PUBLIC
                           "include/"
)

target_link_libraries(${PROJECT_NAME} macros)

set_target_properties(${PROJECT_NAME} PROPERTIES
                      PUBLIC_HEADER ${HEADERS}
)

install(TARGETS ${PROJECT_NAME}
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include
)
