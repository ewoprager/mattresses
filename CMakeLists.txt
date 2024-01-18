cmake_minimum_required(VERSION 3.24)
set(CMAKE_INSTALL_PREFIX "/Users/eprager/local")

file(
GLOB_RECURSE HEADERS
RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
"include/*.h"
)

file(
GLOB_RECURSE SOURCES
RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
"src/*.cpp"
)

project(mattresses CXX)

add_library(${PROJECT_NAME}
            "${CMAKE_CURRENT_SOURCE_DIR}/${SOURCES}"
            "${CMAKE_CURRENT_SOURCE_DIR}/${HEADERS}"
            )

target_include_directories(${PROJECT_NAME} PUBLIC
                           "${CMAKE_CURRENT_SOURCE_DIR}/${HEADERS}"
                           "/Users/eprager/local/include/"
                           )

target_link_libraries(${PROJECT_NAME} macros)

set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER ${HEADERS})

install(TARGETS ${PROJECT_NAME}
		LIBRARY DESTINATION lib
		PUBLIC_HEADER DESTINATION include
)
