cmake_minimum_required(VERSION 3.24)
set(CMAKE_INSTALL_PREFIX "/Users/eprager/local")

file(
GLOB HEADERS
RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
CONFIGURE_DEPENDS
"include/*.h"
)

project(mattresses CXX)

add_library(${PROJECT_NAME} INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/${HEADERS}")

set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 20)

target_include_directories(${PROJECT_NAME} INTERFACE
                           "${CMAKE_CURRENT_SOURCE_DIR}/${HEADERS}"
                           "/Users/eprager/local/include/"
                           )

target_link_libraries(${PROJECT_NAME} INTERFACE macros)

set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER ${HEADERS})

install(TARGETS ${PROJECT_NAME}
		LIBRARY DESTINATION lib
		PUBLIC_HEADER DESTINATION include
)

add_subdirectory("test")
