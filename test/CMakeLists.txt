cmake_minimum_required(VERSION 3.24)

project(mattresses_test CXX)

add_executable(mattresses_test main.cpp)

set_property(TARGET mattresses_test PROPERTY CXX_STANDARD 20)

target_include_directories(mattresses_test PUBLIC
                           "/Users/eprager/local/include/"
						   "${CMAKE_CURRENT_SOURCE_DIR}/../include/"
                           )
