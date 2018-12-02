cmake_minimum_required (VERSION 2.8.11)

project(HW3)

add_executable(HW3 ${SRC} "unitTest.cpp")

target_link_libraries(HW3
    PRIVATE
    nova::tools
    nova::geometry
)


file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})