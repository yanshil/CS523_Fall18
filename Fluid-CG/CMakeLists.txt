cmake_minimum_required (VERSION 2.8.11)

project(simpleFluid)

add_executable(simpleFluid ${SRC} "main.cpp")
add_executable(unitTest ${SRC} "unitTest.cpp")


## Build-in Libraries
find_package(OpenGL)
find_package(GLUT)
include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS})

target_link_libraries(simpleFluid
    PRIVATE
    ${OPENGL_gl_LIBRARY}
    ${OPENGL_glu_LIBRARY}
    ${GLUT_LIBRARIES})
