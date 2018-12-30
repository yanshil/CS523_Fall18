cmake_minimum_required(VERSION 2.6)
project(foo CXX)

############################################################
#
#	Build Options
#
############################################################

set(CMAKE_CXX_STANDARD 11)
add_executable(${PROJECT_NAME} main.cpp)

############################################################
#
#	Auxiliary Libraries
#
############################################################
add_library(RigidBody RigidBody.cpp)
target_link_libraries(${PROJECT_NAME} RigidBody)

############################################################
#
#	Libraries
#
############################################################

### Eigen:
include_directories(include)

### imgui
include_directories(imgui)
add_library(imgui ./imgui/imgui.cpp ./imgui/imgui_demo.cpp ./imgui/imgui_draw.cpp ./imgui/imgui_widgets.cpp ./imgui/example-glut/imgui_impl_freeglut.cpp ./imgui/example-glut/imgui_impl_opengl2.cpp)
target_link_libraries(${PROJECT_NAME} imgui)

## Build-in Libraries
find_package(OpenGL)
find_package(GLUT)
include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS})

target_link_libraries(
    ${PROJECT_NAME}
    ${OPENGL_gl_LIBRARY}
    ${OPENGL_glu_LIBRARY}
    ${GLUT_LIBRARIES})
