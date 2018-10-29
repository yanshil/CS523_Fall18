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

add_library(Contact Contact.cpp)
target_link_libraries(${PROJECT_NAME} Contact)

############################################################
#
#	Libraries
#
############################################################

### Eigen:

include_directories(include)


## Build-in Libraries
find_package(OpenGL)
find_package(GLUT)
include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS})

target_link_libraries(
    ${PROJECT_NAME}
    ${OPENGL_gl_LIBRARY}
    ${OPENGL_glu_LIBRARY}
    ${GLUT_LIBRARIES})
