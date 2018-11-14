add_executable(fluidsimulator
    main.cpp
    #test.cpp
    Fluid.cpp
    #lodepng.cpp
    #UnitTest_Range_Iterator.cpp
    #unitTest_Grid.cpp
    #unitTest_computeVelocity.cpp
    #unitTest_indexAndOffset.cpp
    #unitTest_linearInterpolate.cpp
)

target_link_libraries(fluidsimulator
    PRIVATE
    nova::tools
    nova::geometry

)


file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
 
install(TARGETS fluidsimulator RUNTIME DESTINATION ${REL})

## Build-in Libraries
find_package(OpenGL)
find_package(GLUT)
include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS})

target_link_libraries(fluidsimulator
    PRIVATE
    ${OPENGL_gl_LIBRARY}
    ${OPENGL_glu_LIBRARY}
    ${GLUT_LIBRARIES})