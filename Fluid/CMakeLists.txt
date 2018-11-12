add_executable(fluidsimulator
    main.cpp
    Fluid.cpp
    # UnitTest_Range_Iterator.cpp
)

target_link_libraries(fluidsimulator
    PRIVATE
    nova::tools
    nova::geometry

)

file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
 
install(TARGETS fluidsimulator RUNTIME DESTINATION ${REL})
