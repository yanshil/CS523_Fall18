add_executable(fluidsimulator
    main.cpp
    FluidQuantity.cpp
)

target_link_libraries(fluidsimulator
    PRIVATE
    nova::tools
    nova::geometry

)

file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

include_directories(include)
 
install(TARGETS fluidsimulator RUNTIME DESTINATION ${REL})
