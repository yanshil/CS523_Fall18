add_executable(fluidsimulator
    TestSet.cpp
  )

  target_link_libraries(fluidsimulator
    PRIVATE
    nova::geometry
    nova::tools
    )
  
file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

# include_directories(include)
 
install(TARGETS fluidsimulator RUNTIME DESTINATION ${REL})

