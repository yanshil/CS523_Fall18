add_executable(foo
    hw2.cpp
  )

  target_link_libraries(foo
    PRIVATE
    nova::geometry
    nova::tools
    )
  
file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})


include_directories(include)
 
install(TARGETS foo RUNTIME DESTINATION ${REL})

