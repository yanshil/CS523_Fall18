set(ENABLE_CS523 ON CACHE BOOL "Enable building CS523 HW")

if(ENABLE_CS523)
  #add_subdirectory(HW1)
  #add_subdirectory(HW2)
  #add_subdirectory(Fluid)
  #add_subdirectory(Fluid2)
  #add_subdirectory(0-test)
  add_subdirectory(0-CG)
  # add_subdirectory(HW3)
endif()

