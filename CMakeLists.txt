set(ENABLE_CS523 ON CACHE BOOL "Enable building CS523 HW")

if(ENABLE_CS523)
  #add_subdirectory(HW1)
  #add_subdirectory(HW2)
  #add_subdirectory(HW3)
  add_subdirectory(FluidSimulator_Nova)
  add_subdirectory(FluidSimulator_Array_Based)
endif()

