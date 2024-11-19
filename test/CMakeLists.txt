# Add the test executable
add_executable(test_fluid_solver test_fluid_solver.cpp)

# Link the static library
target_link_libraries(test_fluid_solver PRIVATE AdsorptionSimulator)

# Add the test to CTest
add_test(NAME MyLibraryTest COMMAND test_fluid_solver)
