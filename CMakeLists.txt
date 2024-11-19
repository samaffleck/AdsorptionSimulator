cmake_minimum_required(VERSION 3.15)
project(AdsorptionSimulator VERSION 1.0.0 LANGUAGES CXX)

# Enable testing
enable_testing()

# Add the CoolProp library first
set(COOLPROP_STATIC_LIBRARY ON CACHE BOOL "Build CoolProp as a static library")
add_subdirectory(external/CoolProp)

# Add the source library
add_subdirectory(src)

# Add the Natvis file for Visual Studio debugging
if(MSVC)
    set_property(TARGET AdsorptionSimulator PROPERTY VS_DEBUGGER_VISUALIZERS ${CMAKE_SOURCE_DIR}/external/Eigen/Eigen.natvis)
endif()

# Add tests
add_subdirectory(test)
