# Add the static library
add_library(AdsorptionSimulator STATIC)

# Include public headers
target_include_directories(AdsorptionSimulator PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/external/Eigen)

# Add source files
target_sources(AdsorptionSimulator PRIVATE
    PorousMedia.cpp
    Fluid.cpp
    AdsorptionSystem.cpp
    Isotherm.cpp
    Reactor.cpp
    Wall.cpp
    BoundaryCondition.cpp
    DispersionModel.cpp
    Cycle.cpp
    Solver.cpp
    MassTransferModel.cpp
    ThomasAlgorithm.cpp
)

# Link CoolProp
target_link_libraries(AdsorptionSimulator PUBLIC CoolProp)

# Set C++ standard
target_compile_features(AdsorptionSimulator PUBLIC cxx_std_17)
