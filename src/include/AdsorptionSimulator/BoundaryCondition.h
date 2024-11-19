#pragma once

#include <vector>

enum class BoundaryConditionLocation
{
	TOP, BOTTOM
};


struct Inflow
{
	// Only u_in or Q_in or P_in needs to be specified
	double u_in = 0.1;
	double Q_in = 1.0;
	double P_in = 101325.0;

	double T_in = 293.0;
	std::vector<double> y_in{};

	BoundaryConditionLocation location = BoundaryConditionLocation::BOTTOM;
};


struct Outflow
{
	double P_out = 101325.0;

	BoundaryConditionLocation location = BoundaryConditionLocation::TOP;
};
