#pragma once

#include "DispersionModel.h"
#include "PorousMedia.h"
#include "Wall.h"
//#include "Fluid.h"
#include "BoundaryCondition.h"
#include "Cycle.h"

#include <vector>
#include <unordered_map>
#include <memory>

struct Step;
struct Fluid;

struct InitialCondition
{
	double T0 = 293.0;
	double P0 = 101325.0;
	double u0 = 0.0;
	std::unordered_map<std::string, double> yi0{};

	void addComponent(const std::string& component);
	void removeComponent(const std::string& component);
};

struct Reactor
{
public:
	Reactor(Fluid& fluid);

	void addLayer(const std::string& layerName);
	void removeLayer(const std::string& layerName);
	PorousMedia& getLayer(const std::string& layerName);

	void addComponent(const std::string& component);
	void removeComponent(const std::string& component);

	void updateLength();

	void setDispersionModel(double dispersionCoefficient);

	void initialise(); // Sets the initial conditions to all cell values

	void setBoundaryConditions(const Step& step);
	
	void integrate(double dt);

	void resizeData();

	bool performChecks() const;
	
	Wall wall{};
	
	std::unique_ptr<IDispersionModel> dispersionModel;
	
	InitialCondition initialCondition{};

	FluidData fluidData;
    FluidData fluidDataLastStep;

private:
	Fluid& fluid;
	
	PorousLayers layers{};

	Inflow inflow{};
	Outflow outflow{};

	double length{};

	// Solver matrix elements
	Eigen::VectorXd Ae{};
	Eigen::VectorXd Aw{};
	Eigen::VectorXd Ap{};
	Eigen::VectorXd X{};

private:
	void addIsothermModel(const std::string& component);
	void removeIsothermModel(const std::string& component);
	void updateBoundaryCells();
	void updateIsotherms();
	void updateConstants();
	void updateFlowrates();
	void updateSourceTerms();
	
	void setComponentMolarConcentration(const std::string& component, double dt);
	void updateComponentMolarConcentration(double dt);
	void getComponentMolarConcentrationError(double dt);

	void updateMoleFraction();
	void updateMolarConcentration();
	void updatePressure(double dt);
	void updateVelocity();

};
