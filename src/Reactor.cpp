#include "AdsorptionSimulator/Reactor.h"
#include "AdsorptionSimulator/DispersionModel.h"
#include "AdsorptionSimulator/ThomasAlgoritms.h"

#include "CoolProp.h"

#include <stdexcept>
#include <iostream>
#include <algorithm>  

Reactor::Reactor(Fluid &fluid) : fluid(fluid)
{
    dispersionModel = std::make_unique<DispersionModelConstant>(1e-5);
}

void Reactor::addLayer(const std::string &layerName)
{
    if (layers.find(layerName) != layers.end()) 
    {
        throw std::runtime_error("Layer with name '" + layerName + "' already exists.");
    }
    layers.emplace(layerName, PorousMedia(fluid, *this));
    
    updateLength();
    resizeData();
}

void Reactor::removeLayer(const std::string& layerName)
{
    auto it = layers.find(layerName);
    if (it == layers.end()) 
    {
        throw std::runtime_error("Layer with name '" + layerName + "' does not exist.");
    }
    layers.erase(it);

    updateLength();
    resizeData();
}

PorousMedia& Reactor::getLayer(const std::string& layerName)
{
    auto it = layers.find(layerName);
    if (it == layers.end()) 
    {
        throw std::runtime_error("Layer with name '" + layerName + "' not found.");
    }
    return it->second;
}

void Reactor::updateLength()
{
    length = 0.0;
    for (const auto& [layerName, layer]: layers) 
    {
        length += layer.getLayerLength(); 
    }
}

void Reactor::addComponent(const std::string& component)
{
    addIsothermModel(component);
    initialCondition.addComponent(component);
    resizeData();
}

void Reactor::removeComponent(const std::string& component)
{
    removeIsothermModel(component);
    initialCondition.removeComponent(component);
    resizeData();
}

void Reactor::addIsothermModel(const std::string& component)
{
    for (auto& [layerName, layer] : layers)
    {
        layer.addIsothermModel(component);
    }
}

void Reactor::removeIsothermModel(const std::string& component)
{
    for (auto& [layerName, layer] : layers)
    {
        layer.removeIsothermModel(component);
    }
}

void Reactor::resizeData()
{
    int totalCells = 2; 

    for (auto& [layerName, layer] : layers)
    {
        totalCells += layer.getNumberOfCells();
    }

    fluidData.resize(totalCells, fluid);
    fluidDataLastStep.resize(totalCells, fluid);
}

void Reactor::updateBoundaryCells()
{
    // Get inflow and outflow index. Index = 0 = BOTTOM, Index = N = TOP
    int inIndex = 0;
    int outIndex = 0;
    int di = 1;

    if (inflow.location == BoundaryConditionLocation::BOTTOM)
    {
        inIndex = 0;
        outIndex = int(fluidData.C.size()) - 1;
    }
    else
    {
        inIndex = int(fluidData.C.size()) - 1;
        outIndex = 0;
        di = -1;
    }

    // Set the inflow conditions
    fluidData.u[inIndex] = inflow.u_in;
    fluidData.T[inIndex] = 2 * inflow.T_in - fluidData.T[1];
    fluidData.P[inIndex] = fluidData.P[inIndex + di]; // No pressure drop in the first hald cell
    fluidData.C[inIndex] = fluidData.P[inIndex] / (8.314 * fluidData.T[inIndex]);
    for (const auto& component : fluid.components)
    {
        fluidData.yi[component][inIndex] = 2 * inflow.y_in[component] - fluidData.yi[component][inIndex + di];
        fluidData.Ci[component][inIndex] = 2 * inflow.y_in[component] * fluidData.P[inIndex + di] / (8.314 * inflow.T_in) - fluidData.Ci[component][inIndex + di];
    }
    
    // Stitch the interlayers
    for (auto itt = layers.begin(); itt != layers.end(); ++itt)
    {
        auto nextItt = std::next(itt);
        if (nextItt != layers.end())
        {
            auto& layer = itt->second;
            auto& nextLayer = nextItt->second;
            
            // TODO: what if flow is top down?
            fluidData.u[0] = fluidData.u[fluidData.u.size() - 1];
            fluidData.C[0] = fluidData.C[fluidData.C.size() - 2];
            fluidData.T[0] = fluidData.T[fluidData.T.size() - 2];
            fluidData.C[fluidData.C.size() - 1] = fluidData.C[1];
            fluidData.T[fluidData.T.size() - 1] = fluidData.T[1];
        
            for (const auto& component : fluid.components)
            {
                fluidData.Ci[component][0] = fluidData.Ci[component][fluidData.C.size() - 2];
                fluidData.yi[component][0] = fluidData.yi[component][fluidData.C.size() - 2];

                fluidData.Ci[component][fluidData.C.size() - 1] = fluidData.Ci[component][1];
                fluidData.yi[component][fluidData.C.size() - 1] = fluidData.yi[component][1];
            }
        }
    }

    // Outflow boundary conditions
    fluidData.T[outIndex] = fluidData.T[outIndex - di];
    fluidData.C[outIndex] = fluidData.C[outIndex - di];
    fluidData.P[outIndex] = outflow.P_out;
    for (const auto& component : fluid.components)
    {
        fluidData.Ci[component][outIndex] = fluidData.Ci[component][outIndex - di];
        fluidData.yi[component][outIndex] = fluidData.yi[component][outIndex - di];
    }
}

void Reactor::setDispersionModel(double dispersionCoefficient)
{
    dispersionModel = std::make_unique<DispersionModelConstant>(dispersionCoefficient);
}

void Reactor::initialise()
{    
    resizeData();

    for (int n = 0; n < fluidData.P.size(); ++n)
    {
        fluidData.P[n] = initialCondition.P0;
        fluidData.T[n] = initialCondition.T0;
        fluidData.u[n] = initialCondition.u0;
        fluidData.C[n] = initialCondition.P0 / (8.314 * initialCondition.T0);

        for (const auto& component : fluid.components)
        {
            fluidData.yi[component][n] = initialCondition.yi0[component];
            fluidData.Ci[component][n] = fluidData.C[n] * initialCondition.yi0[component];
        }
    }

    updateIsotherms();

    fluidData.qi = fluidData.qi_sat;
}

void Reactor::updateIsotherms()
{
    int startIndex = 0;
    for (auto& [layerName, layer] : layers)
    {
        int endIndex = startIndex + layer.getNumberOfCells() - 1;
        layer.updateIsotherms(fluidData, startIndex, endIndex); // Updates qi_sat for all components at each cell in the layer
        startIndex = endIndex++;
    }
}

void Reactor::updateSourceTerms()
{
    int startIndex = 0;
    for (auto& [layerName, layer] : layers)
    {
        int endIndex = startIndex + layer.getNumberOfCells() - 1;
        layer.updateSourceTerms(fluidData, startIndex, endIndex); // Updates qi_sat for all components at each cell in the layer
        startIndex = endIndex++;
    }
}

void Reactor::setBoundaryConditions(const Step& step)
{
    inflow = step.inflow;
    outflow = step.outflow;
}

void Reactor::updateConstants()
{
    fluidData.Cp.setConstant(0.);
	fluidData.k.setConstant(0.);
	fluidData.rho.setConstant(0.);
	fluidData.vis.setConstant(0.);
	fluidData.Yt.setConstant(0.);

	for (int n = 0; n < fluidData.C.size(); ++n)
	{
		for (const auto& component : fluid.components)
		{
			fluidData.yi[component][n] = fluidData.Ci[component][n] / fluidData.C[n];
			fluidData.Cp[n] += fluidData.yi[component][n] * CoolProp::PropsSI("CPMOLAR", "P", fluidData.P[n], "T", fluidData.T[n], component);
			fluidData.k[n] += fluidData.yi[component][n] * CoolProp::PropsSI("CONDUCTIVITY", "P", fluidData.P[n], "T", fluidData.T[n], component);
			fluidData.rho[n] += fluidData.Ci[component][n] * CoolProp::Props1SI("molemass", component) * 1e-3;
			fluidData.vis[n] += fluidData.yi[component][n] * CoolProp::PropsSI("V", "P", fluidData.P[n], "T", fluidData.T[n], component);
			fluidData.Yt[n] += fluidData.yi[component][n];
		}
	}
}

void Reactor::updateFlowrates()
{
    double area = wall.getArea();

	for (int n = 0; n < fluidData.u.size(); ++n)
	{
		fluidData.volumeFlow[n] = fluidData.u[n] * area;
		fluidData.massFlow[n] = fluidData.volumeFlow[n] * fluidData.rho[n];
		fluidData.molarFlow[n] = fluidData.volumeFlow[n] * fluidData.C[n];
	}
}

void Reactor::updateMoleFraction(double dt)
{
    Eigen::VectorXd Ae(fluidData.C.size());
	Eigen::VectorXd Aw(fluidData.C.size());
	Eigen::VectorXd Ap(fluidData.C.size());
	Eigen::VectorXd X(fluidData.C.size());

	for (const auto& component : fluid.components)
	{
        Ae.setZero();
        Aw.setZero();
        Ap.setZero();
        X.setZero();

        int startIndex = 0;
		for(auto& [layerName, layer] : layers)
        {
            int endIndex = layer.getNumberOfCells() - 1;
            layer.updateMoleFraction(startIndex, endIndex, component, dt, Ae, Aw, Ap, X);
            startIndex = endIndex + 1;
        }

		LinearSolver::TDMA(Ap, Ae, Aw, X, fluidData.Ci[component]);
	}
}

void Reactor::updatePressure(double dt)
{
    int totalCells = fluidData.P.size();
    Eigen::VectorXd Ae(totalCells);
	Eigen::VectorXd Aw(totalCells);
	Eigen::VectorXd Ap(totalCells);
	Eigen::VectorXd X(totalCells);
    
    Ae.setZero();
    Aw.setZero();
    Ap.setZero();
    X.setZero();

    int n = 0;
    Ap(n) = 1;
    X(n) = 1; // TODO: Inlet pressure

    int startIndex = 0;
    for(auto& [layerName, layer] : layers)
    {
        int endIndex = layer.getNumberOfCells() - 1;
        layer.updatePressure(startIndex, endIndex, dt, Ae, Aw, Ap, X);
        startIndex = endIndex + 1;
    }

    n = totalCells - 1;
    Ap(n) = 1;
    X(n) = 1; // TODO: Outlet pressure

    LinearSolver::TDMA(Ap, Ae, Aw, X, fluidData.P);
}

void Reactor::updateVelocity()
{
    int startIndex = 0;
    for(auto& [layerName, layer] : layers)
    {
        int endIndex = layer.getNumberOfCells() - 1;
        layer.updateVelocity(startIndex, endIndex);
        startIndex = endIndex + 1;
    }
}

void Reactor::integrate(double dt)
{
    updateBoundaryCells(); // Sets the inflow and outflow cells, as well as updating each layers interface

	// Explicit equations
	updateConstants();
	updateFlowrates();
	updateIsotherms();
	updateSourceTerms();
    dispersionModel->update(fluidData);

	// Implicit equations
	updateMoleFraction(dt);
	updatePressure(dt);
	updateVelocity();

	fluidDataLastStep = fluidData;
}

void InitialCondition::addComponent(const std::string& component)
{
    auto itt = yi0.find(component);

    if (itt == yi0.end()) // If it does not exist, add it
    {
        yi0[component] = 0.0;
    }
}

void InitialCondition::removeComponent(const std::string& component)
{
    auto itt = yi0.find(component);
    
    if (itt != yi0.end()) 
    {
        yi0.erase(itt);
    }
}
