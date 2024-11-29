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

    Ae.resize(totalCells);
    Aw.resize(totalCells);
    Ap.resize(totalCells);
    X.resize(totalCells);

    fluidData.resize(totalCells, fluid);
    fluidDataLastStep.resize(totalCells, fluid);
}

bool Reactor::performChecks() const
{
    if (layers.empty())
    {
        return false;
    }

    return true;
}

void Reactor::updateBoundaryCells()
{
    // Get inflow and outflow index. Index = 0 = BOTTOM, Index = N = TOP
    int inIndex = 0;
    PorousMedia* inLayer = &layers.begin()->second;

    int outIndex = 0;
    PorousMedia* outLayer = &std::prev(layers.end())->second;

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
        inLayer = &layers.end()->second;
        outLayer = &layers.begin()->second;
        di = -1;
    }

    // Set the inflow conditions
    fluidData.u[inIndex] = inflow.u_in;
    fluidData.T[inIndex] = 2 * inflow.T_in - fluidData.T[1];
    fluidData.P[inIndex] = fluidData.vis[inIndex] * inflow.u_in * inLayer->getCellWidth() / inLayer->permeability + fluidData.P[inIndex + di];
    fluidData.C[inIndex] = fluidData.P[inIndex] / (8.314 * fluidData.T[inIndex]);
    for (const auto& component : fluid.components)
    {
        fluidData.yi[component][inIndex] = 2 * inflow.y_in[component] - fluidData.yi[component][inIndex + di];
        fluidData.Ci[component][inIndex] = 2 * inflow.y_in[component] * fluidData.P[inIndex + di] / (8.314 * inflow.T_in) - fluidData.Ci[component][inIndex + di];
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

    // Set the start and end index for all of the layers
    int startIndex = 1;
    for(auto& [layerName, layer] : layers)
    {
        int endIndex = startIndex + layer.getNumberOfCells() - 1;
        layer.initialise(startIndex, endIndex);
        startIndex = endIndex + 1;
    }

    for (int n = 0; n < fluidData.u.size(); ++n)
    {
        fluidData.u[n] = initialCondition.u0;
    }

    for (int n = 0; n < fluidData.P.size(); ++n)
    {
        fluidData.P[n] = initialCondition.P0;
        fluidData.T[n] = initialCondition.T0;
        fluidData.C[n] = initialCondition.P0 / (8.314 * initialCondition.T0);

        for (const auto& component : fluid.components)
        {
            fluidData.yi[component][n] = initialCondition.yi0[component];
            fluidData.Ci[component][n] = fluidData.C[n] * initialCondition.yi0[component];
        }
    }

    updateIsotherms();
    updateConstants();

    fluidData.qi = fluidData.qi_sat;
    fluidDataLastStep = fluidData;
}

void Reactor::updateIsotherms()
{
    for (auto& [layerName, layer] : layers)
    {
        layer.updateIsotherms(fluidData); // Updates qi_sat for all components at each cell in the layer
    }
}

void Reactor::updateSourceTerms()
{
    for (auto& [layerName, layer] : layers)
    {
        layer.updateSourceTerms(fluidData); // Updates Smi for all components at each cell in the layer
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
    int totalCells = Ap.size();

	for (const auto& component : fluid.components)
	{
        Ae.setZero();
        Aw.setZero();
        Ap.setZero();
        X.setZero();

        // Inlet ghost cell
        int n = 0;
        Ap(n) = 1;
        X(n) = fluidData.Ci[component][n];

        // Internal nodes
        for(auto& [layerName, layer] : layers)
        {
            layer.updateMoleFraction(component, dt, Ae, Aw, Ap, X);
        }

        n = totalCells - 1;
        Ap(n) = 1;
        X(n) = fluidData.Ci[component][n];

		LinearSolver::TDMA(Ap, Ae, Aw, X, fluidData.Ci[component]);
	}
}

void Reactor::updatePressure(double dt)
{    
    int totalCells = Ap.size();

    Ae.setZero();
    Aw.setZero();
    Ap.setZero();
    X.setZero();

    int n = 0;
    Ap(n) = 1;
    X(n) = fluidData.P[n]; 

    for(auto& [layerName, layer] : layers)
    {
        layer.updatePressure(dt, Ae, Aw, Ap, X);
    }

    n = totalCells - 1;
    Ap(n) = 1;
    X(n) = fluidData.P[n];

    LinearSolver::TDMA(Ap, Ae, Aw, X, fluidData.P);
}

void Reactor::updateVelocity()
{
    for(auto& [layerName, layer] : layers)
    {
        layer.updateVelocity();
    }
}

void Reactor::integrate(double dt)
{
    // Apply boundary conditions
    updateBoundaryCells(); 

	// Explicit equations
	//updateConstants();
	updateFlowrates();
	//updateIsotherms();
	//updateSourceTerms();
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
