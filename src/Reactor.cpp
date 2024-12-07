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
    layers.m_Layers.emplace_back(std::make_unique<PorousMedia>(layerName, &fluid, this));
    
    updateLength();
    resizeData();
}

void Reactor::removeLayer(const std::string& layerName)
{
    layers.remove(layerName);

    updateLength();
    resizeData();
}

PorousMedia& Reactor::getLayer(const std::string& layerName)
{
    return *layers.getLayer(layerName);
}

void Reactor::updateLength()
{
    length = 0.0;
    for (const auto& layer : layers.m_Layers) 
    {
        length += layer->getLayerLength(); 
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
    for (auto& layer : layers.m_Layers)
    {
        layer->addIsothermModel(component);
    }
}

void Reactor::removeIsothermModel(const std::string& component)
{
    for (auto& layer : layers.m_Layers)
    {
        layer->removeIsothermModel(component);
    }
}

void Reactor::resizeData()
{
    int totalCells = 2; 

    for (auto& layer : layers.m_Layers)
    {
        totalCells += layer->getNumberOfCells();
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
    if (layers.m_Layers.empty())
    {
        return false;
    }

    return true;
}

void Reactor::updateBoundaryCells()
{
    // Get inflow and outflow index. Index = 0 = BOTTOM, Index = N = TOP
    int inIndex = 0;
    PorousMedia* inLayer = layers.m_Layers[0].get();;

    int outIndex = 0;
    PorousMedia* outLayer = layers.m_Layers[layers.m_Layers.size() - 1].get();

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
        inLayer = layers.m_Layers[layers.m_Layers.size() - 1].get();
        outLayer = layers.m_Layers[0].get();
        di = -1;
    }

    // Set the inflow conditions
    fluidData.u[inIndex] = inflow.u_in;
    fluidData.T[inIndex] = 2 * inflow.T_in - fluidData.T[1];
    fluidData.P[inIndex] = fluidData.vis[inIndex] * inflow.u_in * inLayer->getCellWidth() / inLayer->permeability + fluidData.P[inIndex + di];
    fluidData.C[inIndex] = fluidData.P[inIndex] / (8.314 * fluidData.T[inIndex]);
    for (const auto& component : fluid.components)
    {
        //fluidData.yi[component][inIndex] = 2 * inflow.y_in[component] - fluidData.yi[component][inIndex + di];
        //fluidData.Ci[component][inIndex] = 2 * inflow.y_in[component] * fluidData.P[inIndex + di] / (8.314 * inflow.T_in) - fluidData.Ci[component][inIndex + di];
    }
    
    // Outflow boundary conditions
    fluidData.T[outIndex] = fluidData.T[outIndex - di];
    fluidData.C[outIndex] = fluidData.C[outIndex - di];
    fluidData.P[outIndex] = outflow.P_out;
    for (const auto& component : fluid.components)
    {
        //fluidData.Ci[component][outIndex] = fluidData.Ci[component][outIndex - di];
        //fluidData.yi[component][outIndex] = fluidData.yi[component][outIndex - di];
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
    for(auto& layer : layers.m_Layers)
    {
        int endIndex = startIndex + layer->getNumberOfCells() - 1;
        layer->initialise(startIndex, endIndex);
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
    for (auto& layer : layers.m_Layers)
    {
        layer->updateIsotherms(fluidData); // Updates qi_sat for all components at each cell in the layer
    }
}

void Reactor::updateSourceTerms()
{
    for (auto& layer : layers.m_Layers)
    {
        layer->updateSourceTerms(fluidData); // Updates Smi for all components at each cell in the layer
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

	for (int n = 0; n < fluidData.C.size(); ++n)
	{
		for (const auto& component : fluid.components)
		{
			fluidData.Cp[n] += fluidData.yi[component][n] * CoolProp::PropsSI("CPMOLAR", "P", fluidData.P[n], "T", fluidData.T[n], component);
			fluidData.k[n] += fluidData.yi[component][n] * CoolProp::PropsSI("CONDUCTIVITY", "P", fluidData.P[n], "T", fluidData.T[n], component);
			fluidData.rho[n] += fluidData.Ci[component][n] * CoolProp::Props1SI("molemass", component) * 1e-3;
			fluidData.vis[n] += fluidData.yi[component][n] * CoolProp::PropsSI("V", "P", fluidData.P[n], "T", fluidData.T[n], component);
		}
	}
}

void Reactor::updateMoleFraction()
{
    fluidData.Yt.setConstant(0.);

    for (int n = 0; n < fluidData.C.size(); ++n)
    {
        for (const auto& component : fluid.components)
        {
            fluidData.yi[component][n] = fluidData.Ci[component][n] / fluidData.C[n];
            fluidData.Yt[n] += fluidData.yi[component][n];
        }
    }

}

void Reactor::updateMolarConcentration()
{
    for (int n = 0; n < fluidData.C.size(); ++n)
    {
        fluidData.C[n] = fluidData.P[n] / (8.314 * fluidData.T[n]);
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

void Reactor::setComponentMolarConcentration(const std::string& component, double dt)
{
    int totalCells = Ap.size();

    Ae.setZero();
    Aw.setZero();
    Ap.setZero();
    X.setZero();

    double alpha = 1;
    double de = 1;
    double dw = 1;
    double dx0 = 1; // Prev layer dx
    double dx1 = 1; // Curr layer dx
    double dx2 = 1; // Next layer dx

    // Inlet ghost cell
    int n = 0;
    Aw(n) = 0;
    Ae(n) = 1;
    Ap(n) = 1;
    X(n) = 2 * inflow.y_in[component] * fluidData.P[n] / (8.314 * inflow.T_in);

    // First and last cell of each layer
    for (int l = 0; l < layers.m_Layers.size(); ++l)
    {
        PorousMedia* layer = layers.m_Layers[l].get();
        const PorousMedia* nextLayer = layer;
        const PorousMedia* previousLayer = layer;
        if (l < layers.m_Layers.size() - 1)
        {
            nextLayer = layers.m_Layers[l + 1].get();
        }
        if (l > 0)
        {
            previousLayer = layers.m_Layers[l - 1].get();
        }

        dx0 = previousLayer->getCellWidth();
        dx1 = layer->getCellWidth();
        dx2 = nextLayer->getCellWidth();

        // First cell of the layer
        n = n + 1;
        alpha = layer->et * dx1 / dt;
        de = layer->eb * 0.5 * (fluidData.Dl[n] + fluidData.Dl[n + 1]) / dx1;
        //dw = layer->eb * 2 * (fluidData.Dl[n]) / (dx0 + dx1);
        dw = 0;

        Aw[n] = -dw - std::max(fluidData.u[n - 1], 0.0);
        Ae[n] = -de - std::max(-fluidData.u[n], 0.0);
        Ap[n] = alpha + std::max(fluidData.u[n], 0.0) + std::max(-fluidData.u[n - 1], 0.0) + de + dw;
        X[n] = alpha * fluidDataLastStep.Ci[component][n] + layer->density * layer->eb * fluidData.Smi[component][n] * dx1;

        // Interior layer nodes
        layer->updateMoleFraction(component, dt, Ae, Aw, Ap, X);

        // Last cell of the layer
        n = n + layer->getNumberOfCells() - 1;
        alpha = layer->et * dx1 / dt;
        //de = layer->eb * 2.0 * fluidData.Dl[n] / (dx1 + dx2);
        de = 0;
        dw = layer->eb * 0.5 * (fluidData.Dl[n] + fluidData.Dl[n - 1]) / dx1;

        Aw[n] = -dw - std::max(fluidData.u[n - 1], 0.0);
        Ae[n] = -de - std::max(-fluidData.u[n], 0.0);
        Ap[n] = alpha + std::max(fluidData.u[n], 0.0) + std::max(-fluidData.u[n - 1], 0.0) + de + dw;
        X[n] = alpha * fluidDataLastStep.Ci[component][n] + layer->density * layer->eb * fluidData.Smi[component][n] * dx1;
    }

    // Outlet ghost cell
    n = totalCells - 1;
    Aw(n) = -1;
    Ae(n) = 0;
    Ap(n) = 1;
    X(n) = 0;
}

void Reactor::updateComponentMolarConcentration(double dt)
{
	for (const auto& component : fluid.components)
	{
        setComponentMolarConcentration(component, dt);
		LinearSolver::TDMA(Ap, Ae, Aw, X, fluidData.Ci[component]);
	}
}

void Reactor::getComponentMolarConcentrationError(double dt)
{
    for (const auto& component : fluid.components)
    {
        setComponentMolarConcentration(component, dt);
        LinearSolver::getError(Ap, Ae, Aw, X, fluidData.Ci[component]);
    }
}

void Reactor::updatePressure(double dt)
{    
    int totalCells = Ap.size();

    Ae.setZero();
    Aw.setZero();
    Ap.setZero();
    X.setZero();

    double alpha = 0;
    double R0 = 0;
    double R1 = 0;
    double R2 = 0;
    double R01 = 0;
    double R12 = 0;
    double dx0 = 0;
    double dx1 = 0;
    double dx2 = 0;

    // Inlet ghost cell
    int n = 0;
    Ap(n) = 1;
    X(n) = fluidData.P[n]; 

    for (int l = 0; l < layers.m_Layers.size(); ++l)
    {
        PorousMedia* layer = layers.m_Layers[l].get();
        const PorousMedia* nextLayer = layer;
        const PorousMedia* previousLayer = layer;
        if (l < layers.m_Layers.size() - 1)
        {
            nextLayer = layers.m_Layers[l + 1].get();
        }
        if (l > 0)
        {
            previousLayer = layers.m_Layers[l - 1].get();
        }

        // First cell of the layer
        n = n + 1;
        dx0 = previousLayer->getCellWidth();
        dx1 = layer->getCellWidth();
        dx2 = nextLayer->getCellWidth();
        R0 = previousLayer->permeability / (dx0 * fluidData.vis[n]);
        R1 = layer->permeability / (dx1 * fluidData.vis[n]);
        R01 = R0 * R1 / (R0 + R1);
        alpha = layer->eb * dx1 / dt;

        Aw(n) = -2 * fluidDataLastStep.P[n] * R01;
        Ae(n) = -fluidDataLastStep.P[n] * R1;
        Ap(n) = 2 * fluidDataLastStep.P[n] * R01 + fluidDataLastStep.P[n] * R1 + alpha;
        X(n) = alpha * fluidDataLastStep.P(n) + dx1 * 8.314 * fluidData.T[n] * fluidData.Sm[n];

        // Update interior nodes
        layer->updatePressure(dt, Ae, Aw, Ap, X);

        // Last cell of the layer
        n = n + layer->getNumberOfCells() - 1;
        R1 = layer->permeability / (dx1 * fluidData.vis[n]);
        R2 = nextLayer->permeability / (dx2 * fluidData.vis[n]);
        R12 = R1 * R2 / (R1 + R2);
        alpha = layer->eb * dx1 / dt;

        Aw(n) = -fluidDataLastStep.P[n] * R1;
        Ae(n) = -2 * fluidDataLastStep.P[n] * R12;
        Ap(n) = 2 * fluidDataLastStep.P[n] * R12  + fluidDataLastStep.P[n] * R1 + alpha;
        X(n) = alpha * fluidDataLastStep.P(n) + dx1 * 8.314 * fluidData.T[n] * fluidData.Sm[n];
    }

    n = totalCells - 1;
    Ap(n) = 1;
    X(n) = fluidData.P[n];

    LinearSolver::TDMA(Ap, Ae, Aw, X, fluidData.P);
}

void Reactor::updateVelocity()
{
    const PorousMedia* firstLayer = layers.m_Layers[0].get();
    const PorousMedia* secondLayer = layers.m_Layers[0].get();

    // Update the first cell in the first layer
    int n = 0;
    double R1 = firstLayer->permeability / (fluidData.vis[n] * firstLayer->getCellWidth());
    double R2 = 0;
    double R12 = 0;

    fluidData.u[n] = R1 * (fluidData.P[n] - fluidData.P[n + 1]);

    // Update the interior cells of each layer
    for (const auto& layer : layers.m_Layers)
    {
        layer->updateVelocity();
    }

    // Update the interior layer interfaces
    for (int l = 0; l < layers.m_Layers.size() - 1; ++l)
    {
        firstLayer = layers.m_Layers[l].get();
        secondLayer = layers.m_Layers[l + 1].get();

        n += firstLayer->getNumberOfCells();

        R1 = firstLayer->permeability / (fluidData.vis[n] * firstLayer->getCellWidth());
        R2 = secondLayer->permeability / (fluidData.vis[n] * secondLayer->getCellWidth());
        R12 = R1 * R2 / (R1 + R2);

        fluidData.u[n] = 2 * R12 * (fluidData.P[n] - fluidData.P[n + 1]);
    }

    // Update the final layer's last cell
    firstLayer = layers.m_Layers[layers.m_Layers.size() - 1].get();
    
    n += firstLayer->getNumberOfCells();
    
    R1 = firstLayer->permeability / (fluidData.vis[n] * firstLayer->getCellWidth());
    R2 = 0;
    R12 = 0;

    fluidData.u[n] = R1 * (fluidData.P[n] - fluidData.P[n + 1]);
}

void Reactor::integrate(double dt)
{
    
    for (int itt = 0; itt < 20; ++itt)
    {
        // Apply boundary conditions
        updateBoundaryCells();

        // Explicit equations
        //updateConstants();
        updateFlowrates();
        updateMolarConcentration();
        updateMoleFraction();
        //updateIsotherms();
        //updateSourceTerms();
        dispersionModel->update(fluidData);

        // Implicit equations
        updateComponentMolarConcentration(dt);
        updatePressure(dt);
        updateVelocity();
    }
    
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
