#include "AdsorptionSimulator/Reactor.h"
#include "AdsorptionSimulator/DispersionModel.h"

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
    for (auto& [layerName, layer] : layers)
    {
        layer.resizeData(); // Assuming PorousMedia has `resizeData`
    }
}

void Reactor::updateBoundaryCells()
{
    // Get inflow and outflow index. Index = 0 = BOTTOM, Index = N = TOP
    int inIndex = 0;
    PorousMedia* inLayer = &layers.begin()->second;

    int outIndex = 0;
    PorousMedia* outLayer = &layers.end()->second;

    int di = 1;

    if (inflow.location == BoundaryConditionLocation::BOTTOM)
    {
        inIndex = 0;
        outIndex = int(inLayer->fluidData.C.size()) - 1;
    }
    else
    {
        inLayer = &layers.end()->second;
        outLayer = &layers.begin()->second;
        inIndex = int(inLayer->fluidData.C.size()) - 1;
        outIndex = 0;
        di = -1;
    }

    // Set the inflow conditions
    inLayer->fluidData.u[inIndex] = inflow.u_in;
    inLayer->fluidData.T[inIndex] = 2 * inflow.T_in - inLayer->fluidData.T[1];
    inLayer->fluidData.P[inIndex] = inLayer->fluidData.P[inIndex + di]; // No pressure drop in the first hald cell
    inLayer->fluidData.C[inIndex] = inLayer->fluidData.P[inIndex] / (8.314 * inLayer->fluidData.T[inIndex]);
    for (const auto& component : fluid.components)
    {
        inLayer->fluidData.yi[component][inIndex] = 2 * inflow.y_in[component] - inLayer->fluidData.yi[component][inIndex + di];
        inLayer->fluidData.Ci[component][inIndex] = 2 * inflow.y_in[component] * inLayer->fluidData.P[inIndex + di] / (8.314 * inflow.T_in) - inLayer->fluidData.Ci[component][inIndex + di];
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
            nextLayer.fluidData.u[0] = layer.fluidData.u[layer.fluidData.u.size() - 1];

            nextLayer.fluidData.C[0] = layer.fluidData.C[layer.fluidData.C.size() - 2];
            nextLayer.fluidData.T[0] = layer.fluidData.T[layer.fluidData.T.size() - 2];
        
            layer.fluidData.C[layer.fluidData.C.size() - 1] = nextLayer.fluidData.C[1];
            layer.fluidData.T[layer.fluidData.T.size() - 1] = nextLayer.fluidData.T[1];
        
            for (const auto& component : fluid.components)
            {
                nextLayer.fluidData.Ci[component][0] = layer.fluidData.Ci[component][layer.fluidData.C.size() - 2];
                nextLayer.fluidData.yi[component][0] = layer.fluidData.yi[component][layer.fluidData.C.size() - 2];

                layer.fluidData.Ci[component][layer.fluidData.C.size() - 1] = nextLayer.fluidData.Ci[component][1];
                layer.fluidData.yi[component][layer.fluidData.C.size() - 1] = nextLayer.fluidData.yi[component][1];
            }
        }
    }

    // Outflow boundary conditions
    outLayer->fluidData.T[outIndex] = outLayer->fluidData.T[outIndex - di];
    outLayer->fluidData.C[outIndex] = outLayer->fluidData.C[outIndex - di];
    outLayer->fluidData.P[outIndex] = outflow.P_out;
    for (const auto& component : fluid.components)
    {
        outLayer->fluidData.Ci[component][outIndex] = outLayer->fluidData.Ci[component][outIndex - di];
        outLayer->fluidData.yi[component][outIndex] = outLayer->fluidData.yi[component][outIndex - di];
    }
}

void Reactor::setDispersionModel(double dispersionCoefficient)
{
    dispersionModel = std::make_unique<DispersionModelConstant>(dispersionCoefficient);
}

void Reactor::initialise()
{
    for (auto& [layerName, layer] : layers)
    {
        for (int n = 0; n < layer.fluidData.P.size(); ++n)
        {
            layer.fluidData.P[n] = initialCondition.P0;
            layer.fluidData.T[n] = initialCondition.T0;
            layer.fluidData.u[n] = initialCondition.u0;
            layer.fluidData.C[n] = initialCondition.P0 / (8.314 * initialCondition.T0);

            for (const auto& component : fluid.components)
            {
                layer.fluidData.yi[component][n] = initialCondition.yi0[component];
                layer.fluidData.Ci[component][n] = layer.fluidData.C[n] * initialCondition.yi0[component];
            }
        }
        layer.updateIsotherms(); // Updates qi_sat for all components at each cell in the layer
        layer.fluidData.qi = layer.fluidData.qi_sat;
    }
}

void Reactor::setBoundaryConditions(const Step& step)
{
    inflow = step.inflow;
    outflow = step.outflow;
}

void Reactor::integrate(double dt)
{
    updateBoundaryCells(); // Sets the inflow and outflow cells, as well as updating each layers interface

    for (auto& [layerName, layer] : layers) 
    {
        dispersionModel->update(layer);
        layer.integrate(dt);
    }
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
