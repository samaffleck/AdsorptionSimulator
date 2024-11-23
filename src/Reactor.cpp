#include "AdsorptionSimulator/Reactor.h"

#include <stdexcept>
#include <iostream>

void Reactor::addLayer(const std::string& layerName)
{
    if (layers.find(layerName) != layers.end()) 
    {
        throw std::runtime_error("Layer with name '" + layerName + "' already exists.");
    }
    layers.try_emplace(layerName, PorousMedia(fluid, *this));

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

void Reactor::updateBoundaryConditions(const Step& step)
{
    inflow = step.inflow;
    outflow = step.outflow;
}

void Reactor::integrate(double dt)
{
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
