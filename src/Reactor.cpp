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
}

void Reactor::removeLayer(const std::string& layerName)
{
    auto it = layers.find(layerName);
    if (it == layers.end()) 
    {
        throw std::runtime_error("Layer with name '" + layerName + "' does not exist.");
    }
    layers.erase(it);
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
        

    }
}
