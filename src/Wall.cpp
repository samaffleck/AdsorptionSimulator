#include "AdsorptionSimulator/Wall.h"

#define _USE_MATH_DEFINES
#include <math.h> 
#include <stdexcept>

Wall::Wall()
{
    updateThickness();
    updateArea();
}

void Wall::setInternalDiameter(double di)
{
    if (di <= 0) 
    {
        throw std::runtime_error("Internal diameter must be positive.");
    }
    d_i = di;
    updateThickness();
    updateArea();
}

// Get the internal diameter
double Wall::getInternalDiameter() const
{
    return d_i;
}

// Set the external diameter of the wall
void Wall::setExternalDiameter(double d_o_)
{
    if (d_o_ <= d_i) 
    {
        throw std::runtime_error("External diameter must be greater than internal diameter.");
    }
    d_o = d_o_;
    updateThickness();
    updateArea();
}

double Wall::getExternalDiameter() const
{
    return d_o;
}

void Wall::updateThickness()
{
    thickness = (d_o - d_i) / 2.0; 
}

void Wall::updateArea()
{
    area = M_PI_4 * d_i * d_i; 
}
