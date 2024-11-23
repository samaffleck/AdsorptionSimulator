#pragma once

#include <vector>
#include <string>

namespace Components 
{
    static std::vector<std::string> componentNames =
    {
        "Air", "Ammonia", "Argon", "Butane", "CarbonDioxide", "CarbonMonoxide",
        "Cyclopentane", "Deuterium", "D5", "D6", "DimethylEther", "Ethane", "Ethanol",
        "Ethylene", "Helium", "Hexane", "Hydrogen", "Isobutane", "Isopentane", "Krypton",
        "Methane", "Methanol", "Neon", "Nitrogen", "NitrousOxide", "Octane", "Oxygen",
        "Pentane", "Propane", "Propylene", "R11", "R12", "R123", "R1234yf", "R124",
        "R125", "R13", "R134a", "R14", "R141b", "R142b", "R143a", "R152a", "R227ea",
        "R23", "R236ea", "R236fa", "R245fa", "R32", "R404A", "R407C", "R410A",
        "R507A", "SES36", "SulfurHexafluoride", "Toluene", "Water", "Xenon"
    };
}

struct Fluid
{
public:
    void addComponent(const std::string& componentName);
    void removeComponent(const std::string& componentName);
    int getNumberOfComponents() { return int(components.size()); }

    std::vector<std::string> components = {};

private:
    bool isComponentInList(const std::string& componentName);
    bool isComponentValid(const std::string& componentName);

};
