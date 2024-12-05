#pragma once

#include "Eigen/Dense"

#include <vector>
#include <string>
#include <unordered_map>


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
    bool isComponentInList(const std::string& componentName) const;
    bool isComponentValid(const std::string& componentName) const;

};

struct FluidData
{
    // Component parameters
    std::unordered_map<std::string, Eigen::VectorXd> yi{};          // Mole fraction of component i [mol/mol]
    std::unordered_map<std::string, Eigen::VectorXd> Ci{};          // Molar concentration of component i [mol/m3]
    std::unordered_map<std::string, Eigen::VectorXd> qi{};          // Solid phase concentration [mol/kg]
    std::unordered_map<std::string, Eigen::VectorXd> qi_sat{};      // Equilibrium solid phase concentration [mol/kg]
    std::unordered_map<std::string, Eigen::VectorXd> Smi{};         // mass source [mol of component i/(kg solid - s)]
    std::unordered_map<std::string, Eigen::VectorXd> Sei{};         // energy source [J/(kg-s)]
    
    // Overall fluid parameters
    Eigen::VectorXd rho{};                      // Gas density [kg / m3]
    Eigen::VectorXd C{};                        // Total molar concentration [mol/m3]
    Eigen::VectorXd T{};                        // Temperature [K]
    Eigen::VectorXd P{};                        // Pressure [Pa]
    Eigen::VectorXd Sm{};                       // Total mass source [mol/(kg-s)]
    Eigen::VectorXd Se{};                       // Total energy source [J/(kg-s)]
    Eigen::VectorXd Yt{};                       // Sum of mole fractions [mol/mol]
    Eigen::VectorXd vis{};                      // Fluid viscisity [Pa-s]
    Eigen::VectorXd Cp{};                       // Fluid heat capcity [J/(mol-K)]
    Eigen::VectorXd Cp_ads{};                   // Adsorbent phase heat capcity [J/(mol-K)]
    Eigen::VectorXd k{};                        // Fluid thermal conductivity [W/(m-K)]
    Eigen::VectorXd Dl{};                       // Axial dispersion coefficient [m/s2]
    Eigen::VectorXd Lam{};                      // Effective axial heat dispersion [m/s2]
    Eigen::VectorXd h_gw{};                     // Heat transfer coefficient between the gas and wall [W/(m2-K)]
    Eigen::VectorXd h_wa{};                     // Heat transfer coefficient between the wall and ambient [W/(m2-K)]
    Eigen::VectorXd Re{};                       // Reynolds number [-]
    Eigen::VectorXd Sc{};                       // Schmidt number [-]
    Eigen::VectorXd Pr{};                       // Prandtl number [-]

    Eigen::VectorXd u{};                        // Fluid velocity [m/s] -> Staggered grid
    Eigen::VectorXd massFlow{};                 // Mass flow rate [kg/s] -> Staggered grid
    Eigen::VectorXd molarFlow{};                // Molar flow rate [mol/s] -> Staggered grid
    Eigen::VectorXd volumeFlow{};               // volumetric flow rate [m3/s] -> Staggered grid

    void resize(int sizeOfVectors, const Fluid& fluid);

};

struct FluidDataLogger
{
    std::vector<FluidData> timeSeriesData{}; 

    void saveDataToCSV(const std::string& fileDirectory) const;

private:
    void logVariableToFile(const Eigen::VectorXd& variable, std::ofstream& variableFile) const;

};
