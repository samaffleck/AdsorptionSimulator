#pragma once

#include "Fluid.h"
#include "Reactor.h"
#include "Isotherm.h"

#include <Eigen/Dense>

#include <string>
#include <unordered_map>


struct Reactor;
struct IsothermModel;

struct FluidData
{
    // Component parameters
    std::unordered_map<std::string, Eigen::VectorXd> yi{};          // Mole fraction of component i [mol/mol]
    std::unordered_map<std::string, Eigen::VectorXd> Ci{};          // Molar concentration of component i [mol/m3]
    std::unordered_map<std::string, Eigen::VectorXd> qi{};          // Solid phase concentration [mol/kg]
    std::unordered_map<std::string, Eigen::VectorXd> qi_sat{};      // Equilibrium solid phase concentration [mol/kg]
    std::unordered_map<std::string, Eigen::VectorXd> Smi{};         // mass source [mol of component i/(kg solid - s)]
    std::unordered_map<std::string, Eigen::VectorXd> Sei{};         // energy source [J/(kg-s)]
    
    std::unordered_map<std::string, Eigen::VectorXd> ki{};          // Mass transfer coefficient [s]
    std::unordered_map<std::string, Eigen::VectorXd> Hads{};        // Heat of adsorption [J/mol]

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


struct PorousMedia
{
public:
    PorousMedia(Fluid& fluid, Reactor& reactor) : fluid(fluid), reactor(reactor) {}

    void addIsothermModel(const std::string& component);
    void removeIsothermModel(const std::string& component);

    void setIsothermModel(const std::string& component, const HenryIsothermParameters& ips);
    void setIsothermModel(const std::string& component, const LangmuirIsothermParameters& ips);
    void setIsothermModel(const std::string& component, const DualSiteLangmuirIsothermParameters& ips);
    void setIsothermModel(const std::string& component);
    
    void updateIsotherms();

    void setMassTransferCoefficient(const std::string& component, double ki);
    void setHeatOfAdsorption(const std::string& component, double Hads);

    void setLayerLength(double length);
    double getLayerLength() const;

    void setNumberOfCells(int numberOfCells);
    int getNumberOfCells() const;

    void resizeData();

    void integrate(double dt);

    double eb = 0.37;                   // Inter-particle bed voidage [m3/m3]
    double et = 0.65;                   // Intra-particle voidage [m3/m3]
    double density = 1200.0;            // Solid density [kg/m3]
    double Cp = 500.0;                  // Specific heat capacity [J/(kg-K)]
    double k = 0.5;                     // Thermal conductivity [W/(m-K)]
    double dp = 1e-3;                   // Particle diameter [m]
    double poreDiameter = 1e-6;         // Pore diameter [m]
    double tau = 0.1;                   // Macropore Void Fraction / tortuosity

    std::string name = "Default Solid";

    FluidData fluidData;
    FluidData fluidDataLastStep;

private:
    Fluid& fluid;
    Reactor& reactor;
    std::unordered_map<IsothermType, IsothermModel> isothermModels = {};  // Stores an isotherm model for each isotherm type
    int numberOfCells = 20;
    double L = 1.0;                             // Length of porous media domain [m]
    double dx = L / numberOfCells;

private:
    void updateConstants();
    void updateSourceTerms();
    void updateFlowrates();

    void updateMoleFraction(double dt);
    void updateVelocity();
};
