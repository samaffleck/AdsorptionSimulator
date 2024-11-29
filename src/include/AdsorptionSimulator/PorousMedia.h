#pragma once

#include "Fluid.h"
#include "Isotherm.h"

#include <Eigen/Dense>

#include <string>
#include <unordered_map>
#include <vector>


struct Reactor;
struct IsothermModel;
struct FluidData;


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
    
    void updateIsotherms(FluidData& fluidData);
    void updateSourceTerms(FluidData& fluidData);

    void updateMoleFraction(const std::string& component, double dt, Eigen::VectorXd& Ae, Eigen::VectorXd& Aw, Eigen::VectorXd& Ap, Eigen::VectorXd& X);
    void updatePressure(double dt, Eigen::VectorXd& Ae, Eigen::VectorXd& Aw, Eigen::VectorXd& Ap, Eigen::VectorXd& X);
    void updateVelocity();

    void setMassTransferCoefficient(const std::string& component, double ki);
    void setHeatOfAdsorption(const std::string& component, double Hads);

    void setLayerLength(double length);
    double getLayerLength() const;

    void setNumberOfCells(int numberOfCells);
    int getNumberOfCells() const;

    double getCellWidth() const;

    void setStartEndIndex(int startIndex, int endIndex);

    double eb = 0.37;                   // Inter-particle bed voidage [m3/m3]
    double et = 0.65;                   // Intra-particle voidage [m3/m3]
    double density = 1200.0;            // Solid density [kg/m3]
    double Cp = 500.0;                  // Specific heat capacity [J/(kg-K)]
    double k = 0.5;                     // Thermal conductivity [W/(m-K)]
    double dp = 1e-3;                   // Particle diameter [m]
    double poreDiameter = 1e-6;         // Pore diameter [m]
    double tau = 0.1;                   // Macropore Void Fraction / tortuosity
    double permeability = 4e-9;         // [1/m2]
    double compressibility = 9.87e-6;   // [1/Pa]

    std::string name = "Default Solid";

private:
    Fluid& fluid;
    Reactor& reactor;
    std::unordered_map<IsothermType, IsothermModel> isothermModels = {};  // Stores an isotherm model for each isotherm type
    int numberOfCells = 20;
    double L = 1.0;                             // Length of porous media domain [m]
    double dx = L / numberOfCells;
    std::unordered_map<std::string, double> massTransferCoefficient{};          // Mass transfer coefficient [s]
    std::unordered_map<std::string, double> heatOfAdsorption{};        // Heat of adsorption [J/mol]

    int m_StartIndex{};
    int m_EndIndex{};

};
