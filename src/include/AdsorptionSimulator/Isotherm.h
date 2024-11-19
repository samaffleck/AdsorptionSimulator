#pragma once

#include <memory>
#include <string>
#include <vector>
#include <cmath> 

// Forward declare
struct IsothermModel;

// Parameter structs
struct HenryIsothermParameters 
{
    double k = 1.0;
    double Q = 1.0;
};

struct LangmuirIsothermParameters 
{
    double m1 = 1.0;
    double b1 = 1.0;
    double q1 = 1.0;
};

struct ExtendedLangmuirIsothermParameters
{
    double m1 = 1.0;
    double b1 = 1.0;
    double q1 = 1.0;
};

struct DualSiteLangmuirIsothermParameters 
{
    double m1 = 1.0;
    double b1 = 1.0;
    double q1 = 1.0;
    double m2 = 1.0;
    double b2 = 1.0;
    double q2 = 1.0;
};

// Namespace for isotherm functions
namespace IsothermFunction 
{
    void Henry(double& qi, double T, double P, double yi, const HenryIsothermParameters& ips);
    void Langmuir(double& qi, double T, double P, double yi, const LangmuirIsothermParameters& ips);
    void DualSiteLangmuir(double& qi, double T, double P, double yi, const DualSiteLangmuirIsothermParameters& ips);
    void ExtendedLangmuir(double& qi, double T, double P, double yi, const std::vector<double>& yis, const ExtendedLangmuirIsothermParameters& ip, const std::vector<ExtendedLangmuirIsothermParameters>& ips);
}

// Abstract isotherm interface
struct IIsotherm 
{
    virtual void update(double& qi, double T, double P, double yi, const std::vector<double>& yis) = 0;
    virtual ~IIsotherm() = default; // Virtual destructor
};

// Concrete isotherm classes
struct HenryIsotherm : public IIsotherm 
{
    HenryIsothermParameters isothermParams{};

    void update(double& qi, double T, double P, double yi, const std::vector<double>& yis) override;
};

struct LangmuirIsotherm : public IIsotherm 
{
    LangmuirIsothermParameters isothermParams{};

    void update(double& qi, double T, double P, double yi, const std::vector<double>& yis) override;
};

struct ExtendedLangmuirIsotherm : public IIsotherm
{
    ExtendedLangmuirIsothermParameters isothermParams{};
    std::vector<IsothermModel>& isothermModel;

    ExtendedLangmuirIsotherm(ExtendedLangmuirIsothermParameters params, std::vector<IsothermModel>& model)
        : isothermParams(params), isothermModel(model) {}

    void update(double& qi, double T, double P, double yi, const std::vector<double>& yis) override;
};

struct DualSiteLangmuirIsotherm : public IIsotherm 
{
    DualSiteLangmuirIsothermParameters isothermParams{};

    void update(double& qi, double T, double P, double yi, const std::vector<double>& yis) override;
};

// Factory functions
std::unique_ptr<IIsotherm> createIsotherm(const HenryIsothermParameters& params);
std::unique_ptr<IIsotherm> createIsotherm(const LangmuirIsothermParameters& params);
std::unique_ptr<IIsotherm> createIsotherm(const ExtendedLangmuirIsothermParameters& params);
std::unique_ptr<IIsotherm> createIsotherm(const DualSiteLangmuirIsothermParameters& params);

// IsothermModel struct
struct IsothermModel 
{
    explicit IsothermModel(std::string component) : component(std::move(component)) {}

    std::string component;
    std::unique_ptr<IIsotherm> isotherm = std::make_unique<LangmuirIsotherm>(); // Default isotherm model
    double ki = 1.0;  // Mass transfer coefficient [s]
    double Hads = 10000;  // Heat of adsorption [J/mol]
};
