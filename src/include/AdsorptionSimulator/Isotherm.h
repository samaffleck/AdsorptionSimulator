#pragma once

#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <cmath> 

// Forward declare
struct IsothermModel;
struct PorousMedia;
struct FluidData;

// Isotherm types
enum class IsothermType
{
    INERT, 
    HENRY,
    LANGMUIR,
    DUALSITELANGMUIR
};


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
    void Inert(double& qi);
    void Henry(double& qi, double T, double P, double yi, const HenryIsothermParameters& ips);
    void Langmuir(double& qi, double T, double P, double yi, const LangmuirIsothermParameters& ips);
    void DualSiteLangmuir(double& qi, double T, double P, double yi, const DualSiteLangmuirIsothermParameters& ips);
}

struct IIsotherm 
{
    virtual void update(double& qi, double T, double P, double yi) = 0;
    virtual ~IIsotherm() = default; 
};

struct InertIsotherm : public IIsotherm
{
    void update(double& qi, double T, double P, double yi) override;
};

struct HenryIsotherm : public IIsotherm 
{
    HenryIsotherm(HenryIsothermParameters ips) : isothermParams(ips) {}
    HenryIsothermParameters isothermParams;

    void update(double& qi, double T, double P, double yi) override;
};

struct LangmuirIsotherm : public IIsotherm 
{
    LangmuirIsotherm(LangmuirIsothermParameters ips) : isothermParams(ips) {}
    LangmuirIsothermParameters isothermParams{};

    void update(double& qi, double T, double P, double yi) override;
};

struct DualSiteLangmuirIsotherm : public IIsotherm 
{
    DualSiteLangmuirIsotherm(DualSiteLangmuirIsothermParameters ips) : isothermParams(ips) {}
    DualSiteLangmuirIsothermParameters isothermParams{};

    void update(double& qi, double T, double P, double yi) override;
};

// Factory functions
std::unique_ptr<IIsotherm> createIsotherm();
std::unique_ptr<IIsotherm> createIsotherm(const HenryIsothermParameters& params);
std::unique_ptr<IIsotherm> createIsotherm(const LangmuirIsothermParameters& params);
std::unique_ptr<IIsotherm> createIsotherm(const DualSiteLangmuirIsothermParameters& params);

// IsothermModel struct
struct IsothermModel 
{
public:
    void addComponent(std::string component);
    void removeComponent(const std::string& component);
    void updateIsotherm(FluidData& fluidData);
    void updateSourceTerms(FluidData& fluidData);
    IIsotherm* getIsotherm(const std::string& component);
    
    std::unique_ptr<IIsotherm> isotherm = std::make_unique<InertIsotherm>();    // Default isotherm model is inert
    
private:
    std::vector<std::string> components{};

};
