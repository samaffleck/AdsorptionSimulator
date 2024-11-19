#include "AdsorptionSimulator/Isotherm.h"

void IsothermFunction::Henry(double& qi, double T, double P, double yi, const HenryIsothermParameters& ips)
{
}

void IsothermFunction::Langmuir(double& qi, double T, double P, double yi, const LangmuirIsothermParameters& ips)
{
}

void IsothermFunction::DualSiteLangmuir(double& qi, double T, double P, double yi, const DualSiteLangmuirIsothermParameters& ips)
{
}

void IsothermFunction::ExtendedLangmuir(double& qi, double T, double P, double yi, const std::vector<double>& yis, const ExtendedLangmuirIsothermParameters& ip, const std::vector<ExtendedLangmuirIsothermParameters>& ips)
{
}

void HenryIsotherm::update(double& qi, double T, double P, double yi, const std::vector<double>& yis)
{
}

void LangmuirIsotherm::update(double& qi, double T, double P, double yi, const std::vector<double>& yis)
{
}

void ExtendedLangmuirIsotherm::update(double& qi, double T, double P, double yi, const std::vector<double>& yis)
{
}

void DualSiteLangmuirIsotherm::update(double& qi, double T, double P, double yi, const std::vector<double>& yis)
{
}

std::unique_ptr<IIsotherm> createIsotherm(const HenryIsothermParameters& params)
{
	return std::unique_ptr<IIsotherm>();
}

std::unique_ptr<IIsotherm> createIsotherm(const LangmuirIsothermParameters& params)
{
	return std::unique_ptr<IIsotherm>();
}

std::unique_ptr<IIsotherm> createIsotherm(const ExtendedLangmuirIsothermParameters& params)
{
	return std::unique_ptr<IIsotherm>();
}

std::unique_ptr<IIsotherm> createIsotherm(const DualSiteLangmuirIsothermParameters& params)
{
	return std::unique_ptr<IIsotherm>();
}
