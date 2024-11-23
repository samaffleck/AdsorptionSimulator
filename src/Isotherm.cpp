#include "AdsorptionSimulator/Isotherm.h"

#include <algorithm>
#include <stdexcept>

void IsothermFunction::Inert(double& qi)
{
	qi = 0;
}

void IsothermFunction::Henry(double& qi, double T, double P, double yi, const HenryIsothermParameters& ips)
{
	qi = ips.k * P * yi * std::exp(ips.Q / (8.314 * T));
}

void IsothermFunction::Langmuir(double& qi, double T, double P, double yi, const LangmuirIsothermParameters& ips)
{
	qi = ips.m1 * ips.b1 * std::exp(ips.q1 / (8.314 * T) * P * yi) / (1 + ips.b1 * std::exp(ips.q1 / (8.314 * T) * P * yi));
}

void IsothermFunction::DualSiteLangmuir(double& qi, double T, double P, double yi, const DualSiteLangmuirIsothermParameters& ips)
{
	qi = (ips.m1 * ips.b1 * std::exp(ips.q1 / (8.314 * T) * P * yi) / (1 + ips.b1 * std::exp(ips.q1 / (8.314 * T) * P * yi))) + 
		(ips.m2 * ips.b2 * std::exp(ips.q2 / (8.314 * T) * P * yi) / (1 + ips.b2 * std::exp(ips.q2 / (8.314 * T) * P * yi)));
}

void InertIsotherm::update(double& qi, double T, double P, double yi)
{
	IsothermFunction::Inert(qi);
}

void HenryIsotherm::update(double& qi, double T, double P, double yi)
{
	IsothermFunction::Henry(qi, T, P, yi, isothermParams);
}

void LangmuirIsotherm::update(double& qi, double T, double P, double yi)
{
	IsothermFunction::Langmuir(qi, T, P, yi, isothermParams);
}

void DualSiteLangmuirIsotherm::update(double& qi, double T, double P, double yi)
{
	IsothermFunction::DualSiteLangmuir(qi, T, P, yi, isothermParams);
}

std::unique_ptr<IIsotherm> createIsotherm()
{
	return std::unique_ptr<InertIsotherm>();
}

std::unique_ptr<IIsotherm> createIsotherm(const HenryIsothermParameters& params)
{
	return std::make_unique<HenryIsotherm>(params);
}

std::unique_ptr<IIsotherm> createIsotherm(const LangmuirIsothermParameters& params)
{
	return std::make_unique<LangmuirIsotherm>(params);
}

std::unique_ptr<IIsotherm> createIsotherm(const DualSiteLangmuirIsothermParameters& params)
{
	return std::make_unique<DualSiteLangmuirIsotherm>(params);
}

void IsothermModel::addComponent(std::string component)
{
	if (std::find(components.begin(), components.end(), component) != components.end()) 
	{
		return;
	}

	components.push_back(std::move(component));
	ki[components.back()] = 1.0;   
	Hads[components.back()] = 1e5;
}

void IsothermModel::removeComponent(const std::string& component)
{
	auto it = std::find(components.begin(), components.end(), component);
	if (it == components.end()) 
	{
		return;
	}

	components.erase(it);
	ki.erase(component);
	Hads.erase(component);
}

void IsothermModel::updateIsotherm(FluidData& fluidData)
{
	for (const auto& component : components)
	{
		for (int n = 1; n < fluidData.C.size() - 1; ++n) // Loop through central cells only
		{
			isotherm->update(fluidData.qi_sat[component][n], fluidData.T[n], fluidData.P[n], fluidData.yi[component][n]);
		}
	}
}

IIsotherm* IsothermModel::getIsotherm(const std::string& component)
{
	for (int i = 0; i < components.size(); ++i)
	{
		if (components[i] == component)
		{
			return isotherm.get();
		}
	}
	return nullptr;
}

