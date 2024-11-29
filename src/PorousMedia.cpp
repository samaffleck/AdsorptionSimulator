#include "AdsorptionSimulator/PorousMedia.h"
#include "AdsorptionSimulator/Isotherm.h"
#include "AdsorptionSimulator/ThomasAlgoritms.h"
#include "AdsorptionSimulator/Fluid.h"
#include "AdsorptionSimulator/Reactor.h"

#include <memory>


void PorousMedia::addIsothermModel(const std::string& component)
{
	isothermModels[IsothermType::INERT].addComponent(component);
}

void PorousMedia::removeIsothermModel(const std::string& component)
{
	for (auto itt = isothermModels.begin(); itt != isothermModels.end(); itt++)
	{
		itt->second.removeComponent(component);
	}
}

void PorousMedia::setIsothermModel(const std::string& component, const HenryIsothermParameters& ips)
{
	// Remove the component from all existing isotherm models
	for (auto& [type, model] : isothermModels) 
	{
		model.removeComponent(component);
	}

	IsothermModel& henryModel = isothermModels[IsothermType::HENRY];
	henryModel.addComponent(component);

	henryModel.isotherm = std::make_unique<HenryIsotherm>(ips);
}

void PorousMedia::setIsothermModel(const std::string& component, const LangmuirIsothermParameters& ips)
{
	// Remove the component from all existing isotherm models
	for (auto& [type, model] : isothermModels) 
	{	
		model.removeComponent(component);
	}
	
	IsothermModel& langmuirModel = isothermModels[IsothermType::LANGMUIR];
	langmuirModel.addComponent(component);

	langmuirModel.isotherm = std::make_unique<LangmuirIsotherm>(ips);
}

void PorousMedia::setIsothermModel(const std::string& component, const DualSiteLangmuirIsothermParameters& ips)
{
	// Remove the component from all existing isotherm models
	for (auto& [type, model] : isothermModels) 
	{
		model.removeComponent(component);
	}

	IsothermModel& dualSiteLangmuirModel = isothermModels[IsothermType::DUALSITELANGMUIR];
	dualSiteLangmuirModel.addComponent(component);

	dualSiteLangmuirModel.isotherm = std::make_unique<DualSiteLangmuirIsotherm>(ips);
}

void PorousMedia::setIsothermModel(const std::string& component)
{
	// Remove the component from all existing isotherm models
	for (auto& [type, model] : isothermModels)
	{
		model.removeComponent(component);
	}

	IsothermModel& inertIsothermModel = isothermModels[IsothermType::INERT];
	inertIsothermModel.addComponent(component);

	inertIsothermModel.isotherm = std::make_unique<InertIsotherm>();
}

void PorousMedia::setMassTransferCoefficient(const std::string& component, double ki)
{
	massTransferCoefficient[component] = ki;
}

void PorousMedia::setHeatOfAdsorption(const std::string& component, double Hads)
{
	heatOfAdsorption[component] = Hads;
}

void PorousMedia::setLayerLength(double length) 
{
	L = length;
	reactor.updateLength();
}

double PorousMedia::getLayerLength() const
{
	return L;
}

void PorousMedia::setNumberOfCells(int numberOfCells_)
{
	numberOfCells = numberOfCells_;
	dx = L / numberOfCells;
	reactor.resizeData();
}

int PorousMedia::getNumberOfCells() const
{
	return numberOfCells;
}

double PorousMedia::getCellWidth() const
{
	return dx;
}

void PorousMedia::setStartEndIndex(int startIndex, int endIndex)
{
	m_StartIndex = startIndex;
	m_EndIndex = endIndex;
}

void PorousMedia::updateMoleFraction(const std::string &component, double dt, Eigen::VectorXd &Ae, Eigen::VectorXd &Aw, Eigen::VectorXd &Ap, Eigen::VectorXd &X)
{
	double alpha = et * dx / dt;
	double de = 0;
	double dw = 0;
	
	for (int n = m_StartIndex; n <= m_EndIndex; ++n)
	{
		de = eb * 0.5 * (reactor.fluidData.Dl[n] + reactor.fluidData.Dl[n + 1]) / dx;
		dw = eb * 0.5 * (reactor.fluidData.Dl[n] + reactor.fluidData.Dl[n - 1]) / dx;

		Aw[n] = -dw - std::max(reactor.fluidData.u[n - 1], 0.0);
		Ae[n] = -de - std::max(-reactor.fluidData.u[n], 0.0);
		Ap[n] = alpha + std::max(reactor.fluidData.u[n], 0.0) + std::max(-reactor.fluidData.u[n - 1], 0.0) + de + dw;
		X[n] = alpha * reactor.fluidDataLastStep.Ci[component][n] + density * eb * reactor.fluidData.Smi[component][n] * dx;
	}
}

void PorousMedia::updatePressure(double dt, Eigen::VectorXd& Ae, Eigen::VectorXd& Aw, Eigen::VectorXd& Ap, Eigen::VectorXd& X)
{
    // Constants
    double D = permeability / (reactor.fluidData.vis[0] / reactor.fluidData.rho[0]);   // Diffusivity
    double S = eb * compressibility;                 						// Storage coefficient
    double alpha = S * dx * dx / dt;        								// Coefficient for time discretization
	
    for (int i = m_StartIndex; i <= m_EndIndex; ++i)
    {
		D = k / (reactor.fluidData.vis[i] / reactor.fluidData.rho[i]);

        Aw(i) = -D;
        Ae(i) = -D;
        Ap(i) = 2 * D + alpha;
        X(i) = alpha * reactor.fluidDataLastStep.P(i);
    }
}

void PorousMedia::updateVelocity()
{
	for (int n = m_StartIndex; n <= m_EndIndex; ++n)
	{
        reactor.fluidData.u[n] = - (permeability / reactor.fluidData.vis[n]) * (reactor.fluidData.P[n + 1] - reactor.fluidData.P[n]) / dx;
	}
}

void PorousMedia::updateIsotherms(FluidData& fluidData)
{
	for (auto& [type, model] : isothermModels)
	{
		model.updateIsotherm(fluidData, m_StartIndex, m_EndIndex);
	}
}

void PorousMedia::updateSourceTerms(FluidData& fluidData)
{
	for (int n = m_StartIndex; n <= m_EndIndex; ++n) // Loop through all central cells
	{
		fluidData.Sm[n] = 0;
		fluidData.Se[n] = 0;
		for (const auto& component : fluid.components) // Sub-Component list
		{
			fluidData.Smi[component][n] = massTransferCoefficient[component] * (fluidData.qi_sat[component][n] - fluidData.qi[component][n]);
			fluidData.Sei[component][n] = fluidData.Smi[component][n] * heatOfAdsorption[component];
			fluidData.Sm[n] += fluidData.Smi[component][n];
			fluidData.Se[n] += fluidData.Sei[component][n];
		}
	}
}

void FluidData::resize(int sizeOfVectors, const Fluid& fluid)
{
	rho.resize(sizeOfVectors);
	rho.setConstant(0.);
	C.resize(sizeOfVectors);
	C.setConstant(0.);
	T.resize(sizeOfVectors);
	T.setConstant(0.);
	P.resize(sizeOfVectors);
	P.setConstant(0.);
	Sm.resize(sizeOfVectors);
	Sm.setConstant(0.);
	Se.resize(sizeOfVectors);
	Se.setConstant(0.);
	Yt.resize(sizeOfVectors);
	Yt.setConstant(0.);
	vis.resize(sizeOfVectors);
	vis.setConstant(0.);
	Cp.resize(sizeOfVectors);
	Cp.setConstant(0.);
	Cp_ads.resize(sizeOfVectors);
	Cp_ads.setConstant(0.);
	k.resize(sizeOfVectors);
	k.setConstant(0.);
	Dl.resize(sizeOfVectors);
	Dl.setConstant(0.);
	Lam.resize(sizeOfVectors);
	Lam.setConstant(0.);
	h_gw.resize(sizeOfVectors);
	h_gw.setConstant(0.);
	h_wa.resize(sizeOfVectors);
	h_wa.setConstant(0.);
	Re.resize(sizeOfVectors);
	Re.setConstant(0.);
	Sc.resize(sizeOfVectors);
	Sc.setConstant(0.);
	Pr.resize(sizeOfVectors);
	Pr.setConstant(0.);

	// Staggered grid is used for velocity, so 1 less cell is needed
	u.resize(sizeOfVectors - 1);
	u.setConstant(0.);
	massFlow.resize(sizeOfVectors - 1);
	massFlow.setConstant(0.);
	molarFlow.resize(sizeOfVectors - 1);
	molarFlow.setConstant(0.);
	volumeFlow.resize(sizeOfVectors - 1);
	volumeFlow.setConstant(0.);

	// Component parameters
	yi.clear();
	Ci.clear();
	qi.clear();
	qi_sat.clear();
	Smi.clear();
	Sei.clear();
	for (const auto& component : fluid.components)
	{
		yi[component].resize(sizeOfVectors);
		yi[component].setConstant(0.);

		Ci[component].resize(sizeOfVectors);
		Ci[component].setConstant(0.);

		qi[component].resize(sizeOfVectors);
		qi[component].setConstant(0.);

		qi_sat[component].resize(sizeOfVectors);
		qi_sat[component].setConstant(0.);

		Smi[component].resize(sizeOfVectors);
		Smi[component].setConstant(0.);

		Sei[component].resize(sizeOfVectors);
		Sei[component].setConstant(0.);
	}
}
