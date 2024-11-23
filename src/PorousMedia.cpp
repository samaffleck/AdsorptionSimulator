#include "AdsorptionSimulator/PorousMedia.h"
#include "AdsorptionSimulator/Isotherm.h"
#include "AdsorptionSimulator/ThomasAlgoritms.h"

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
	IsothermModel* isothermModel = getIsothermModel(component);
	if (isothermModel != nullptr)
	{
		isothermModel->isotherm = createIsotherm(ips);
	}
}


void PorousMedia::setIsothermModel(const std::string& component, const LangmuirIsothermParameters& ips)
{
	IsothermModel* isothermModel = getIsothermModel(component);
	if (isothermModel != nullptr)
	{
		isothermModel->isotherm = createIsotherm(ips);
	}
}

void PorousMedia::setMassTransferCoefficient(const std::string& component, double ki)
{
	IsothermModel* isothermModel = getIsothermModel(component);
	if (isothermModel != nullptr)
	{
		isothermModel->ki = ki;
	}
}

void PorousMedia::setHeatOfAdsorption(const std::string& component, double Hads)
{
	IsothermModel* isothermModel = getIsothermModel(component);
	if (isothermModel != nullptr)
	{
		isothermModel->Hads = Hads;
	}
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
	resizeData();
}

int PorousMedia::getNumberOfCells() const
{
	return numberOfCells;
}

void PorousMedia::resizeData()
{
	int sizeOfVectors = numberOfCells + 2; // 2 ghost cells are used either size of the domain
	fluidData.resize(sizeOfVectors, fluid.getNumberOfComponents());
	fluidDataLastStep.resize(sizeOfVectors, fluid.getNumberOfComponents());
}

void PorousMedia::updateConstants()
{
	updateIsotherm();
	updateSourceTerms();
}

void PorousMedia::updateMoleFraction(double dt)
{
	Eigen::VectorXd Ae(fluidData.C.size());
	Eigen::VectorXd Aw(fluidData.C.size());
	Eigen::VectorXd Ap(fluidData.C.size());
	Eigen::VectorXd X(fluidData.C.size());

	double alpha = et * dx / dt;

	for (int i = 0; i < fluidData.Ci.size(); ++i)
	{
		Ae.setConstant(0.);
		Aw.setConstant(0.);
		Ap.setConstant(0.);
		X.setConstant(0.);
		
		double de = 0;
		double dw = 0;
		
		for (int n = 1; n < fluidData.C.size() - 1; ++n)
		{
			de = eb * 0.5 * (fluidData.Dl[n] + fluidData.Dl[n + 1]) / dx;
			dw = eb * 0.5 * (fluidData.Dl[n] + fluidData.Dl[n - 1]) / dx;

			Aw[n] = -dw - std::max(fluidData.u[n], 0.0);
			Ae[n] = -de - std::max(-fluidData.u[n + 1], 0.0);
			Ap[n] = alpha + std::max(fluidData.u[n + 1], 0.0) + std::max(-fluidData.u[n], 0.0) + de + dw;
			X[n] = alpha * fluidDataLastStep.Ci[i][n] + fluidData.Smi[i][n] * dx;

			LinearSolver::TDMA(Ap, Ae, Aw, X, fluidData.Ci[i]);
		}
	}
}


void PorousMedia::setIsothermModel(const std::string& component, const DualSiteLangmuirIsothermParameters& ips)
{
	IsothermModel* isothermModel = getIsothermModel(component);
	if (isothermModel != nullptr)
	{
		isothermModel->isotherm = createIsotherm(ips);
	}
}

void PorousMedia::setIsothermModel(const std::string& component)
{
	IIsotherm* isotherm= getIsothermModel(component);
	if (isotherm != nullptr)
	{
		isotherm = createIsotherm();
	}
}

IIsotherm* PorousMedia::getIsothermModel(const std::string& component)
{
	if (isothermModels.size() < 1) return nullptr;
	
	IIsotherm* model = nullptr;

	for (auto& [type, isothermModel] : isothermModels)
	{
		model = isothermModel.getIsotherm(component);
	}

	return model;
}

void PorousMedia::updateIsotherms()
{
	for (int i = 0; i < fluidData.qi_sat.size(); ++i)
	{
		for (auto& [type, model] : isothermModels)
		{
			model.updateIsotherm(fluidData);
		}
	}
}

void PorousMedia::updateSourceTerms()
{
	for (int n = 1; n < fluidData.C.size() - 1; ++n) // Loop through all central cells
	{
		fluidData.Sm[n] = 0;
		fluidData.Se[n] = 0;
		for (int i = 0; i < isothermModels.size(); ++i) // Number of components
		{
			fluidData.Smi[i][n] = density * eb * isothermModels[i].ki * (fluidData.qi_sat[i][n] - fluidData.qi[i][n]);
			fluidData.Sei[i][n] = fluidData.Smi[i][n] * isothermModels[i].Hads;
			fluidData.Sm[n] += fluidData.Smi[i][n];
			fluidData.Se[n] += fluidData.Sei[i][n];
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
	for (auto& component : componentNames)
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
