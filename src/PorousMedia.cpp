#include "AdsorptionSimulator/PorousMedia.h"
#include "AdsorptionSimulator/Isotherm.h"
#include "AdsorptionSimulator/ThomasAlgoritms.h"

#include <memory>


void PorousMedia::addIsothermModel(const std::string& component)
{
	isothermModels.emplace_back(component);
}


void PorousMedia::removeIsothermModel(const std::string& component)
{
	for (auto itt = isothermModels.begin(); itt != isothermModels.end(); itt++)
	{
		if (itt->component == component)
		{
			isothermModels.erase(itt);
			break;
		}
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


void PorousMedia::setIsothermModel(const std::string& component, const ExtendedLangmuirIsothermParameters& ips)
{
	IsothermModel* isothermModel = getIsothermModel(component);
	if (isothermModel != nullptr)
	{
		//isothermModel->isotherm = createIsotherm(ips, isothermModels);
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


IsothermModel* PorousMedia::getIsothermModel(const std::string& component)
{
	if (isothermModels.size() < 1) return nullptr;
	
	for (auto& isothermModel : isothermModels)
	{
		if (isothermModel.component == component)
		{
			return &isothermModel;
		}
	}

	return nullptr;
}


void PorousMedia::updateIsotherm()
{
	for (int i = 1; i < fluidData.C.size() - 1; ++i) // Loop through all central cells
	{
		std::vector<double> yis;
		for (int j = 0; j < isothermModels.size(); ++j) // Number of components
		{
			yis.emplace_back(fluidData.yi[j][i]);
		}

		for (int j = 0; j < isothermModels.size(); ++j) // Number of components
		{
			isothermModels[j].isotherm->update(fluidData.qi_sat[j][i], fluidData.T[i], fluidData.P[i], fluidData.yi[j][i], yis);
		}
	}
}


void PorousMedia::updateSourceTerms()
{
	for (int i = 1; i < fluidData.C.size() - 1; ++i) // Loop through all central cells
	{
		fluidData.Sm[i] = 0;
		fluidData.Se[i] = 0;
		for (int j = 0; j < isothermModels.size(); ++j) // Number of components
		{
			fluidData.Smi[j][i] = density * eb * isothermModels[j].ki * (fluidData.qi_sat[j][i] - fluidData.qi[j][i]);
			fluidData.Sei[j][i] = fluidData.Smi[j][i] * isothermModels[j].Hads;
			fluidData.Sm[i] += fluidData.Smi[j][i];
			fluidData.Se[i] += fluidData.Sei[j][i];
		}
	}
}

void FluidData::resize(int sizeOfVectors, int numberOfComponents)
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

	u.resize(sizeOfVectors - 1);			// Staggered grid is used for velocity, so 1 less cell is needed
	u.setConstant(0.);
	massFlow.resize(sizeOfVectors - 1);	// Staggered grid is used for velocity, so 1 less cell is needed
	massFlow.setConstant(0.);
	molarFlow.resize(sizeOfVectors - 1);	// Staggered grid is used for velocity, so 1 less cell is needed
	molarFlow.setConstant(0.);
	volumeFlow.resize(sizeOfVectors - 1);	// Staggered grid is used for velocity, so 1 less cell is needed
	volumeFlow.setConstant(0.);

	// Component parameters
	yi.resize(numberOfComponents);
	Ci.resize(numberOfComponents);
	qi.resize(numberOfComponents);
	qi_sat.resize(numberOfComponents);
	Smi.resize(numberOfComponents);
	Sei.resize(numberOfComponents);
	for (int i = 0; i < numberOfComponents; ++i)
	{
		yi[i].resize(sizeOfVectors);
		yi[i].setConstant(0.);
		Ci[i].resize(sizeOfVectors);
		Ci[i].setConstant(0.);
		qi[i].resize(sizeOfVectors);
		qi[i].setConstant(0.);
		qi_sat[i].resize(sizeOfVectors);
		qi_sat[i].setConstant(0.);
		Smi[i].resize(sizeOfVectors);
		Smi[i].setConstant(0.);
		Sei[i].resize(sizeOfVectors);
		Sei[i].setConstant(0.);
	}
}
