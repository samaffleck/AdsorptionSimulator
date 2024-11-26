#include "AdsorptionSimulator/PorousMedia.h"
#include "AdsorptionSimulator/Isotherm.h"
#include "AdsorptionSimulator/ThomasAlgoritms.h"
#include "AdsorptionSimulator/Fluid.h"
#include "AdsorptionSimulator/Reactor.h"

#include "CoolProp.h"

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
	fluidData.ki[component].setConstant(ki);
}

void PorousMedia::setHeatOfAdsorption(const std::string& component, double Hads)
{
	fluidData.Hads[component].setConstant(Hads);
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
	fluidData.resize(sizeOfVectors, fluid);
	fluidDataLastStep.resize(sizeOfVectors, fluid);
}

void PorousMedia::integrate(double dt)
{
	// Explicit equations
	updateConstants();
	updateFlowrates();
	updateIsotherms();
	updateSourceTerms();
	
	// Implicit equations
	updateMoleFraction(dt);
	updatePressure(dt);
	updateVelocity();

	fluidDataLastStep = fluidData;
}

void PorousMedia::updateMoleFraction(double dt)
{
	Eigen::VectorXd Ae(fluidData.C.size());
	Eigen::VectorXd Aw(fluidData.C.size());
	Eigen::VectorXd Ap(fluidData.C.size());
	Eigen::VectorXd X(fluidData.C.size());

	double alpha = et * dx / dt;

	for (const auto& component : fluid.components)
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

			Aw[n] = -dw - std::max(fluidData.u[n - 1], 0.0);
			Ae[n] = -de - std::max(-fluidData.u[n], 0.0);
			Ap[n] = alpha + std::max(fluidData.u[n], 0.0) + std::max(-fluidData.u[n - 1], 0.0) + de + dw;
			X[n] = alpha * fluidDataLastStep.Ci[component][n] + density * eb * fluidData.Smi[component][n] * dx;

			LinearSolver::TDMA(Ap, Ae, Aw, X, fluidData.Ci[component]);
		}
	}
}

void PorousMedia::updatePressure(double dt)
{
    int N_total = fluidData.P.size();       // Total number of nodes including ghost nodes

    Eigen::VectorXd Ae(N_total);
    Eigen::VectorXd Aw(N_total);
    Eigen::VectorXd Ap(N_total);
    Eigen::VectorXd X(N_total);

    // Constants
    double D = k / (fluidData.vis[0] / fluidData.rho[0]);   // Diffusivity
    double S = eb * compressibility;                 		// Storage coefficient
    double alpha = S * dx * dx / dt;        				// Coefficient for time discretization

    Ae.setZero();
    Aw.setZero();
    Ap.setZero();
    X.setZero();

	int i = 0;
	Aw(i) = 0;
	Ae(i) = 0;
	Ap(i) = 1;
	X(i) = fluidData.P(i);

    // Interior Nodes (i = 2 to N_physical - 1)
    for (i = 1; i < N_total - 1; ++i)
    {
		D = k / (fluidData.vis[i] / fluidData.rho[i]);

        Aw(i) = -D;
        Ae(i) = -D;
        Ap(i) = 2 * D + alpha;
        X(i) = alpha * fluidDataLastStep.P(i);
    }

	i = N_total - 1;
	D = k / (fluidData.vis[i] / fluidData.rho[i]);
	Aw(i) = 0;
	Ae(i) = 0;
	Ap(i) = 1;
	X(i) = fluidData.P(i);

    LinearSolver::TDMA(Ap, Ae, Aw, X, fluidData.P);
}

void PorousMedia::updateVelocity()
{
	for (int n = 0; n < fluidData.u.size(); ++n)
	{
        fluidData.u[n] = - (k / (fluidData.vis[0] / fluidData.rho[0])) * (fluidData.P[n + 1] - fluidData.P[n]) / dx;
	}
}

void PorousMedia::updateConstants()
{
	fluidData.Cp.setConstant(0.);
	fluidData.k.setConstant(0.);
	fluidData.rho.setConstant(0.);
	fluidData.vis.setConstant(0.);
	fluidData.Yt.setConstant(0.);

	for (int n = 0; n < fluidData.C.size(); ++n)
	{
		for (const auto& component : fluid.components)
		{
			fluidData.yi[component][n] = fluidData.Ci[component][n] / fluidData.C[n];
			fluidData.Cp[n] += fluidData.yi[component][n] * CoolProp::PropsSI("CPMOLAR", "P", fluidData.P[n], "T", fluidData.T[n], component);
			fluidData.k[n] += fluidData.yi[component][n] * CoolProp::PropsSI("CONDUCTIVITY", "P", fluidData.P[n], "T", fluidData.T[n], component);
			fluidData.rho[n] += fluidData.Ci[component][n] * CoolProp::Props1SI("molemass", component) * 1e-3;
			fluidData.vis[n] += fluidData.yi[component][n] * CoolProp::PropsSI("V", "P", fluidData.P[n], "T", fluidData.T[n], component);
			fluidData.Yt[n] += fluidData.yi[component][n];
		}
	}
}

void PorousMedia::updateIsotherms()
{
	for (auto& [type, model] : isothermModels)
	{
		model.updateIsotherm(fluidData);
	}
}

void PorousMedia::updateSourceTerms()
{
	for (auto& [type, model] : isothermModels)
	{
		model.updateSourceTerms(fluidData);
	}
}

void PorousMedia::updateFlowrates()
{
	double area = reactor.wall.getArea();

	for (int n = 0; n < fluidData.u.size(); ++n)
	{
		fluidData.volumeFlow[n] = fluidData.u[n] * area;
		fluidData.massFlow[n] = fluidData.volumeFlow[n] * fluidData.rho[n];
		fluidData.molarFlow[n] = fluidData.volumeFlow[n] * fluidData.C[n];
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
	ki.clear();
	Hads.clear();
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

		ki[component].resize(sizeOfVectors);
		ki[component].setConstant(0.);

		Hads[component].resize(sizeOfVectors);
		Hads[component].setConstant(0.);
	}
}
