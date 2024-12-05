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
	reactor->updateLength();
}

double PorousMedia::getLayerLength() const
{
	return L;
}

void PorousMedia::setNumberOfCells(int numberOfCells_)
{
	numberOfCells = numberOfCells_;
	dx = L / numberOfCells;
	reactor->resizeData();
}

int PorousMedia::getNumberOfCells() const
{
	return numberOfCells;
}

double PorousMedia::getCellWidth() const
{
	return dx;
}

void PorousMedia::initialise(int startIndex, int endIndex)
{
	setStartEndIndex(startIndex, endIndex);
	permeability = eb * eb * eb * dp * dp / (180 * (1 - eb) * (1 - eb));
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
	
	for (int n = m_StartIndex + 1; n < m_EndIndex; ++n)
	{
		de = eb * 0.5 * (reactor->fluidData.Dl[n] + reactor->fluidData.Dl[n + 1]) / dx;
		dw = eb * 0.5 * (reactor->fluidData.Dl[n] + reactor->fluidData.Dl[n - 1]) / dx;

		Aw[n] = -dw - std::max(reactor->fluidData.u[n - 1], 0.0);
		Ae[n] = -de - std::max(-reactor->fluidData.u[n], 0.0);
		Ap[n] = alpha + std::max(reactor->fluidData.u[n], 0.0) + std::max(-reactor->fluidData.u[n - 1], 0.0) + de + dw;
		X[n] = alpha * reactor->fluidDataLastStep.Ci[component][n] + density * eb * reactor->fluidData.Smi[component][n] * dx;
	}
}

void PorousMedia::updatePressure(double dt, Eigen::VectorXd& Ae, Eigen::VectorXd& Aw, Eigen::VectorXd& Ap, Eigen::VectorXd& X)
{
    // Constants
	double D = (permeability * reactor->fluidDataLastStep.P[0]) / (reactor->fluidData.vis[0] * dx);	// Diffusivity
    double alpha = eb * dx / dt;        															// Coefficient for time discretization
	
    for (int i = m_StartIndex + 1; i < m_EndIndex; ++i) // Interior nodes
    {
		D = (permeability * reactor->fluidDataLastStep.P[i]) / (reactor->fluidData.vis[i] * dx);

        Aw(i) = -D;
        Ae(i) = -D;
        Ap(i) = 2 * D + alpha;
        X(i) = alpha * reactor->fluidDataLastStep.P(i) + dx * 8.314 * reactor->fluidData.T[i] * reactor->fluidData.Sm[i];
    }
}

void PorousMedia::updateVelocity()
{
	double R = 1;
	for (int n = m_StartIndex; n < m_EndIndex; ++n)
	{
		R = permeability / (reactor->fluidData.vis[n] * dx);
        reactor->fluidData.u[n] = R * (reactor->fluidData.P[n] - reactor->fluidData.P[n + 1]);
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
		for (const auto& component : fluid->components) // Sub-Component list
		{
			fluidData.Smi[component][n] = massTransferCoefficient[component] * (fluidData.qi_sat[component][n] - fluidData.qi[component][n]);
			fluidData.Sei[component][n] = fluidData.Smi[component][n] * heatOfAdsorption[component];
			fluidData.Sm[n] += fluidData.Smi[component][n];
			fluidData.Se[n] += fluidData.Sei[component][n];
		}
	}
}

void PorousLayers::remove(std::string name)
{
	auto it = std::remove_if(m_Layers.begin(), m_Layers.end(),
		[&name](const std::unique_ptr<PorousMedia>& layer) {
		return layer->name == name;
	});

	if (it != m_Layers.end()) {
		m_Layers.erase(it, m_Layers.end());
	}
}

PorousMedia* PorousLayers::getLayer(const std::string& layerName)
{
	auto it = std::find_if(m_Layers.begin(), m_Layers.end(),
		[&layerName](const std::unique_ptr<PorousMedia>& layer) {
		return layer->name == layerName;
	});

	if (it != m_Layers.end()) {
		return it->get();
	}
	throw std::invalid_argument("Layer with the specified name does not exist.");
}
