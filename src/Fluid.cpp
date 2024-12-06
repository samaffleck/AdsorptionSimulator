#include "AdsorptionSimulator/Fluid.h"

#include <algorithm>  
#include <filesystem>
#include <fstream>


void Fluid::addComponent(const std::string& componentName)
{
	if (!isComponentValid(componentName)) return;
	if (isComponentInList(componentName)) return;
	components.emplace_back(componentName);
}

void Fluid::removeComponent(const std::string& componentName)
{
	if (!isComponentValid(componentName)) return;
	if (!isComponentInList(componentName)) return;
	components.erase(std::find(components.begin(), components.end(), componentName));
}

bool Fluid::isComponentValid(const std::string& componentName) const
{
	for (const auto& component : Components::componentNames)
	{
		if (component == componentName)
		{
			return true;
		}
	}

	return false;
}

bool Fluid::isComponentInList(const std::string& componentName) const
{
	for (const auto& component : components)
	{
		if (component == componentName)
		{
			return true;
		}
	}

	return false;
}

void FluidDataLogger::saveDataToCSV(const std::string& fileDirectory) const
{
	std::filesystem::create_directory(fileDirectory);

	auto NC = int(timeSeriesData[0].Ci.size()); // Number of componenets

	// Create all of the files
	std::ofstream pressureFile;
	std::ofstream velocityFile;
	std::vector<std::ofstream> moleFractionFiles(NC);

	pressureFile.open(fileDirectory + "/pressure.csv");
	velocityFile.open(fileDirectory + "/velocity.csv");

	int i = 0;
	for (const auto& [componentName, data] : timeSeriesData[0].Ci)
	{
		moleFractionFiles[i].open(fileDirectory + "/" + componentName + "_mole_fraction.csv");
		i++;
	}

	for (const auto& data : timeSeriesData)
	{
		logVariableToFile(data.P, pressureFile);
		logVariableToFile(data.u, velocityFile);
		
		i = 0;
		for (const auto& [componentName, componentData] : data.yi)
		{
			logVariableToFile(componentData, moleFractionFiles[i]);
			i++;
		}
	}

	pressureFile.close();
	velocityFile.close();

	for (i = 0; i < NC; ++i)
	{
		moleFractionFiles[i].close();
	}
}

void FluidDataLogger::logVariableToFile(const Eigen::VectorXd& var, std::ofstream& variableFile) const
{
	for (int n = 0; n < var.size(); ++n)
	{
		variableFile << var[n] << ",";
	}
	variableFile << std::endl;
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
