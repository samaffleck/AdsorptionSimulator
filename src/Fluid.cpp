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

	// Create all of the files
	std::ofstream pressureFile;
	std::ofstream velocityFile;

	pressureFile.open(fileDirectory + "/pressure.csv");
	velocityFile.open(fileDirectory + "/velocity.csv");

	for (const auto& data : timeSeriesData)
	{
		logVariableToFile(data.P, pressureFile);
		logVariableToFile(data.u, velocityFile);
	}

	pressureFile.close();
	velocityFile.close();
}

void FluidDataLogger::logVariableToFile(const Eigen::VectorXd& var, std::ofstream& variableFile) const
{
	for (int n = 0; n < var.size(); ++n)
	{
		variableFile << var[n] << ",";
	}
	variableFile << std::endl;
}
