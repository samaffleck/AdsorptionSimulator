#include "AdsorptionSimulator/Fluid.h"


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


bool Fluid::isComponentValid(const std::string& componentName)
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


bool Fluid::isComponentInList(const std::string& componentName)
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
