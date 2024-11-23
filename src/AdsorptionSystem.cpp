#include "AdsorptionSimulator/AdsorptionSystem.h"

void AdsorptionSystem::addComponent(const std::string& component)
{
	fluid.addComponent(component);
	reactor.addComponent(component);
}

void AdsorptionSystem::removeComponent(const std::string& component)
{
	fluid.removeComponent(component);
	reactor.removeComponent(component);
}

PorousMedia& AdsorptionSystem::getAdsorbent(const std::string& adsorbent)
{
	return reactor.getLayer(adsorbent);
}
