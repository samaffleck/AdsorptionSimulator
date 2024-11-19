#include "AdsorptionSimulator/AdsorptionSystem.h"

void AdsorptionSystem::addComponent(const std::string& component)
{
	fluid.addComponent(component);
	reactor.addIsothermModel(component);
	reactor.resizeData();
}

void AdsorptionSystem::removeComponent(const std::string& component)
{
	fluid.removeComponent(component);
	reactor.removeIsothermModel(component);
	reactor.resizeData();
}

PorousMedia& AdsorptionSystem::getAdsorbent(const std::string& adsorbent)
{
	return reactor.getLayer(adsorbent);
}
