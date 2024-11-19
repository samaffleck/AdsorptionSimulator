#pragma once 

#include "Fluid.h"
#include "Reactor.h"
#include "Cycle.h"

struct AdsorptionSystem
{
public:
	AdsorptionSystem() : fluid(), reactor(fluid), cycle(reactor) {}

	void addComponent(const std::string& component);
	void removeComponent(const std::string& component);

	PorousMedia& getAdsorbent(const std::string& adsorbent);
	Reactor& getReactor() { return reactor; }
	Fluid& getFluid() { return fluid; }
	Cycle& getCycle() { return cycle; }

private:
	Fluid fluid;
	Reactor reactor;
	Cycle cycle;

};
