#pragma once

#include "BoundaryCondition.h"
#include "Reactor.h"

#include <unordered_map>


struct Reactor;

struct Step
{
	Inflow inflow;
	Outflow outflow;
	double stepTime = 100.0;	// [s]
};

struct Cycle
{
	Cycle(Reactor& reactor) : reactor(reactor) {}

	void run();
	bool performChecks() const;

	void addStep(const std::string& stepName);
	void removeStep(const std::string& stepName);
	Step& getStep(const std::string& stepName);
	
	double dt = 0.1;				// Time step [s]
	int numberOfCycles = 1;			// Number of cycles

private:
	void updateCycleTime();

private:
	Reactor& reactor;
	std::unordered_map<std::string, Step> steps;
	double cycleTime{};

};
