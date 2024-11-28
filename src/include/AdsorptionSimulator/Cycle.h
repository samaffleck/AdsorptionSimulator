#pragma once

#include "AdsorptionSimulator/BoundaryCondition.h"

#include <unordered_map>
#include <string>

struct Reactor;
struct Inflow;
struct Outflow;

struct Step
{
	Inflow inflow{};
	Outflow outflow{};
	double stepTime = 100.0;	// [s]
};

struct Cycle
{
	Cycle(Reactor& reactor) : reactor(reactor) {}

	void run();
	bool performChecks() const;

	void addStep(std::string stepName);
	void removeStep(std::string stepName);
	Step& getStep(std::string stepName);
	
	double dt = 0.1;				// Time step [s]
	int numberOfCycles = 1;			// Number of cycles
	bool logData = true;

private:
	void updateCycleTime();
	void logSimulationData();

private:
	Reactor& reactor;
	std::unordered_map<std::string, Step> steps;
	double cycleTime{};

};
