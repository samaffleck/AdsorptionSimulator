#include "AdsorptionSimulator/Cycle.h"
#include "AdsorptionSimulator/Fluid.h"
#include "AdsorptionSimulator/Reactor.h"

#include <stdexcept>
#include <iostream>


void Cycle::run()
{
    reactor.initialise();

    for (int cycleNumber = 0; cycleNumber < numberOfCycles; ++cycleNumber) 
    {
        std::cout << "Running cycle " << cycleNumber + 1 << " of " << numberOfCycles << std::endl;
        runCycle();
    }

    std::cout << "Simulation Finished!";

    if (logData) fluidDataLogger.saveDataToCSV("C:/Users/samaf/OneDrive/Desktop/AdsResults");
}

void Cycle::runCycle()
{
    for (const auto& [stepName, step] : steps)
    {
        std::cout << "Executing step: " << stepName << std::endl;

        reactor.setBoundaryConditions(step);

        double elapsedTime = 0.0;
        while (elapsedTime < step.stepTime)
        {
            reactor.integrate(dt);
            if (logData) logSimulationData();
            elapsedTime += dt;
            std::cout << "Time: " << elapsedTime << " [s]\n";
        }
    }
}

void Cycle::logSimulationData()
{
    // Log the data
    fluidDataLogger.timeSeriesData.emplace_back(reactor.fluidData);
}

bool Cycle::performChecks() const
{
    if (!reactor.performChecks())
    {
        return false;
    }

    if (steps.empty()) 
    {
        std::cerr << "Cycle check failed: No steps defined." << std::endl;
        return false;
    }

    for (const auto& [stepName, step] : steps) 
    {
        if (step.stepTime <= 0.0) {
            std::cerr << "Cycle check failed: Step '" << stepName << "' has invalid step time." << std::endl;
            return false;
        }
    }

    std::cout << "Cycle checks passed successfully." << std::endl;
    return true;
}

void Cycle::addStep(std::string stepName)
{
    if (steps.find(stepName) != steps.end()) 
    {
        return;
    }

    steps.emplace(stepName, Step());
}

void Cycle::removeStep(std::string stepName)
{
    auto it = steps.find(stepName);
    if (it == steps.end()) 
    {
        return;
    }

    steps.erase(it);
}

Step& Cycle::getStep(std::string stepName)
{
    auto it = steps.find(stepName);
    if (it == steps.end()) 
    {
        throw std::runtime_error("Step with name '" + stepName + "' not found.");
    }

    return it->second;
}

void Cycle::updateCycleTime()
{
    cycleTime = 0.0;
    for (const auto& [stepName, step] : steps) 
    {
        cycleTime += step.stepTime;
    }

    std::cout << "Total cycle time updated to " << cycleTime << " seconds." << std::endl;
}
