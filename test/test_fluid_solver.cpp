#include "AdsorptionSimulator/AdsorptionSystem.h"

#include <iostream>


int main() 
{    
    AdsorptionSystem system;

    system.addComponent("Oxygen");
    system.addComponent("Nitrogen");

    system.getReactor().addLayer("AA");
    system.getReactor().addLayer("13X");
    
    system.getAdsorbent("13X").setIsothermModel("Oxygen", LangmuirIsothermParameters({2, 3, 4}));
    system.getAdsorbent("13X").setIsothermModel("Nitrogen", HenryIsothermParameters({0, 1}));

    system.getAdsorbent("13X").setMassTransferCoefficient("Oxygen", 1.0);
    system.getAdsorbent("13X").setMassTransferCoefficient("Nitrogen", 1.0);
    
    system.getReactor().setDispersionModel(0);

    system.getAdsorbent("13X").setNumberOfCells(10);

    system.getCycle().addStep("Feed");
    system.getCycle().getStep("Feed").inflow.location = BoundaryConditionLocation::BOTTOM;
    system.getCycle().getStep("Feed").inflow.u_in = 0.05;
    system.getCycle().getStep("Feed").inflow.T_in = 298;
    system.getCycle().getStep("Feed").inflow.y_in = {0.1, 0.9};
    system.getCycle().getStep("Feed").outflow.location = BoundaryConditionLocation::TOP;
    system.getCycle().getStep("Feed").outflow.P_out = 101325.0;
    system.getCycle().getStep("Feed").stepTime = 100.0;

    system.getCycle().run();

    return 0;
}
