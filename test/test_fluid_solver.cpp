#include "AdsorptionSimulator/AdsorptionSystem.h"

#include <iostream>


int main() 
{    
    AdsorptionSystem system;

    system.addComponent("Oxygen");
    system.addComponent("Nitrogen");
    system.addComponent("Water");
    system.removeComponent("Water");

    system.getReactor().addLayer("AA");
    system.getReactor().addLayer("13X");
    
    system.getReactor().getLayer("AA").setLayerLength(0.25);
    system.getReactor().getLayer("AA").dp = 0.0035;
    system.getReactor().getLayer("13X").setLayerLength(0.45);
    system.getReactor().getLayer("13X").dp = 0.0021;

    system.getAdsorbent("13X").setIsothermModel("Oxygen", LangmuirIsothermParameters({2, 3, 4}));
    system.getAdsorbent("13X").setIsothermModel("Nitrogen", HenryIsothermParameters({0, 1}));

    system.getAdsorbent("13X").setMassTransferCoefficient("Oxygen", 1.0);
    system.getAdsorbent("13X").setMassTransferCoefficient("Nitrogen", 1.0);
    
    system.getReactor().setDispersionModel(0);

    system.getAdsorbent("13X").setNumberOfCells(10);

    system.getReactor().initialCondition.P0 = 101325;
    system.getReactor().initialCondition.T0 = 298;
    system.getReactor().initialCondition.yi0["Oxygen"] = 0.0;
    system.getReactor().initialCondition.yi0["Nitrogen"] = 1.0;

    system.getCycle().addStep("Feed");
    system.getCycle().getStep("Feed").inflow.location = BoundaryConditionLocation::BOTTOM;
    system.getCycle().getStep("Feed").inflow.u_in = 0.05;
    system.getCycle().getStep("Feed").inflow.T_in = 298;
    system.getCycle().getStep("Feed").inflow.y_in["Oxygen"] = 0.21;
    system.getCycle().getStep("Feed").inflow.y_in["Nitrogen"] = 0.79;
    system.getCycle().getStep("Feed").outflow.location = BoundaryConditionLocation::TOP;
    system.getCycle().getStep("Feed").outflow.P_out = 101325.0;
    system.getCycle().getStep("Feed").stepTime = 100.0;

    system.getCycle().run();

    return 0;
}
