#pragma once

#include "Reactor.h"


namespace Solver
{
	void updateMoleFraction(PorousMedia& adsorbent, double dt);
	void updateVelocity(PorousMedia& adsorbent, double dt);
	void updateTemperature(PorousMedia& adsorbent, double dt);
	void updateMolarDensity(PorousMedia& adsorbent, double dt);
	void updateSolidConcentration(PorousMedia& adsorbent, double dt);
	void updatePressure(PorousMedia& adsorbent, double dt);
	void updateWallTemperature(PorousMedia& adsorbent, double dt);
} // End Solver namespace
