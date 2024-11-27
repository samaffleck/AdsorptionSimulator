#include "AdsorptionSimulator/DispersionModel.h"
#include "AdsorptionSimulator/PorousMedia.h"


void DispersionModelConstant::update(FluidData& fluidData)
{
	for (int n = 0; n < fluidData.Dl.size(); ++n)
	{
		fluidData.Dl[n] = constantDl;
	}
}
