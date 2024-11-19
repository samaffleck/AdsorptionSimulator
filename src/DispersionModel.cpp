#include "AdsorptionSimulator/DispersionModel.h"
#include "AdsorptionSimulator/PorousMedia.h"


void DispersionModelConstant::update(PorousMedia& adsorbent)
{
	for (int n = 0; n < adsorbent.fluidData.Dl.size(); ++n)
	{
		adsorbent.fluidData.Dl[n] = constantDl;
	}
}
