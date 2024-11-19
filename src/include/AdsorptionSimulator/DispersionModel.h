#pragma once

struct PorousMedia;

struct IDispersionModel
{
	virtual ~IDispersionModel() = default;
	virtual void update(PorousMedia& adsorbent) = 0;
};


struct DispersionModelConstant : public IDispersionModel
{
	DispersionModelConstant(double Dl) : constantDl(Dl) {}

	void update(PorousMedia& adsorbent);

	double constantDl;
};
