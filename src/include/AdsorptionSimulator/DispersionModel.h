#pragma once

struct FluidData;

struct IDispersionModel
{
	virtual ~IDispersionModel() = default;
	virtual void update(FluidData& fluidData) = 0;
};


struct DispersionModelConstant : public IDispersionModel
{
	DispersionModelConstant(double Dl) : constantDl(Dl) {}

	void update(FluidData& fluidData);

	double constantDl;
};
