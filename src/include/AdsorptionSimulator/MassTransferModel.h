#pragma once

struct PorousMedia;

struct IMassTransferModel
{
	virtual ~IMassTransferModel() = default;
	virtual void update(PorousMedia& adsorbent) = 0;
};


struct LinearDrivingForce : public IMassTransferModel
{
	LinearDrivingForce(double ki) : ki(ki) {}

	void update(PorousMedia& adsorbent);

	double ki;
};
