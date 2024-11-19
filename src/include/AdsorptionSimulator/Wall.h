#pragma once

struct Wall
{
public:
	Wall();
	~Wall() = default;

	void setInternalDiameter(double di);
	double getInternalDiameter() const;

	void setExternalDiameter(double di);
	double getExternalDiameter() const;

	void updateThickness();
	void updateArea();
	double getArea() const { return area; }

	double Cp = 450;			// Specific heat capacity [J/(kg-K)]
	double k = 16.3;			// Thermal conductivity [W/(m-K)]

private:
	double d_i = 0.1;			// Internal diameter [m]
	double d_o = 0.12;			// External diameter [m]
	double density = 7800;		// Material density [kg/m3]
	double thickness{};			// Wall thickness [m]
	double area{};

};
