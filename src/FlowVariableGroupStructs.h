/*
 * FlowVariableScalarStructs.h
 *
 *  Created on: Oct 10, 2022
 *      Author: frederk
 */

#ifndef SRC_FLOWVARIABLEGROUPSTRUCTS_H_
#define SRC_FLOWVARIABLEGROUPSTRUCTS_H_


struct ConservedVariablesScalars
{
public:
	double rho;		// Mass density
	double rho_u;	// Momentum density in x-direction
	double rho_v;	// Momentum density in y-direction
	double rho_w;	// Momentum density in z-direction
	double rho_E;		// Total specific energy per volume

	ConservedVariablesScalars(double rho, double rho_u, double rho_v, double rho_w, double rho_E) :
		rho{rho}, rho_u{rho_u}, rho_v{rho_v}, rho_w{rho_w}, rho_E{rho_E}
	{}
};

struct ConservedVariablesArrayGroup
{
	Array3D_d rho;		// Mass density
	Array3D_d rho_u;	// Momentum density in x-direction
	Array3D_d rho_v;	// Momentum density in y-direction
	Array3D_d rho_w;	// Momentum density in z-direction
	Array3D_d rho_E;	// Total specific energy per volume

	ConservedVariablesArrayGroup(uint gridSizeX, uint gridSizeY, uint gridSizeZ) :
	rho  (gridSizeX, gridSizeY, gridSizeZ),
	rho_u(gridSizeX, gridSizeY, gridSizeZ),
	rho_v(gridSizeX, gridSizeY, gridSizeZ),
	rho_w(gridSizeX, gridSizeY, gridSizeZ),
	rho_E(gridSizeX, gridSizeY, gridSizeZ)
	{}
};

struct PrimitiveVariablesArrayGroup
{
	Array3D_d u;	// Velocity component in x-direction
	Array3D_d v;	// Velocity component in y-direction
	Array3D_d w;	// Velocity component in z-direction
	Array3D_d p;	// Pressure
	Array3D_d T;	// Temperature

	PrimitiveVariablesArrayGroup(uint gridSizeX, uint gridSizeY, uint gridSizeZ) :
	u(gridSizeX, gridSizeY, gridSizeZ),
	v(gridSizeX, gridSizeY, gridSizeZ),
	w(gridSizeX, gridSizeY, gridSizeZ),
	p(gridSizeX, gridSizeY, gridSizeZ),
	T(gridSizeX, gridSizeY, gridSizeZ)
	{}
};

struct TransportPropertiesArrayGroup
{
	Array3D_d mu;		// Dynamic viscosity
	Array3D_d kappa;	// Thermal conductivity

	TransportPropertiesArrayGroup(uint gridSizeX, uint gridSizeY, uint gridSizeZ) :
	mu   (gridSizeX, gridSizeY, gridSizeZ),
	kappa(gridSizeX, gridSizeY, gridSizeZ)
	{}
};

struct RK4slopesArrayGroup
{
	ConservedVariablesArrayGroup k1, k2, k3, k4; // 4 slopes per conserved variable

	RK4slopesArrayGroup(uint gridSizeX, uint gridSizeY, uint gridSizeZ) :
	k1(gridSizeX, gridSizeY, gridSizeZ),
	k2(gridSizeX, gridSizeY, gridSizeZ),
	k3(gridSizeX, gridSizeY, gridSizeZ),
	k4(gridSizeX, gridSizeY, gridSizeZ)
	{}
};


#endif /* SRC_FLOWVARIABLEGROUPSTRUCTS_H_ */




