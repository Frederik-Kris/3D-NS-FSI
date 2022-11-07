/*
 * FlowVariableScalarStructs.h
 *
 *  Created on: Oct 10, 2022
 *      Author: frederk
 */

#ifndef SRC_FLOWVARIABLEGROUPSTRUCTS_H_
#define SRC_FLOWVARIABLEGROUPSTRUCTS_H_

#include "ConfigSettings.h"

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

struct PrimitiveVariablesScalars
{
public:
	double u;	// Velocity component in x-direction
	double v;	// Velocity component in y-direction
	double w;	// Velocity component in z-direction
	double p;	// Pressure
	double T;	// Temperature

	PrimitiveVariablesScalars(double u, double v, double w, double p, double T) :
		u{u}, v{v}, w{w}, p{p}, T{T}
	{}
};

struct TransportPropertiesScalars
{
public:
	double mu;		// Dynamic viscosity
	double kappa;	// Thermal conductivity

	TransportPropertiesScalars(double mu, double kappa) :
		mu{mu}, kappa{kappa}
	{}
};

struct ConservedVariablesArrayGroup
{
	Array3D_d rho;		// Mass density
	Array3D_d rho_u;	// Momentum density in x-direction
	Array3D_d rho_v;	// Momentum density in y-direction
	Array3D_d rho_w;	// Momentum density in z-direction
	Array3D_d rho_E;	// Total specific energy per volume

	ConservedVariablesArrayGroup(size_t gridSizeX, size_t gridSizeY, size_t gridSizeZ) :
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

	PrimitiveVariablesArrayGroup(size_t gridSizeX, size_t gridSizeY, size_t gridSizeZ) :
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

	TransportPropertiesArrayGroup(size_t gridSizeX, size_t gridSizeY, size_t gridSizeZ) :
	mu   (gridSizeX, gridSizeY, gridSizeZ),
	kappa(gridSizeX, gridSizeY, gridSizeZ)
	{}
};

struct RK4slopesArrayGroup
{
	ConservedVariablesArrayGroup k1, k2, k3, k4; // 4 slopes per conserved variable

	RK4slopesArrayGroup(size_t gridSizeX, size_t gridSizeY, size_t gridSizeZ) :
	k1(gridSizeX, gridSizeY, gridSizeZ),
	k2(gridSizeX, gridSizeY, gridSizeZ),
	k3(gridSizeX, gridSizeY, gridSizeZ),
	k4(gridSizeX, gridSizeY, gridSizeZ)
	{}
};

// Computes scalar values for conserved variables, based on the primitive variables.
// Intent: Can use this when applying initial conditions (IC) or boundary conditions (BC).
// E.g., decide on the primitives and use this function to get the values for the conserved variables.
inline ConservedVariablesScalars deriveConservedVariables(const PrimitiveVariablesScalars& primitiveVariables, const ConfigSettings& params)
{
	const double u{primitiveVariables.u}, v{primitiveVariables.v}, w{primitiveVariables.w}, p{primitiveVariables.p}, T{primitiveVariables.T};
	double rho   = ( params.Gamma * p - T ) / ( 1 + T );
	double rho_u = (1 + rho) * u;
	double rho_v = (1 + rho) * v;
	double rho_w = (1 + rho) * w;
	double rho_E = p / ( params.Gamma - 1 ) + (1 + rho)/2 * ( u*u + v*v + w*w );
	return ConservedVariablesScalars(rho, rho_u, rho_v, rho_w, rho_E);
}

// Computes scalar values for primitive variables, based on the conserved variables.
// Intent: Can use this when applying initial conditions (IC) or boundary conditions (BC).
// E.g., decide on the conserved and use this function to get the values for the primitive variables.
inline PrimitiveVariablesScalars derivePrimitiveVariables(const ConservedVariablesScalars& conservedVariables, const ConfigSettings& params)
{
	const double rho{conservedVariables.rho}, rho_u{conservedVariables.rho_u}, rho_v{conservedVariables.rho_v}, rho_w{conservedVariables.rho_w}, rho_E{conservedVariables.rho_E};
	double u = rho_u / rho;
	double v = rho_v / rho;
	double w = rho_w / rho;
	double p = ( params.Gamma - 1 )*( rho_E - (1 + rho)/2 * ( u*u + v*v + w*w ));
	double T = ( params.Gamma * p - rho ) / ( 1 + rho );
	return PrimitiveVariablesScalars(u, v, w, p, T);
}

// Computes scalar values for transport properties, based on the primitive variables.
// Intent: Can use this when applying initial conditions (IC) or boundary conditions (BC).
inline TransportPropertiesScalars deriveTransportProperties(const PrimitiveVariablesScalars& primitiveVariables, const ConfigSettings& params)
{
	const double T{primitiveVariables.T};
	double ScPlusOne = 1 + params.T_0 / params.sutherlands_C2;
	double mu = pow( 1+T, 1.5 ) * ScPlusOne / ( params.Re*( T + ScPlusOne ) );
	double kappa = mu / ( (params.Gamma - 1) * params.Pr );
	return TransportPropertiesScalars(mu, kappa);
}


#endif /* SRC_FLOWVARIABLEGROUPSTRUCTS_H_ */




