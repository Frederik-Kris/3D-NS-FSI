/*
 * Simulation.cpp
 *
 *  Created on: Oct 6, 2022
 *      Author: frederk
 */

#include "Simulation.h"

Simulation::Simulation() :
params(), t{0}, timeLevel{0}
{

}

// Set the time step size as large as possible, within the stability criterion.
// Computes two time step sizes, using the inviscid CFL condition, and the viscous von Neumann condition
// and assigns the smallest one to dt (strictest criterion).
void Solver::updateTimeStepSize()
{
	// First, find the inviscid time step limit (CFL condition):
	uint nNodes = params.NI * params.NJ * params.NK;
	double maxSpectralRadiusX{0}, maxSpectralRadiusY{0}, maxSpectralRadiusZ{0};
	for (uint i{0}; i<nNodes; ++i)
	{
		double c_i = sqrt( (1 + params.Gamma * p(i)) / (1 + rho(i)) );    //<- Speed of sound at node i
		maxSpectralRadiusX = max( maxSpectralRadiusX, c_i + fabs(u(i)) );
		maxSpectralRadiusY = max( maxSpectralRadiusY, c_i + fabs(v(i)) );
		maxSpectralRadiusZ = max( maxSpectralRadiusZ, c_i + fabs(w(i)) );
	}
	double dt_conv = fabs(params.convStabilityLimit) / ( maxSpectralRadiusX / dx
			                                           + maxSpectralRadiusY / dy
												       + maxSpectralRadiusZ / dz );

	// Then find the viscous time step limit (von Neumann condition):
	double viscosityModifier = max( 4./3, params.Gamma / params.Pr );
	double maxViscosityFactor{0};   //<- Modified viscosity 'nu' used in the stability criterion
	for (uint i{0}; i<nNodes; ++i)
	{
		double nu = mu(i) / (rho(i)+1) * viscosityModifier;
		maxViscosityFactor = max( maxViscosityFactor, nu );
	}
	double dt_visc = fabs(params.viscStabilityLimit) / ( maxViscosityFactor * (1/(dx*dx) + 1/(dy*dy) + 1/(dz*dz)) );

	// Choose the strictest criterion. If that will take us past the end-time, then adapt dt to hit t_end exactly:
	dt = min( dt_conv, dt_visc );
	if ( params.stopCriterion == StopCriterionEnum::end_time && t + dt > params.t_end )
		dt = params.t_end - t;
}

// Check if computing another timestep will violate the stopping criterion. Returns true if simulation should stop.
// If the the criterion is to stop after a given number of timesteps, no more timesteps are required after this.
// To stop at an exact time, one last timestep with an adapted dt is required to hit the end-time.
bool Solver::checkStoppingCriterion()
{
	bool stopSimulating = true;		// until the opposite is proven.
	switch (params.stopCriterion)
	{
	case StopCriterionEnum::timesteps:
		if ( timeLevel < params.stopTimeLevel )
			stopSimulating = false;
		else
			stopSimulating = true;
		break;
	case StopCriterionEnum::end_time:
		if (t < params.t_end)
			stopSimulating = false;
		else
			stopSimulating = true;
		break;
	}
	return stopSimulating;
}


