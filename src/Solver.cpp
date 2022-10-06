/*
 * Solver.cpp
 *
 *  Created on: Apr 28, 2021
 *      Author: frederik
 */


#include "Solver.h"

// Default constructor.
// Initializes the ConfigSettings 'params', which loads settings and parameters from 'ConfigFile'.
// Initializes all 3D arrays (but not their elements) containing flow variables and transport properties.
// Starts the simulation.
Solver::Solver() :
params(),
rho     (params.NI, params.NJ, params.NK),
rho_u   (params.NI, params.NJ, params.NK),
rho_v   (params.NI, params.NJ, params.NK),
rho_w   (params.NI, params.NJ, params.NK),
E       (params.NI, params.NJ, params.NK),
u       (params.NI, params.NJ, params.NK),
v       (params.NI, params.NJ, params.NK),
w       (params.NI, params.NJ, params.NK),
p       (params.NI, params.NJ, params.NK),
T       (params.NI, params.NJ, params.NK),
mu      (params.NI, params.NJ, params.NK),
kappa   (params.NI, params.NJ, params.NK),
k1_rho  (params.NI, params.NJ, params.NK),
k2_rho  (params.NI, params.NJ, params.NK),
k3_rho  (params.NI, params.NJ, params.NK),
k4_rho  (params.NI, params.NJ, params.NK),
k1_rho_u(params.NI, params.NJ, params.NK),
k2_rho_u(params.NI, params.NJ, params.NK),
k3_rho_u(params.NI, params.NJ, params.NK),
k4_rho_u(params.NI, params.NJ, params.NK),
k1_rho_v(params.NI, params.NJ, params.NK),
k2_rho_v(params.NI, params.NJ, params.NK),
k3_rho_v(params.NI, params.NJ, params.NK),
k4_rho_v(params.NI, params.NJ, params.NK),
k1_rho_w(params.NI, params.NJ, params.NK),
k2_rho_w(params.NI, params.NJ, params.NK),
k3_rho_w(params.NI, params.NJ, params.NK),
k4_rho_w(params.NI, params.NJ, params.NK),
k1_E    (params.NI, params.NJ, params.NK),
k2_E    (params.NI, params.NJ, params.NK),
k3_E    (params.NI, params.NJ, params.NK),
k4_E    (params.NI, params.NJ, params.NK),
interm_rho  (params.NI, params.NJ, params.NK),
interm_rho_u(params.NI, params.NJ, params.NK),
interm_rho_v(params.NI, params.NJ, params.NK),
interm_rho_w(params.NI, params.NJ, params.NK),
interm_E    (params.NI, params.NJ, params.NK),
savedSolutions{ 0 },
timeLevel{ 0 }, t{ 0 }
{
	setGridSpacings();
	applyStagnation_IC();
	if ( params.save_IC )
		storeCurrentSolution_csv();
	updateTimeStepSize();
	Clock statusReportTimer;
	while ( !checkStoppingCriterion() )
	{
		marchTimeStep();
		processOutput(statusReportTimer);
	}
	if ( params.save_final )
		storeCurrentSolution_csv();
	writeStatusReport_toScreen();
	writeOutputTimes();
	writeNormHistoryFiles();
}

// Calculate the space between nodes in the grid, based on domain size and no. of nodes.
void Solver::setGridSpacings()
{
	dx = params.L_x / (params.NI - 1);
	dy = params.L_y / (params.NJ - 1);
	dz = params.L_z / (params.NK - 1);
	cout << "Grid spacings set: dx = " << dx << " , dy = " << dy << " , dz = " << dz << endl;
}

// Applies stagnation initial condition (IC) by setting the values of flow variables at every node
void Solver::applyStagnation_IC()
{
	cout << "Applying stagnation initial condition... ";
	double u_IC{0}, v_IC{0}, w_IC{0}, p_IC{0}, T_IC{0};                  // <- Decide primitive variables
	double rho_IC, rho_u_IC, rho_v_IC, rho_w_IC, E_IC, mu_IC, kappa_IC;  // <- Other variables follow equations.
	getDerivedVariables_atPoint(u_IC, v_IC, w_IC, p_IC, T_IC, rho_IC, rho_u_IC, rho_v_IC, rho_w_IC, E_IC, mu_IC, kappa_IC);
	u    .setAll(u_IC);
	v    .setAll(v_IC);
	w    .setAll(w_IC);
	p    .setAll(p_IC);
	T    .setAll(T_IC);
	rho  .setAll(rho_IC);
	rho_u.setAll(rho_u_IC);
	rho_v.setAll(rho_v_IC);
	rho_w.setAll(rho_w_IC);
	E	 .setAll(E_IC);
	mu   .setAll(mu_IC);
	kappa.setAll(kappa_IC);
	cout << "Done" << endl << endl;
}

// Applies uniform flow initial condition (IC) by setting the values of flow variables at every node
void Solver::applyUniformFlow_IC()
{
	double velocity = params.M_0 / sqrt(3);
	double u_IC{velocity}, v_IC{velocity}, w_IC{velocity}, p_IC{0}, T_IC{0}; // <- Decide primitive variables
	double rho_IC, rho_u_IC, rho_v_IC, rho_w_IC, E_IC, mu_IC, kappa_IC;      // <- Other variables follow equations.
	getDerivedVariables_atPoint(u_IC, v_IC, w_IC, p_IC, T_IC, rho_IC, rho_u_IC, rho_v_IC, rho_w_IC, E_IC, mu_IC, kappa_IC);
	u    .setAll(u_IC);
	v    .setAll(v_IC);
	w    .setAll(w_IC);
	p    .setAll(p_IC);
	T    .setAll(T_IC);
	rho  .setAll(rho_IC);
	rho_u.setAll(rho_u_IC);
	rho_v.setAll(rho_v_IC);
	rho_w.setAll(rho_w_IC);
	E	 .setAll(E_IC);
	mu   .setAll(mu_IC);
	kappa.setAll(kappa_IC);
}

// Computes scalar values for conserved variables and transport properties, based on the primitive variables.
// Intent: Use this when applying initial conditions (IC) or boundary conditions (BC). Decide on the primitives and
// use this function to get the values for the other flow variables. Subscript _P denotes scalar values in one point.
void Solver::getDerivedVariables_atPoint(double u_P,double v_P, double w_P, double p_P, double T_P,
		double& rho_P, double& rho_u_P, double& rho_v_P, double& rho_w_P, double& E_P, double& mu_P, double& kappa_P)
{
	rho_P   = ( params.Gamma * p_P - T_P ) / ( 1 + T_P );
	rho_u_P = (1 + rho_P) * u_P;
	rho_v_P = (1 + rho_P) * v_P;
	rho_w_P = (1 + rho_P) * w_P;
	E_P     = p_P / ( params.Gamma - 1 ) + (1 + rho_P)/2 * ( u_P*u_P + v_P*v_P + w_P*w_P );
	double ScPlusOne = 1 + params.T_0 / params.sutherlands_C2;
	mu_P    = pow( 1+T_P, 1.5 ) * ScPlusOne / ( params.Re*( T_P + ScPlusOne ) );
	kappa_P = mu_P / ( (params.Gamma - 1) * params.Pr );
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

// Advances the conserved variables, primitive variables and transport properties from time t to t + dt, using RK4.
// Updates time step size per the stability criterion after advancing the solution.
// TODO: This function is long. Consider splitting it up into smaller ones.
void Solver::marchTimeStep()
{
	applyFilter_ifAppropriate(rho,   interm_rho  );

	compute_RK4_step_continuity(rho_u, rho_v, rho_w, k1_rho  );  // Compute step 1 (k1), i.e. the slopes at
	compute_RK4_step_xMomentum (rho_u,               k1_rho_u);  // time t, using Euler's method.
	compute_RK4_step_yMomentum (rho_v,               k1_rho_v);
	compute_RK4_step_zMomentum (rho_w,               k1_rho_w);
	compute_RK4_step_energy    (E,                   k1_E    );

	computeIntermediateSolution(rho,   k1_rho,   interm_rho,   dt/2);   // Compute intermediate solutions at
	computeIntermediateSolution(rho_u, k1_rho_u, interm_rho_u, dt/2);   // time t + dt/2, using the slopes k1.
	computeIntermediateSolution(rho_v, k1_rho_v, interm_rho_v, dt/2);
	computeIntermediateSolution(rho_w, k1_rho_w, interm_rho_w, dt/2);
	computeIntermediateSolution(E,     k1_E,     interm_E,     dt/2);

	updatePrimitiveVariables(interm_rho, interm_rho_u, interm_rho_v, interm_rho_w, interm_E);
	applyInjectionBC_risingInletVelocityChannelFlow(interm_rho, interm_rho_u, interm_rho_v, interm_rho_w, interm_E, t+dt/2);

	compute_RK4_step_continuity(interm_rho_u, interm_rho_v, interm_rho_w, k2_rho  );  // k2: The slopes at time t + dt/2
	compute_RK4_step_xMomentum (interm_rho_u,                             k2_rho_u);
	compute_RK4_step_yMomentum (interm_rho_v,                             k2_rho_v);
	compute_RK4_step_zMomentum (interm_rho_w,                             k2_rho_w);
	compute_RK4_step_energy    (interm_E,                                 k2_E    );

	computeIntermediateSolution(rho,   k2_rho,   interm_rho,   dt/2);   // Compute intermediate solutions at
	computeIntermediateSolution(rho_u, k2_rho_u, interm_rho_u, dt/2);   // time t + dt/2, using the slopes k2.
	computeIntermediateSolution(rho_v, k2_rho_v, interm_rho_v, dt/2);
	computeIntermediateSolution(rho_w, k2_rho_w, interm_rho_w, dt/2);
	computeIntermediateSolution(E,     k2_E,     interm_E,     dt/2);

	updatePrimitiveVariables(interm_rho, interm_rho_u, interm_rho_v, interm_rho_w, interm_E);
	applyInjectionBC_risingInletVelocityChannelFlow(interm_rho, interm_rho_u, interm_rho_v, interm_rho_w, interm_E, t+dt/2);

	compute_RK4_step_continuity(interm_rho_u, interm_rho_v, interm_rho_w, k3_rho  );  // k3: Again, the slopes at time
	compute_RK4_step_xMomentum (interm_rho_u,                             k3_rho_u);  //     t + dt/2, but now using k2
	compute_RK4_step_yMomentum (interm_rho_v,                             k3_rho_v);
	compute_RK4_step_zMomentum (interm_rho_w,                             k3_rho_w);
	compute_RK4_step_energy    (interm_E,                                 k3_E    );

	computeIntermediateSolution(rho,   k3_rho,   interm_rho,   dt);   // Compute intermediate solutions at
	computeIntermediateSolution(rho_u, k3_rho_u, interm_rho_u, dt);   // time t + dt, using the slopes k3.
	computeIntermediateSolution(rho_v, k3_rho_v, interm_rho_v, dt);
	computeIntermediateSolution(rho_w, k3_rho_w, interm_rho_w, dt);
	computeIntermediateSolution(E,     k3_E,     interm_E,     dt);

	updatePrimitiveVariables(interm_rho, interm_rho_u, interm_rho_v, interm_rho_w, interm_E);
	applyInjectionBC_risingInletVelocityChannelFlow(interm_rho, interm_rho_u, interm_rho_v, interm_rho_w, interm_E, t+dt);

	compute_RK4_step_continuity(interm_rho_u, interm_rho_v, interm_rho_w, k4_rho  );  // k4: The slopes at time t + dt
	compute_RK4_step_xMomentum (interm_rho_u,                             k4_rho_u);
	compute_RK4_step_yMomentum (interm_rho_v,                             k4_rho_v);
	compute_RK4_step_zMomentum (interm_rho_w,                             k4_rho_w);
	compute_RK4_step_energy    (interm_E,                                 k4_E    );

	compute_RK4_final_step(k1_rho,   k2_rho,   k3_rho,   k4_rho,   rho  , interm_rho  );  // Use a weighted average of the four slopes
	compute_RK4_final_step(k1_rho_u, k2_rho_u, k3_rho_u, k4_rho_u, rho_u, interm_rho_u);  // to advance the solution from time t, to
	compute_RK4_final_step(k1_rho_v, k2_rho_v, k3_rho_v, k4_rho_v, rho_v, interm_rho_v);  // t + dt. The solutions at the new time level
	compute_RK4_final_step(k1_rho_w, k2_rho_w, k3_rho_w, k4_rho_w, rho_w, interm_rho_w);  // are stored in the 'interm' arrays, so the
	compute_RK4_final_step(k1_E,     k2_E,     k3_E,     k4_E,     E    , interm_E    );  // norm of the change can be computed.

	updatePrimitiveVariables(interm_rho, interm_rho_u, interm_rho_v, interm_rho_w, interm_E);
	applyInjectionBC_risingInletVelocityChannelFlow(interm_rho, interm_rho_u, interm_rho_v, interm_rho_w, interm_E, t+dt);

	// At this stage, the conserved variables at the next time level are stored in the intermediate arrays.
	computeNorms_conservedVariables();	// Compute norm of the change between old and new time levels, and add to history.
	swapConservedVariables();			// Swap conserved variables to their main arrays by move-semantics.

	t += dt;
	++timeLevel;
	updateTimeStepSize();
}

// Filters one variable field, i.e., one solution array, 'filterVariable' and stores the filtered
// result in 'variableTemporaryStorage'. Then, the arrays are swapped by move-semantics.
// Only filters if the modulo of time level plus one, by the filter interval is zero.
// This causes the first filtering to happen as late as possible.
void Solver::applyFilter_ifAppropriate(Array3D_d& filterVariable, Array3D_d& variableTemporaryStorage)
{
	if(params.filterInterval > 0)
		if( (timeLevel+1) % params.filterInterval == 0 )
		{
			uint iMax{params.NI-1}, jMax{params.NJ-1}, kMax{params.NK-1};
			// Copy boundary nodes:
			for(uint i{0}; i<=iMax; ++i)
				for(uint j{0}; j<=jMax; ++j)
				{
					variableTemporaryStorage(i,j,0   ) = filterVariable(i,j,0   );
					variableTemporaryStorage(i,j,kMax) = filterVariable(i,j,kMax);
				}
			for(uint i{0}; i<=iMax; ++i)
				for(uint k{1}; k<=kMax-1; ++k)
				{
					variableTemporaryStorage(i,0,   k) = filterVariable(i,0,   k);
					variableTemporaryStorage(i,jMax,k) = filterVariable(i,jMax,k);
				}
			for(uint j{1}; j<=jMax-1; ++j)
				for(uint k{1}; k<=kMax-1; ++k)
				{
					variableTemporaryStorage(0,   j,k) = filterVariable(0,   j,k);
					variableTemporaryStorage(iMax,j,k) = filterVariable(iMax,j,k);
				}
			// Apply filter to interior nodes:
			for(uint i{1}; i<=iMax-1; ++i)
				for(uint j{1}; j<=jMax-1; ++j)
					for(uint k{1}; k<=kMax-1; ++k)
					{
						variableTemporaryStorage(i,j,k) = 1./2.  *   filterVariable(i,j,k)
														+ 1./12. * ( filterVariable(i+1,j,k) + filterVariable(i-1,j,k)
																   + filterVariable(i,j+1,k) + filterVariable(i,j-1,k)
																   + filterVariable(i,j,k+1) + filterVariable(i,j,k-1) );
					}
			filterVariable.dataSwap(variableTemporaryStorage);	// Swap the arrays using move-sematics (super-fast)
		}
}

// Computes an RK4 step (slope: k1, ..., k4) for the mass density. The previous state of the momentum densities, rho_u,
// rho_v and rho_w can be supplied at time t (pass the member 'rho_u' as argument 'rho_u', etc.) or as intermediate
// solutions computed in the calling function etc. Result is stored in argument 'RK4_slope'.
void Solver::compute_RK4_step_continuity(const Array3D_d& rho_u, const Array3D_d& rho_v, const Array3D_d& rho_w,
		                                    Array3D_d& RK4_slope)
{
	for(uint i{1}; i<params.NI-1; ++i)
	{
		for(uint j{1}; j<params.NJ-1; ++j)
		{
			for(uint k{1}; k<params.NK-1; ++k)
			{
				RK4_slope(i,j,k) = - ( rho_u(i+1, j  , k  )-rho_u(i-1, j  , k  ) ) / (2*dx)
							       - ( rho_v(i  , j+1, k  )-rho_v(i  , j-1, k  ) ) / (2*dy)
							       - ( rho_w(i  , j  , k+1)-rho_w(i  , j  , k-1) ) / (2*dz);
			}
		}
	}
}

// Computes an RK4 step (slope: k1, ..., k4) for the x-momentum density. The previous rho_u can be
// supplied at time t (pass the member 'rho_u' as argument 'rho_u') or as an intermediate solution
// computed in the calling function etc. Result is stored in argument 'RK4_slope'.
void Solver::compute_RK4_step_xMomentum(const Array3D_d& rho_u, Array3D_d& RK4_slope)
{
	double neg_1_div_2dx   = -1. / ( 2 * dx );       // <- The parts of the flux/residual evaluation
	double neg_1_div_2dy   = -1. / ( 2 * dy );       //    that don't change from node to node are
	double neg_1_div_2dz   = -1. / ( 2 * dz );       //    pre-computed, outside the loop.
	double pos_2_div_3dxdx =  2. / ( 3 * dx * dx );
	double pos_1_div_2dydy =  1. / ( 2 * dy * dy );
	double pos_1_div_2dzdz =  1. / ( 2 * dz * dz );
	double neg_1_div_6dxdy = -1. / ( 6 * dx * dy );
	double pos_1_div_4dxdy =  1. / ( 4 * dx * dy );
	double neg_1_div_6dxdz = -1. / ( 6 * dx * dz );
	double pos_1_div_4dxdz =  1. / ( 4 * dx * dz );

	for(uint i{1}; i<params.NI-1; ++i)
	{
		for(uint j{1}; j<params.NJ-1; ++j)
		{
			for(uint k{1}; k<params.NK-1; ++k)
			{
				RK4_slope(i,j,k) = neg_1_div_2dx   * ( rho_u(i+1,j,k) * u(i+1,j,k) + p(i+1,j,k) - rho_u(i-1,j,k) * u(i-1,j,k) - p(i-1,j,k) )
						         + neg_1_div_2dy   * ( rho_u(i,j+1,k) * v(i,j+1,k) - rho_u(i,j-1,k) * v(i,j-1,k) )
								 + neg_1_div_2dz   * ( rho_u(i,j,k+1) * w(i,j,k+1) - rho_u(i,j,k-1) * w(i,j,k-1) )
								 + pos_2_div_3dxdx * ( (mu(i+1,j,k) + mu(i,j,k)) * (u(i+1,j,k) - u(i  ,j,k))
										             - (mu(i-1,j,k) + mu(i,j,k)) * (u(i  ,j,k) - u(i-1,j,k)) )
								 + pos_1_div_2dydy * ( (mu(i,j+1,k) + mu(i,j,k)) * (u(i,j+1,k) - u(i,j  ,k))
										             - (mu(i,j-1,k) + mu(i,j,k)) * (u(i,j  ,k) - u(i,j-1,k)) )
								 + pos_1_div_2dzdz * ( (mu(i,j,k+1) + mu(i,j,k)) * (u(i,j,k+1) - u(i,j,k  ))
										             - (mu(i,j,k-1) + mu(i,j,k)) * (u(i,j,k  ) - u(i,j,k-1)) )
								 + neg_1_div_6dxdy * ( mu(i+1,j,k) * (v(i+1,j+1,k) - v(i+1,j-1,k)) - mu(i-1,j,k) * (v(i-1,j+1,k) - v(i-1,j-1,k)) )
								 + pos_1_div_4dxdy * ( mu(i,j+1,k) * (v(i+1,j+1,k) - v(i-1,j+1,k)) - mu(i,j-1,k) * (v(i+1,j-1,k) - v(i-1,j-1,k)) )
								 + neg_1_div_6dxdz * ( mu(i+1,j,k) * (w(i+1,j,k+1) - w(i+1,j,k-1)) - mu(i-1,j,k) * (w(i-1,j,k+1) - w(i-1,j,k-1)) )
								 + pos_1_div_4dxdz * ( mu(i,j,k+1) * (w(i+1,j,k+1) - w(i-1,j,k+1)) - mu(i,j,k-1) * (w(i+1,j,k-1) - w(i-1,j,k-1)) ) ;
			}
		}
	}
}

// Computes an RK4 step (slope: k1, ..., k4) for the y-momentum density. The previous rho_v can be
// supplied at time t (pass the member 'rho_v' as argument 'rho_v') or as an intermediate solution
// computed in the calling function etc. Result is stored in argument 'RK4_slope'.
void Solver::compute_RK4_step_yMomentum(const Array3D_d& rho_v, Array3D_d& RK4_slope)
{
	double neg_1_div_2dx   = -1. / ( 2 * dx );       // <- The parts of the flux/residual evaluation
	double neg_1_div_2dy   = -1. / ( 2 * dy );       //    that don't change from node to node are
	double neg_1_div_2dz   = -1. / ( 2 * dz );       //    pre-computed, outside the loop.
	double pos_1_div_2dxdx =  1. / ( 2 * dx * dx );
	double pos_2_div_3dydy =  2. / ( 3 * dy * dy );
	double pos_1_div_2dzdz =  1. / ( 2 * dz * dz );
	double neg_1_div_6dxdy = -1. / ( 6 * dx * dy );
	double pos_1_div_4dydx =  1. / ( 4 * dy * dx );
	double neg_1_div_6dzdy = -1. / ( 6 * dz * dy );
	double pos_1_div_4dydz =  1. / ( 4 * dy * dz );

	for(uint i{1}; i<params.NI-1; ++i)
	{
		for(uint j{1}; j<params.NJ-1; ++j)
		{
			for(uint k{1}; k<params.NK-1; ++k)
			{
				RK4_slope(i,j,k) = neg_1_div_2dx   * ( rho_v(i+1,j,k) * u(i+1,j,k) - rho_v(i-1,j,k) * u(i-1,j,k) )
						         + neg_1_div_2dy   * ( rho_v(i,j+1,k) * v(i,j+1,k) + p(i,j+1,k) - rho_v(i,j-1,k) * v(i,j-1,k) - p(i,j-1,k) )
								 + neg_1_div_2dz   * ( rho_v(i,j,k+1) * w(i,j,k+1) - rho_v(i,j,k-1) * w(i,j,k-1) )
								 + pos_1_div_2dxdx * ( (mu(i+1,j,k) + mu(i,j,k)) * (v(i+1,j,k) - v(i  ,j,k))
										             - (mu(i-1,j,k) + mu(i,j,k)) * (v(i  ,j,k) - v(i-1,j,k)) )
								 + pos_2_div_3dydy * ( (mu(i,j+1,k) + mu(i,j,k)) * (v(i,j+1,k) - v(i,j  ,k))
										             - (mu(i,j-1,k) + mu(i,j,k)) * (v(i,j  ,k) - v(i,j-1,k)) )
								 + pos_1_div_2dzdz * ( (mu(i,j,k+1) + mu(i,j,k)) * (v(i,j,k+1) - v(i,j,k  ))
										             - (mu(i,j,k-1) + mu(i,j,k)) * (v(i,j,k  ) - v(i,j,k-1)) )
								 + neg_1_div_6dxdy * ( mu(i,j+1,k) * (u(i+1,j+1,k) - u(i-1,j+1,k)) - mu(i,j-1,k) * (u(i+1,j-1,k) - u(i-1,j-1,k)) )
								 + pos_1_div_4dydx * ( mu(i+1,j,k) * (u(i+1,j+1,k) - u(i+1,j-1,k)) - mu(i-1,j,k) * (u(i-1,j+1,k) - u(i-1,j-1,k)) )
								 + neg_1_div_6dzdy * ( mu(i,j+1,k) * (w(i,j+1,k+1) - w(i,j+1,k-1)) - mu(i,j-1,k) * (w(i,j-1,k+1) - w(i,j-1,k-1)) )
								 + pos_1_div_4dydz * ( mu(i,j,k+1) * (w(i,j+1,k+1) - w(i,j-1,k+1)) - mu(i,j,k-1) * (w(i,j+1,k-1) - w(i,j-1,k-1)) ) ;
			}
		}
	}
}

// Computes an RK4 step (slope: k1, ..., k4) for the z-momentum density. The previous rho_w can be
// supplied at time t (pass the member 'rho_w' as argument 'rho_w') or as an intermediate solution
// computed in the calling function etc. Result is stored in argument 'RK4_slope'.
void Solver::compute_RK4_step_zMomentum(const Array3D_d& rho_w, Array3D_d& RK4_slope)
{
	double neg_1_div_2dx   = -1. / ( 2 * dx );       // <- The parts of the flux/residual evaluation
	double neg_1_div_2dy   = -1. / ( 2 * dy );       //    that don't change from node to node are
	double neg_1_div_2dz   = -1. / ( 2 * dz );       //    pre-computed, outside the loop.
	double pos_1_div_2dxdx =  1. / ( 2 * dx * dx );
	double pos_1_div_2dydy =  1. / ( 2 * dy * dy );
	double pos_2_div_3dzdz =  2. / ( 3 * dz * dz );
	double pos_1_div_4dzdx =  1. / ( 4 * dz * dx );
	double neg_1_div_6dxdz = -1. / ( 6 * dx * dz );
	double pos_1_div_4dzdy =  1. / ( 4 * dz * dy );
	double neg_1_div_6dydz = -1. / ( 6 * dy * dz );

	for(uint i{1}; i<params.NI-1; ++i)
	{
		for(uint j{1}; j<params.NJ-1; ++j)
		{
			for(uint k{1}; k<params.NK-1; ++k)
			{
				RK4_slope(i,j,k) = neg_1_div_2dx   * ( rho_w(i+1,j,k) * u(i+1,j,k) - rho_w(i-1,j,k) * u(i-1,j,k) )
						         + neg_1_div_2dy   * ( rho_w(i,j+1,k) * v(i,j+1,k) - rho_w(i,j-1,k) * v(i,j-1,k) )
								 + neg_1_div_2dz   * ( rho_w(i,j,k+1) * w(i,j,k+1) + p(i,j,k+1) - rho_w(i,j,k-1) * w(i,j,k-1) - p(i,j,k-1) )
								 + pos_1_div_2dxdx * ( (mu(i+1,j,k) + mu(i,j,k)) * (w(i+1,j,k) - w(i  ,j,k))
										             - (mu(i-1,j,k) + mu(i,j,k)) * (w(i  ,j,k) - w(i-1,j,k)) )
								 + pos_1_div_2dydy * ( (mu(i,j+1,k) + mu(i,j,k)) * (w(i,j+1,k) - w(i,j  ,k))
										             - (mu(i,j-1,k) + mu(i,j,k)) * (w(i,j  ,k) - w(i,j-1,k)) )
								 + pos_2_div_3dzdz * ( (mu(i,j,k+1) + mu(i,j,k)) * (w(i,j,k+1) - w(i,j,k  ))
										             - (mu(i,j,k-1) + mu(i,j,k)) * (w(i,j,k  ) - w(i,j,k-1)) )
								 + pos_1_div_4dzdx * ( mu(i+1,j,k) * (u(i+1,j,k+1) - u(i+1,j,k-1)) - mu(i-1,j,k) * (u(i-1,j,k+1) - u(i-1,j,k-1)) )
								 + neg_1_div_6dxdz * ( mu(i,j,k+1) * (u(i+1,j,k+1) - u(i-1,j,k+1)) - mu(i,j,k-1) * (u(i+1,j,k-1) - u(i-1,j,k-1)) )
								 + pos_1_div_4dzdy * ( mu(i,j+1,k) * (v(i,j+1,k+1) - v(i,j+1,k-1)) - mu(i,j-1,k) * (v(i,j-1,k+1) - v(i,j-1,k-1)) )
								 + neg_1_div_6dydz * ( mu(i,j,k+1) * (v(i,j+1,k+1) - v(i,j-1,k+1)) - mu(i,j,k-1) * (v(i,j+1,k-1) - v(i,j-1,k-1)) ) ;
			}
		}
	}
}

// Computes an RK4 step (slope: k1, ..., k4) for the total specific energy. The previous energy can
// be supplied at time t (pass the member 'E' as argument 'E') or as an intermediate solution
// computed in the calling function etc. Result is stored in argument 'RK4_slope'.
void Solver::compute_RK4_step_energy(const Array3D_d& E, Array3D_d& RK4_slope)
{
	double neg_1_div_2dx   = -1. / ( 2 * dx );       // <- The parts of the flux/residual evaluation
	double neg_1_div_2dy   = -1. / ( 2 * dy );       //    that don't change from node to node are
	double neg_1_div_2dz   = -1. / ( 2 * dz );       //    pre-computed, outside the loop.
	double pos_2_div_3dxdx =  2. / ( 3 * dx * dx );
	double pos_1_div_2dxdx =  1. / ( 2 * dx * dx );
	double neg_1_div_6dydx = -1. / ( 6 * dy * dx );
	double neg_1_div_6dzdx = -1. / ( 6 * dz * dx );
	double pos_1_div_4dydx =  1. / ( 4 * dy * dx );
	double pos_1_div_4dzdx =  1. / ( 4 * dz * dx );
	double pos_1_div_2dydy =  1. / ( 2 * dy * dy );
	double pos_2_div_3dydy =  2. / ( 3 * dy * dy );
	double neg_1_div_6dzdy = -1. / ( 6 * dz * dy );
	double pos_1_div_4dzdy =  1. / ( 4 * dz * dy );
	double pos_1_div_2dzdz =  1. / ( 2 * dz * dz );
	double pos_2_div_3dzdz =  2. / ( 3 * dz * dz );
	double rho_H_0 = 1. / ( params.Gamma - 1 );     // <- Dimensionless reference enthalpy

	for(uint i{1}; i<params.NI-1; ++i)
	{
		for(uint j{1}; j<params.NJ-1; ++j)
		{
			for(uint k{1}; k<params.NK-1; ++k)
			{
				double mu_P{mu(i,j,k)}, mu_W{mu(i-1,j,k)}, mu_E{mu(i+1,j,k)}, mu_S{mu(i,j-1,k)}, mu_N{mu(i,j+1,k)}, mu_D{mu(i,j,k-1)}, mu_U{mu(i,j,k+1)};
				double  u_P{ u(i,j,k)},  u_W{ u(i-1,j,k)},  u_E{ u(i+1,j,k)},  u_S{ u(i,j-1,k)},  u_N{ u(i,j+1,k)},  u_D{ u(i,j,k-1)},  u_U{ u(i,j,k+1)};
				double  v_P{ v(i,j,k)},  v_W{ v(i-1,j,k)},  v_E{ v(i+1,j,k)},  v_S{ v(i,j-1,k)},  v_N{ v(i,j+1,k)},  v_D{ v(i,j,k-1)},  v_U{ v(i,j,k+1)};
				double  w_P{ w(i,j,k)},  w_W{ w(i-1,j,k)},  w_E{ w(i+1,j,k)},  w_S{ w(i,j-1,k)},  w_N{ w(i,j+1,k)},  w_D{ w(i,j,k-1)},  w_U{ w(i,j,k+1)};
				double kappa_P{kappa(i,j,k)}, T_P{T(i,j,k)};

				double inviscidFluxes = neg_1_div_2dx * ( (rho_H_0 + E(i+1,j,k) + p(i+1,j,k)) * u_E - (rho_H_0 + E(i-1,j,k) + p(i-1,j,k)) * u_W )
						              + neg_1_div_2dy * ( (rho_H_0 + E(i,j+1,k) + p(i,j+1,k)) * v_N - (rho_H_0 + E(i,j-1,k) + p(i,j-1,k)) * v_S )
									  + neg_1_div_2dz * ( (rho_H_0 + E(i,j,k+1) + p(i,j,k+1)) * w_U - (rho_H_0 + E(i,j,k-1) + p(i,j,k-1)) * w_D );

				double viscousFluxX = pos_2_div_3dxdx * ( (mu_E*u_E + mu_P*u_P) * (u_E - u_P) - (mu_W*u_W + mu_P*u_P) * (u_P - u_W) )
						            + pos_1_div_2dxdx * ( (mu_E*v_E + mu_P*v_P) * (v_E - v_P) - (mu_W*v_W + mu_P*v_P) * (v_P - v_W) )
									+ pos_1_div_2dxdx * ( (mu_E*w_E + mu_P*w_P) * (w_E - w_P) - (mu_W*w_W + mu_P*w_P) * (w_P - w_W) )
									+ neg_1_div_6dydx * ( mu_E * u_E * (v(i+1,j+1,k) - v(i+1,j-1,k)) - mu_W * u_W * (v(i-1,j+1,k) - v(i-1,j-1,k)) )
									+ neg_1_div_6dzdx * ( mu_E * u_E * (w(i+1,j,k+1) - w(i+1,j,k-1)) - mu_W * u_W * (w(i-1,j,k+1) - w(i-1,j,k-1)) )
									+ pos_1_div_4dydx * ( mu_E * v_E * (u(i+1,j+1,k) - u(i+1,j-1,k)) - mu_W * v_W * (u(i-1,j+1,k) - u(i-1,j-1,k)) )
									+ pos_1_div_4dzdx * ( mu_E * w_E * (u(i+1,j,k+1) - u(i+1,j,k-1)) - mu_W * w_W * (u(i-1,j,k+1) - u(i-1,j,k-1)) );

				double viscousFluxY = pos_1_div_2dydy * ( (mu_N*u_N + mu_P*u_P) * (u_N - u_P) - (mu_S*u_S + mu_P*u_P) * (u_P - u_S) )
						            + pos_2_div_3dydy * ( (mu_N*v_N + mu_P*v_P) * (v_N - v_P) - (mu_S*v_S + mu_P*v_P) * (v_P - v_S) )
									+ pos_1_div_2dydy * ( (mu_N*w_N + mu_P*w_P) * (w_N - w_P) - (mu_S*w_S + mu_P*w_P) * (w_P - w_S) )
									+ neg_1_div_6dydx * ( mu_N * v_N * (u(i+1,j+1,k) - u(i-1,j+1,k)) - mu_S * v_S * (u(i+1,j-1,k) - u(i-1,j-1,k)) )
									+ neg_1_div_6dzdy * ( mu_N * v_N * (w(i,j+1,k+1) - w(i,j+1,k-1)) - mu_S * v_S * (w(i,j-1,k+1) - w(i,j-1,k-1)) )
									+ pos_1_div_4dydx * ( mu_N * u_N * (v(i+1,j+1,k) - v(i-1,j+1,k)) - mu_S * u_S * (v(i+1,j-1,k) - v(i-1,j-1,k)) )
									+ pos_1_div_4dzdy * ( mu_N * w_N * (v(i,j+1,k+1) - v(i,j+1,k-1)) - mu_S * w_S * (v(i,j-1,k+1) - v(i,j-1,k-1)) );

				double viscousFluxZ = pos_1_div_2dzdz * ( (mu_U*u_U + mu_P*u_P) * (u_U - u_P) - (mu_D*u_D + mu_P*u_P) * (u_P - u_D) )
						            + pos_1_div_2dzdz * ( (mu_U*v_U + mu_P*v_P) * (v_U - v_P) - (mu_D*v_D + mu_P*v_P) * (v_P - v_D) )
									+ pos_2_div_3dzdz * ( (mu_U*w_U + mu_P*w_P) * (w_U - w_P) - (mu_D*w_D + mu_P*w_P) * (w_P - w_D) )
									+ neg_1_div_6dzdx * ( mu_U * w_U * (u(i+1,j,k+1) - u(i-1,j,k+1)) - mu_D * w_D * (u(i+1,j,k-1) - u(i-1,j,k-1)) )
									+ neg_1_div_6dzdy * ( mu_U * w_U * (v(i,j+1,k+1) - v(i,j-1,k+1)) - mu_D * w_D * (v(i,j+1,k-1) - v(i,j-1,k-1)) )
									+ pos_1_div_4dzdx * ( mu_U * u_U * (w(i+1,j,k+1) - w(i-1,j,k+1)) - mu_D * u_D * (w(i+1,j,k-1) - w(i-1,j,k-1)) )
									+ pos_1_div_4dzdy * ( mu_U * v_U * (w(i,j+1,k+1) - w(i,j-1,k+1)) - mu_D * v_D * (w(i,j+1,k-1) - w(i,j-1,k-1)) );

				double heatFluxes = pos_1_div_2dxdx * ( (kappa(i+1,j,k) + kappa_P) * (T(i+1,j,k) - T_P) - (kappa(i-1,j,k) + kappa_P) * (T_P - T(i-1,j,k)) )
						          + pos_1_div_2dydy * ( (kappa(i,j+1,k) + kappa_P) * (T(i,j+1,k) - T_P) - (kappa(i,j-1,k) + kappa_P) * (T_P - T(i,j-1,k)) )
								  + pos_1_div_2dzdz * ( (kappa(i,j,k+1) + kappa_P) * (T(i,j,k+1) - T_P) - (kappa(i,j,k-1) + kappa_P) * (T_P - T(i,j,k-1)) );

				RK4_slope(i,j,k) = inviscidFluxes + viscousFluxX + viscousFluxY + viscousFluxZ + heatFluxes;
			}
		}
	}
}

// Evaluates an intermediate solution of the given conserved variable, using the given RK4 step (slope)
// and the appropriate time increment (dt/2 if k1 or k2 is given, dt if k3 is given).
// Result is stored in 'intermSolution'. Only writes interior nodes.
void Solver::computeIntermediateSolution(const Array3D_d& conservedVar, const Array3D_d& RK4_slope,
		                                       Array3D_d& intermSolution, double timeIncrement)
{
	for (uint i{1}; i<params.NI-1; ++i)
		for (uint j{1}; j<params.NJ-1; ++j)
			for (uint k{1}; k<params.NK-1; ++k)
				intermSolution(i,j,k) = conservedVar(i,j,k) + timeIncrement * RK4_slope(i,j,k);
}

// Updates the primitive variables and transport properties, using the specified conserved variables.
// The conserved variables can be given at time t (just pass the data-members rho, rho_u, ..., E),
// or as intermediate solutions computed in the calling function etc. Only writes interior nodes!
void Solver::updatePrimitiveVariables(const Array3D_d& rho, const Array3D_d& rho_u, const Array3D_d& rho_v,
		                                 const Array3D_d& rho_w, const Array3D_d& E)
{
	double gammaMinusOne = params.Gamma - 1;
	double sutherlands_Sc = params.T_0 / params.sutherlands_C2;
	double ScPlusOne = 1 + sutherlands_Sc;
	double prandtlFactor = 1 / ( gammaMinusOne * params.Pr );

	for (uint i{1}; i<params.NI-1; ++i)
		for (uint j{1}; j<params.NJ-1; ++j)
			for (uint k{1}; k<params.NK-1; ++k)
				u(i,j,k) = rho_u(i,j,k) / ( rho(i,j,k) + 1 );
	for (uint i{1}; i<params.NI-1; ++i)
		for (uint j{1}; j<params.NJ-1; ++j)
			for (uint k{1}; k<params.NK-1; ++k)
				v(i,j,k) = rho_v(i,j,k) / ( rho(i,j,k) + 1 );
	for (uint i{1}; i<params.NI-1; ++i)
		for (uint j{1}; j<params.NJ-1; ++j)
			for (uint k{1}; k<params.NK-1; ++k)
				w(i,j,k) = rho_w(i,j,k) / ( rho(i,j,k) + 1 );
	for (uint i{1}; i<params.NI-1; ++i)
		for (uint j{1}; j<params.NJ-1; ++j)
			for (uint k{1}; k<params.NK-1; ++k)
				p(i,j,k) = gammaMinusOne * ( E(i,j,k) - (rho(i,j,k)+1)/2 * (u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k)) );
	for (uint i{1}; i<params.NI-1; ++i)
		for (uint j{1}; j<params.NJ-1; ++j)
			for (uint k{1}; k<params.NK-1; ++k)
				T(i,j,k) = ( params.Gamma * p(i,j,k) - rho(i,j,k) ) / ( 1+rho(i,j,k) );
	for (uint i{1}; i<params.NI-1; ++i)
		for (uint j{1}; j<params.NJ-1; ++j)
			for (uint k{1}; k<params.NK-1; ++k)
				mu(i,j,k) = pow( 1+T(i,j,k), 1.5 ) * ScPlusOne / ( params.Re*( T(i,j,k) + ScPlusOne ) );
	for (uint i{1}; i<params.NI-1; ++i)
		for (uint j{1}; j<params.NJ-1; ++j)
			for (uint k{1}; k<params.NK-1; ++k)
				kappa(i,j,k) = mu(i,j,k) * prandtlFactor;
}

// Applies channel flow boundary conditions, by setting the variables in the outer nodes.
// Two walls with no-slip condition, no heat flux and no pressure gradient. Periodic BC in one direction.
// Inlet: Parabolic velocity profile from incompressible analytic solution is gradually applied from stagnation. Zero temperature perturbation, pressure extrapolated.
// Outlet: Zero pressure perturbation. Density and momentums extrapolated.
void Solver::applyInjectionBC_risingInletVelocityChannelFlow(Array3D_d& rho, Array3D_d& rho_u, Array3D_d& rho_v, Array3D_d& rho_w, Array3D_d& E, double time)
{
	uint iMax{params.NI - 1}, jMax{params.NJ - 1}, kMax{params.NK - 1};
	double gammaMinusOne = params.Gamma - 1;
	double sutherlands_Sc = params.T_0 / params.sutherlands_C2;
	double ScPlusOne = 1 + sutherlands_Sc;
	double prandtlFactor = 1 / ( gammaMinusOne * params.Pr );
	for (uint j{1}; j<=jMax-1; ++j)
		for (uint k{1}; k<=kMax-1; ++k)
		{
			double z = k * dz;
			double u_parabolic = 4. * params.M_0 * ( z - z*z );
			u(0,j,k) = min( 0.2*t*u_parabolic, u_parabolic );
			v(0,j,k) = 0;
			w(0,j,k) = 0;
			p(0,j,k) = 2*p(1,j,k) - p(2,j,k);
			T(0,j,k) = 0;
			getDerivedVariables_atPoint(u(0,j,k), v(0,j,k), w(0,j,k), p(0,j,k), T(0,j,k), rho(0,j,k), rho_u(0,j,k), rho_v(0,j,k), rho_w(0,j,k), E(0,j,k), mu(0,j,k), kappa(0,j,k));

			rho  (iMax, j, k) = 2*rho  (iMax-1, j, k) - rho  (iMax-2, j, k);
			rho_u(iMax, j, k) = 2*rho_u(iMax-1, j, k) - rho_u(iMax-2, j, k);
			rho_v(iMax, j, k) = 2*rho_v(iMax-1, j, k) - rho_v(iMax-2, j, k);
			rho_w(iMax, j, k) = 2*rho_w(iMax-1, j, k) - rho_w(iMax-2, j, k);
			p(iMax,j,k) = 0;
			u(iMax,j,k) = rho_u(iMax,j,k) / (1+rho(iMax,j,k));
			v(iMax,j,k) = rho_v(iMax,j,k) / (1+rho(iMax,j,k));
			w(iMax,j,k) = rho_w(iMax,j,k) / (1+rho(iMax,j,k));
			T(iMax,j,k) = (params.Gamma * p(iMax,j,k) - rho(iMax,j,k))/(1 + rho(iMax,j,k));
			E(iMax,j,k) = p(iMax,j,k) / gammaMinusOne + (1 + rho(iMax,j,k))/2 * ( pow(u(iMax,j,k),2) + pow(v(iMax,j,k),2) + pow(w(iMax,j,k),2) );
			mu(iMax,j,k) = pow( 1+T(iMax,j,k), 1.5 ) * ScPlusOne / ( params.Re*( T(iMax,j,k) + ScPlusOne ) );
			kappa(iMax,j,k) = mu(iMax,j,k) * prandtlFactor;
		}
	for (uint i{0}; i<=iMax; ++i)
		for (uint k{1}; k<=kMax-1; ++k)
		{
			u  (i,0,k) = u(i,jMax-1,k);
			v  (i,0,k) = v(i,jMax-1,k);
			w  (i,0,k) = w(i,jMax-1,k);
			p  (i,0,k) = p(i,jMax-1,k);
			T  (i,0,k) = T(i,jMax-1,k);
			getDerivedVariables_atPoint(u(i,0,k), v(i,0,k), w(i,0,k), p(i,0,k), T(i,0,k), rho(i,0,k), rho_u(i,0,k), rho_v(i,0,k), rho_w(i,0,k), E(i,0,k), mu(i,0,k), kappa(i,0,k));

			u  (i, jMax, k) = u(i, 1, k);
			v  (i, jMax, k) = v(i, 1, k);
			w  (i, jMax, k) = w(i, 1, k);
			p  (i, jMax, k) = p(i, 1, k);
			T  (i, jMax, k) = T(i, 1, k);
			getDerivedVariables_atPoint(u(i, jMax, k), v(i, jMax, k), w(i, jMax, k), p(i, jMax, k), T(i, jMax, k),
					rho(i, jMax, k), rho_u(i, jMax, k), rho_v(i, jMax, k), rho_w(i, jMax, k), E(i, jMax, k), mu(i, jMax, k), kappa(i, jMax, k));
		}
	for (uint i{0}; i<=iMax; ++i)
		for (uint j{0}; j<=jMax; ++j)
		{

			u(i,j,0) = 0;	// No-slip
			v(i,j,0) = 0;
			w(i,j,0) = 0;
			p(i,j,0) = (4*p(i,j,1) - p(i,j,2))/3;	// 2nd order Neumann BC (zero gradient)
			T(i,j,0) = (4*T(i,j,1) - T(i,j,2))/3;
			getDerivedVariables_atPoint(u(i,j,0), v(i,j,0), w(i,j,0), p(i,j,0), T(i,j,0), rho(i,j,0), rho_u(i,j,0), rho_v(i,j,0), rho_w(i,j,0), E(i,j,0), mu(i,j,0), kappa(i,j,0));

			u(i,j,kMax) = 0;	// No-slip
			v(i,j,kMax) = 0;
			w(i,j,kMax) = 0;
			p(i,j,kMax) = (4*p(i,j,kMax-1) - p(i,j,kMax-2))/3;	// 2nd order Neumann BC (zero gradient)
			T(i,j,kMax) = (4*T(i,j,kMax-1) - T(i,j,kMax-2))/3;
			getDerivedVariables_atPoint(u(i,j,kMax), v(i,j,kMax), w(i,j,kMax), p(i,j,kMax), T(i,j,kMax),
					rho(i,j,kMax), rho_u(i,j,kMax), rho_v(i,j,kMax), rho_w(i,j,kMax), E(i,j,kMax), mu(i,j,kMax), kappa(i,j,kMax));
		}
}

// Uses all the intermediate solutions for a variable, k1,2,3,4, to advance the variable from t to t+dt
void Solver::compute_RK4_final_step(const Array3D_d& k1, const Array3D_d& k2,
									const Array3D_d& k3, const Array3D_d& k4,
									const Array3D_d& conservedVar_old, Array3D_d& conservedVar_new)
{
	double dtFactor{dt / 6};
	for (uint i{1}; i<params.NI-1; ++i)
		for (uint j{1}; j<params.NJ-1; ++j)
			for (uint k{1}; k<params.NK-1; ++k)
			{
				double residualValue = dtFactor*( k1(i,j,k) + 2*k2(i,j,k) + 2*k3(i,j,k) + k4(i,j,k) );
				conservedVar_new(i,j,k) = conservedVar_old(i,j,k) + residualValue;
			}
}

// Compute the norm of change in the conserved variables, and store in the history vectors.
void Solver::computeNorms_conservedVariables()
{
	normHistory_rho  .push_back( getNormOfChange(rho  , interm_rho  ) );
	normHistory_rho_u.push_back( getNormOfChange(rho_u, interm_rho_u) );
	normHistory_rho_v.push_back( getNormOfChange(rho_v, interm_rho_v) );
	normHistory_rho_w.push_back( getNormOfChange(rho_w, interm_rho_w) );
	normHistory_E    .push_back( getNormOfChange(E    , interm_E    ) );
}

// Compute the 2-norm of the difference between two arrys. Intended to monitor the change between two consecutive time levels.
// E.g. to check convergence of solution.
double Solver::getNormOfChange(const Array3D_d& oldValue, const Array3D_d& newValue)
{
	double sumOfSquaredChanges = 0;
	uint numberOfMeshNodes = params.NI * params.NJ * params.NK;
	for(uint i=0; i<numberOfMeshNodes; ++i)
		sumOfSquaredChanges += pow(oldValue(i)-newValue(i), 2);
	double normOfChange = sqrt( dx*dy*dz * sumOfSquaredChanges );
	return normOfChange;
}

// Swap the contents of all the arrays of conserved variables and the intermediate arrays, by move-semantics.
// This operation is super fast and needs no extra copy. Only the ownership of the data is changed.
void Solver::swapConservedVariables()
{
	rho  .dataSwap(interm_rho  );
	rho_u.dataSwap(interm_rho_u);
	rho_v.dataSwap(interm_rho_v);
	rho_w.dataSwap(interm_rho_w);
	E    .dataSwap(interm_E    );
}

// Checks whether it is time to store the solution to disk, or write out a status report to screen and does it if appropriate.
// It will save solution if the next timestep would take the solution past the save-time.
// Thus, in general it saves too early, but the deviation is dt at most.
void Solver::processOutput(Clock& statusReportTimer)
{
	double timeSinceSave = fmod(t, params.save_period);
	if ( timeSinceSave + dt >= params.save_period && params.save_intervals )
		storeCurrentSolution_csv();
	double timeSinceStatusReport = statusReportTimer.getElapsedTime().asSeconds();
	if ( timeSinceStatusReport >= params.statusReportInterval )
	{
		writeStatusReport_toScreen();
		statusReportTimer.restart();
	}
}

// Store selected variables from the solution at current time level, using the format(s) specified
void Solver::storeCurrentSolution_csv()
{
	if(params.saveForParaview)
		storeCurrentSolution_csv_paraview();
	if (params.saveForMatlab)
		storeCurrentSolution_csv_matlab();
	++savedSolutions;
	outputTimes.push_back(t);
}

// Writes a .csv file with the specified flow variables at the current time level.
// The file is formatted to fit how ParaView wants to import it, with coordinates in the leftmost column.
// Only the flow variables selected in the ConfigFile are saved.
void Solver::storeCurrentSolution_csv_paraview()
{
	vector<Array3D_d*> flowVariables = getPlotVariables();
	if ( flowVariables.empty() )
		return;

	ofstream outputFile;
	string filename = "output/out.csv." + to_string(savedSolutions);
	outputFile.open( filename );
	if ( !outputFile )
	{
		cout << "Could not open file:  " + filename << endl
		     << "Solution was not saved. You may have to move the 'output' folder to the location of the source files or to the executable itself, depending on whether you run the executable through an IDE or by itself." << endl;
		return;
	}
	outputFile << get_csvHeaderString();

	for (uint i{0}; i<params.NI; ++i)
		for (uint j{0}; j<params.NJ; ++j)
			for (uint k{0}; k<params.NK; ++k)
			{
				outputFile << endl;
				double x{ i*dx }, y{ j*dy }, z{ k*dz};
				outputFile << x << ", " << y << ", " << z;
				for (Array3D_d* flowVar : flowVariables)
					outputFile << ", " << (*flowVar)(i,j,k);
			}
	outputFile.close();
}

// Writes multiple .csv files with the specified flow variables at the current time level.
// The file is formatted as a table or matrix, with the values from a specified 2D plane. Matlab
// can read this file to a matrix. Only the flow variables selected in the ConfigFile are saved.
void Solver::storeCurrentSolution_csv_matlab()
{
	vector<Array3D_d*> flowVariables = getPlotVariables();		// Get pointers to the arrays with data to save
	vector<string> variableFileNames = getVariableFileNames();	// Get a vector with the names of the variables

	for(uint i=0; i<flowVariables.size(); ++i)	// Let i go from zero to no. of variables to save.
	{
		ofstream outputFile;
		string filename = "output/" + variableFileNames.at(i) + "_" + to_string(savedSolutions) + ".csv";
		outputFile.open( filename );
		if ( !outputFile )
		{
			cout << "Could not open file:  " + filename << endl
			     << "Solution was not saved. You may have to move the 'output' folder to the location of the source files or to the executable itself, depending on whether you run the executable through an IDE or by itself." << endl;
			return;
		}
		writePlaneTo_csv(outputFile, flowVariables.at(i));
	}
}

// Get a vector with pointers to the flow variables that should be saved.
vector<Array3D_d*> Solver::getPlotVariables()
{
	vector<Array3D_d*> flowVariables;
	if ( params.save_rho )
		flowVariables.push_back(&rho);
	if ( params.save_rho_u )
		flowVariables.push_back(&rho_u);
	if ( params.save_rho_v )
		flowVariables.push_back(&rho_v);
	if ( params.save_rho_w )
		flowVariables.push_back(&rho_w);
	if ( params.save_E )
		flowVariables.push_back(&E);
	if ( params.save_u )
		flowVariables.push_back(&u);
	if ( params.save_v )
		flowVariables.push_back(&v);
	if ( params.save_w )
		flowVariables.push_back(&w);
	if ( params.save_p )
		flowVariables.push_back(&p);
	if ( params.save_T )
		flowVariables.push_back(&T);
	if ( params.save_mu )
		flowVariables.push_back(&mu);
	if ( params.save_kappa )
		flowVariables.push_back(&kappa);
	return flowVariables;
}

// Returns a camma-separated string with the names of the flow variables to save.
// For ParaView to read a .csv file it needs the headers on the first line.
string Solver::get_csvHeaderString()
{
	string headers = "x, y, z";
	if ( params.save_rho )
		headers += ", density";
	if ( params.save_rho_u )
		headers += ", x-momentum";
	if ( params.save_rho_v )
		headers += ", y-momentum";
	if ( params.save_rho_w )
		headers += ", z-momentum";
	if ( params.save_E )
		headers += ", total energy";
	if ( params.save_u )
		headers += ", velocity comp u";
	if ( params.save_v )
		headers += ", velocity comp v";
	if ( params.save_w )
		headers += ", velocity comp w";
	if ( params.save_p )
		headers += ", pressure";
	if ( params.save_T )
		headers += ", temperature";
	if ( params.save_mu )
		headers += ", dynamic viscosity";
	if ( params.save_kappa )
		headers += ", thermal conductivity";
	return headers;
}

// Returns a vector of strings, which are the names of the variables to save.
vector<string> Solver::getVariableFileNames()
{
	vector<string> variableNames;
	if ( params.save_rho )
		variableNames.push_back("rho");
	if ( params.save_rho_u )
		variableNames.push_back("rho_u");
	if ( params.save_rho_v )
		variableNames.push_back("rho_v");
	if ( params.save_rho_w )
		variableNames.push_back("rho_w");
	if ( params.save_E )
		variableNames.push_back("E");
	if ( params.save_u )
		variableNames.push_back("u");
	if ( params.save_v )
		variableNames.push_back("v");
	if ( params.save_w )
		variableNames.push_back("w");
	if ( params.save_p )
		variableNames.push_back("p");
	if ( params.save_T )
		variableNames.push_back("T");
	if ( params.save_mu )
		variableNames.push_back("mu");
	if ( params.save_kappa )
		variableNames.push_back("kappa");
	return variableNames;
}

// Write the values from a plane, defined in ConfigFile, of one flow variable. It's written like a comma separated table into 'outputFile'.
void Solver::writePlaneTo_csv(ofstream& outputFile, Array3D_d* flowVariable)
{
	switch(params.saveNormalAxis)
	{
	case saveNormalAxisEnum::x:
		for(uint j=0; j<params.NJ; ++j)
		{
			outputFile << (*flowVariable)(params.saveConstantIndex, j, 0);
			for(uint k=1; k<params.NK; ++k)
				outputFile << ", " << (*flowVariable)(params.saveConstantIndex, j, k);
			outputFile << endl;
		}
		break;
	case saveNormalAxisEnum::y:
		for(uint i=0; i<params.NI; ++i)
		{
			outputFile << (*flowVariable)(i, params.saveConstantIndex, 0);
			for(uint k=1; k<params.NK; ++k)
				outputFile << ", " << (*flowVariable)(i, params.saveConstantIndex, k);
			outputFile << endl;
		}
		break;
	case saveNormalAxisEnum::z:
		for(uint i=0; i<params.NI; ++i)
		{
			outputFile << (*flowVariable)(i, 0, params.saveConstantIndex);
			for(uint j=1; j<params.NJ; ++j)
				outputFile << ", " << (*flowVariable)(i, j, params.saveConstantIndex);
			outputFile << endl;
		}
		break;
	}
}

// Writes a brief report with progression, timestep size and wall clock time elapsed, etc.
// 'setprecision' is used to control no. of significant digits, default is 6.
void Solver::writeStatusReport_toScreen()
{
	cout << "Simulated time: t = " << t;
	if (params.stopCriterion == StopCriterionEnum::end_time)
	{
		double progressPercentage = t / params.t_end * 100.;
		cout << " , t_end = " << params.t_end << " ( " << setprecision(3) << progressPercentage << setprecision(6) << " % )" << endl;
	}
	else
		cout << endl;

	cout << "Time level: n = " << timeLevel;
	if (params.stopCriterion == StopCriterionEnum::timesteps)
	{
		double progressPercentage = (double)timeLevel / (double)params.stopTimeLevel * 100.;  // Cast to double to avoid integer division.
		cout << " , n_max = " << params.stopTimeLevel << " ( " << setprecision(3) << progressPercentage << setprecision(6) << " % )" << endl;
	}
	else
		cout << endl;
	cout << "Timestep size: dt = " << dt << endl;
	double inletMassFlux=0, outletMassFlux=0;
	checkMassConservation(inletMassFlux, outletMassFlux);
	cout << "Inlet mass flux: " << inletMassFlux << " , Outlet mass flux: " << outletMassFlux << " , Difference: " << inletMassFlux-outletMassFlux << endl;
	cout << "Wall clock time: " << setprecision(3) << wallClockTimer.getElapsedTime().asSeconds() << setprecision(6) << " sec" << endl << endl;
}

// Write a file with a list of the actual times when the solution was saved. For debugging.
void Solver::writeOutputTimes()
{
	ofstream timeFile;
	string filename = "output/times.dat";
	timeFile.open( filename );
	if ( !timeFile )
	{
		cout << "Could not open file:  " + filename << endl
		     << "Output times were not written to .dat file. You may have to move the 'output' folder to the location of the source files or to the executable itself, depending on whether you run the executable through an IDE or by itself." << endl;
		return;
	}
	for (double time : outputTimes)
		timeFile << time << endl;
	timeFile.close();
}

// Write files with lists of the norm of change for the conserved variables.
// TODO: This function should be customizable to choose in ConfigFile whether to log norms of particular variables. OR, find out which variable is the best indicator for convergence, and only log that one.
void Solver::writeNormHistoryFiles()
{
	vector<string> filenames = {"output/norm_rho.dat", "output/norm_rho_u.dat", "output/norm_rho_v.dat", "output/norm_rho_w.dat", "output/norm_E.dat"};
	vector<vector<double>*> normVectorPointers = {&normHistory_rho, &normHistory_rho_u, &normHistory_rho_v, &normHistory_rho_w, &normHistory_E};
	for(uint i=0; i<filenames.size(); ++i)
	{
		ofstream normFile;
		normFile.open( filenames.at(i) );
		if ( !normFile )
		{
			cout << "Could not open file:  " + filenames.at(i) << endl
				 << "Norm history was not written to .dat file. You may have to move the 'output' folder to the location of the source files or to the executable itself, depending on whether you run the executable through an IDE or by itself." << endl;
		}
		else
			for (double norm : *normVectorPointers.at(i))
				normFile << norm << endl;
		normFile.close();
	}
}

// Check whether the influx and outflux of mass are equal, by integrating momentum over inlet and outlet boundaries.
// Integration is done by double trapezoidal rule, i.e., a double loop summing the nodes in the inlet and outlet planes, wighting the outer nodes with 1/2.
void Solver::checkMassConservation(double& inFluxSum, double& outFluxSum)
{
	inFluxSum = outFluxSum = 0;
	for (uint j=0; j<params.NJ; ++j)
	{
		double inletColumnSum = 0, outletColumnSum = 0;
		inletColumnSum  += rho_u(0          ,j,0)/2 + rho_u(0          ,j,params.NK-1)/2;	// Add the boundary nodes with weigth 1/2
		outletColumnSum += rho_u(params.NI-1,j,0)/2 + rho_u(params.NI-1,j,params.NK-1)/2;
		for(uint k=1; k<params.NK-1; ++k)
		{
			inletColumnSum  += rho_u(0,j,k);
			outletColumnSum += rho_u(params.NI-1,j,k);
		}
		inletColumnSum  *= dz;
		outletColumnSum *= dz;
		double weightMultiplier = (j==0 || j==params.NJ-1)? 0.5 : 1.;
		inFluxSum  += weightMultiplier * inletColumnSum;
		outFluxSum += weightMultiplier * outletColumnSum;
	}
}






