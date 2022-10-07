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






