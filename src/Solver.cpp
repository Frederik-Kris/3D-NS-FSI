/*
 * Solver.cpp
 *
 *  Created on: Apr 28, 2021
 *      Author: frederik
 */


#include "Solver.h"

// Constructor. Takes ConfigSettings 'params', with necessary settings from config file
// Initializes mesh with 3D arrays (but not their elements)
Solver::Solver(const ConfigSettings& params) :
params(params),
mesh(params),
dt{0}
{
}

// Apply initial condition (IC) and set the time-step size based on the IC.
void Solver::initialize()
{
	mesh.setupBoundaries(params);
	mesh.categorizeNodes(params);
	applyStagnation_IC();
	updateTimeStepSize(0);
}

// Applies stagnation initial condition (IC) by setting the values of flow variables at every node
void Solver::applyStagnation_IC()
{
	cout << "Applying stagnation initial condition... ";
	double u_IC{0}, v_IC{0}, w_IC{0}, p_IC{0}, T_IC{0};                  	// <- Decide primitive variables
	PrimitiveVariablesScalars  decidedPrimitiveVariables(u_IC, v_IC, w_IC, p_IC, T_IC); // and derive the others
	ConservedVariablesScalars  derivedConservedVariables  = deriveConservedVariables(decidedPrimitiveVariables, params);
	TransportPropertiesScalars derivedTransportProperties = deriveTransportProperties(decidedPrimitiveVariables, params);
	mesh.primitiveVariables.u.setAll(decidedPrimitiveVariables.u);
	mesh.primitiveVariables.v.setAll(decidedPrimitiveVariables.v);
	mesh.primitiveVariables.w.setAll(decidedPrimitiveVariables.w);
	mesh.primitiveVariables.p.setAll(decidedPrimitiveVariables.p);
	mesh.primitiveVariables.T.setAll(decidedPrimitiveVariables.T);
	mesh.conservedVariables.rho  .setAll(derivedConservedVariables.rho);
	mesh.conservedVariables.rho_u.setAll(derivedConservedVariables.rho_u);
	mesh.conservedVariables.rho_v.setAll(derivedConservedVariables.rho_v);
	mesh.conservedVariables.rho_w.setAll(derivedConservedVariables.rho_w);
	mesh.conservedVariables.rho_E.setAll(derivedConservedVariables.rho_E);
	mesh.transportProperties.mu   .setAll(derivedTransportProperties.mu);
	mesh.transportProperties.kappa.setAll(derivedTransportProperties.kappa);
	cout << "Done" << endl << endl;
}

// Applies uniform flow initial condition (IC) by setting the values of flow variables at every node
void Solver::applyUniformFlow_IC()
{
	cout << "Applying oblique uniform flow initial condition... ";
	double velocity = params.M_0 / sqrt(3);
	double u_IC{velocity}, v_IC{velocity}, w_IC{velocity}, p_IC{0}, T_IC{0}; // <- Decide primitive variables
	PrimitiveVariablesScalars  decidedPrimitiveVariables(u_IC, v_IC, w_IC, p_IC, T_IC);// and derive the others
	ConservedVariablesScalars  derivedConservedVariables  = deriveConservedVariables(decidedPrimitiveVariables, params);
	TransportPropertiesScalars derivedTransportProperties = deriveTransportProperties(decidedPrimitiveVariables, params);
	mesh.primitiveVariables.u.setAll(decidedPrimitiveVariables.u);
	mesh.primitiveVariables.v.setAll(decidedPrimitiveVariables.v);
	mesh.primitiveVariables.w.setAll(decidedPrimitiveVariables.w);
	mesh.primitiveVariables.p.setAll(decidedPrimitiveVariables.p);
	mesh.primitiveVariables.T.setAll(decidedPrimitiveVariables.T);
	mesh.conservedVariables.rho  .setAll(derivedConservedVariables.rho);
	mesh.conservedVariables.rho_u.setAll(derivedConservedVariables.rho_u);
	mesh.conservedVariables.rho_v.setAll(derivedConservedVariables.rho_v);
	mesh.conservedVariables.rho_w.setAll(derivedConservedVariables.rho_w);
	mesh.conservedVariables.rho_E.setAll(derivedConservedVariables.rho_E);
	mesh.transportProperties.mu   .setAll(derivedTransportProperties.mu);
	mesh.transportProperties.kappa.setAll(derivedTransportProperties.kappa);
	cout << "Done" << endl << endl;
}

// Set the time step size (dt) as large as possible, within the stability criterion.
// Computes two time step sizes, using the inviscid CFL condition, and the viscous von Neumann condition
// and assigns the smallest one to dt (strictest criterion).
void Solver::updateTimeStepSize(double t)
{
	double dtConvective = getInviscidTimeStepLimit();
	double dtInviscid   = getViscousTimeStepLimit();
	// Choose the strictest criterion. If that will take us past the end-time, then adapt dt to hit t_end exactly:
	dt = min( dtConvective, dtInviscid );
	if ( params.stopCriterion == StopCriterionEnum::end_time && t + dt > params.t_end )
		dt = params.t_end - t;
}

// Find the inviscid time step limit (CFL condition):
double Solver::getInviscidTimeStepLimit()
{
	const Array3D_d& p{mesh.primitiveVariables.p}, &rho{mesh.conservedVariables.rho},
		&u{mesh.primitiveVariables.u}, &v{mesh.primitiveVariables.v}, &w{mesh.primitiveVariables.w};
	double maxSpectralRadiusX{0}, maxSpectralRadiusY{0}, maxSpectralRadiusZ{0};
	for (size_t i{0}; i<mesh.nNodesTotal; ++i)
	{
		double c_i = sqrt( (1 + params.Gamma * p(i)) / (1 + rho(i)) );    //<- Speed of sound at node i
		maxSpectralRadiusX = max( maxSpectralRadiusX, c_i + fabs(u(i)) );
		maxSpectralRadiusY = max( maxSpectralRadiusY, c_i + fabs(v(i)) );
		maxSpectralRadiusZ = max( maxSpectralRadiusZ, c_i + fabs(w(i)) );
	}
	return fabs(params.convStabilityLimit) / ( maxSpectralRadiusX / mesh.dx
			                                 + maxSpectralRadiusY / mesh.dy
											 + maxSpectralRadiusZ / mesh.dz );
}

// Find the viscous time step limit (von Neumann condition):
double Solver::getViscousTimeStepLimit()
{
	const Array3D_d& mu{mesh.transportProperties.mu}, &rho{mesh.conservedVariables.rho};
	double viscosityModifier = max( 4./3, params.Gamma / params.Pr );
	double maxViscosityFactor{0};   //<- Modified viscosity 'nu' used in the stability criterion
	for (size_t i{0}; i<mesh.nNodesTotal; ++i)
	{
		double nu = mu(i) / (rho(i)+1) * viscosityModifier;
		maxViscosityFactor = max( maxViscosityFactor, nu );
	}
	double dx_2{mesh.dx*mesh.dx}, dy_2{mesh.dy*mesh.dy}, dz_2{mesh.dz*mesh.dz};
	return fabs(params.viscStabilityLimit) / ( maxViscosityFactor * (1/dx_2 + 1/dy_2 + 1/dz_2) );
}


// Advances the conserved variables, primitive variables and transport properties from time t to t + dt, using RK4.
// Updates time step size per the stability criterion after advancing the solution.
void Solver::marchTimeStep(double& t, ulong& timeLevel)
{
	mesh.checkFinity(mesh.conservedVariables);
	mesh.applyFilter_ifAppropriate(mesh.conservedVariables.rho, mesh.conservedVariablesOld.rho, params.filterInterval, timeLevel);
	mesh.swapConservedVariables();	// Swap conserved variables to temporary storage by move-semantics.
	// Now, the previously calculated variables (current time level) are stored in conservedVariablesOld.
	// We will use the "main" conserved variable arrays for the intermediate solutions between RK4 steps,
	// and eventually the new solution at the end of this function.
	ConservedVariablesArrayGroup& k1{mesh.RK4slopes.k1}, &k2{mesh.RK4slopes.k2}, &k3{mesh.RK4slopes.k3}, &k4{mesh.RK4slopes.k4};

	computeRK4slopes(mesh.conservedVariablesOld, k1);	// Compute step 1 (k1), i.e. the slopes at time t, using Euler's method
	computeAllIntermediateSolutions(k1, dt/2);		// Compute intermediate solutions at time t + dt/2, using the slopes k1
	updatePrimitiveVariables();
	mesh.applyAllBoundaryConditions(t+dt/2, params);

	computeRK4slopes(mesh.conservedVariables, k2);  // k2: The slopes at time t + dt/2
	computeAllIntermediateSolutions(k2, dt/2);	// Compute intermediate solutions at time t + dt/2, using the slopes k2.
	updatePrimitiveVariables();
	mesh.applyAllBoundaryConditions(t+dt/2, params);

	computeRK4slopes(mesh.conservedVariables, k3);	// k3: Again, the slopes at time t + dt/2, but now using k2
	computeAllIntermediateSolutions(k3, dt);	// Compute intermediate solutions at time t + dt, using the slopes k3.
	updatePrimitiveVariables();
	mesh.applyAllBoundaryConditions(t+dt, params);

	computeRK4slopes(mesh.conservedVariables, k4);	// k4: The slopes at time t + dt
	computeRK4finalStepAllVariables();	// Use slopes k1, ..., k4 to compute solution at t+dt
	updatePrimitiveVariables();
	mesh.applyAllBoundaryConditions(t+dt, params);

	// At this stage, the conserved variables at the next time level are stored in the intermediate arrays.
	storeNormsOfChange();			// Compute norm of the change between old and new time levels, and add to history.

	t += dt;
	++timeLevel;
	updateTimeStepSize(t);
}

void Solver::computeRK4slopes(const ConservedVariablesArrayGroup& conservedVariables, ConservedVariablesArrayGroup& RK4slopes)
{
	compute_RK4_step_continuity(conservedVariables.rho_u, conservedVariables.rho_v, conservedVariables.rho_w, RK4slopes.rho  );
	compute_RK4_step_xMomentum (conservedVariables.rho_u,               									  RK4slopes.rho_u);
	compute_RK4_step_yMomentum (conservedVariables.rho_v,               									  RK4slopes.rho_v);
	compute_RK4_step_zMomentum (conservedVariables.rho_w,               									  RK4slopes.rho_w);
	compute_RK4_step_energy    (conservedVariables.rho_E,                   								  RK4slopes.rho_E);
}

// Computes an RK4 step (slope: k1, ..., k4) for the mass density. The previous state of the momentum densities, rho_u,
// rho_v and rho_w can be supplied at time t (pass the member 'rho_u' as argument 'rho_u', etc.) or as intermediate
// solutions computed in the calling function etc. Result is stored in argument 'RK4_slope'.
void Solver::compute_RK4_step_continuity(const Array3D_d& rho_u, const Array3D_d& rho_v, const Array3D_d& rho_w,
		                                    Array3D_d& RK4_slope)
{
	const Vector3_u nMeshNodes(mesh.NI, mesh.NJ, mesh.NK);
	for(size_t index1D : mesh.indexByType.fluidInterior)
	{
		Vector3_u indices3D = getIndices3D(index1D, nMeshNodes);
		size_t i = indices3D.i, j = indices3D.j, k = indices3D.k;
		RK4_slope(i,j,k) = - ( rho_u(i+1, j  , k  )-rho_u(i-1, j  , k  ) ) / (2*mesh.dx)
							    		   - ( rho_v(i  , j+1, k  )-rho_v(i  , j-1, k  ) ) / (2*mesh.dy)
										   - ( rho_w(i  , j  , k+1)-rho_w(i  , j  , k-1) ) / (2*mesh.dz);
	}
}

// Computes an RK4 step (slope: k1, ..., k4) for the x-momentum density. The previous rho_u can be
// supplied at time t (pass the member 'rho_u' as argument 'rho_u') or as an intermediate solution
// computed in the calling function etc. Result is stored in argument 'RK4_slope'.
void Solver::compute_RK4_step_xMomentum(const Array3D_d& rho_u, Array3D_d& RK4_slope)
{
	double neg_1_div_2dx   = -1. / ( 2 * mesh.dx );       // <- The parts of the flux/residual evaluation
	double neg_1_div_2dy   = -1. / ( 2 * mesh.dy );       //    that don't change from node to node are
	double neg_1_div_2dz   = -1. / ( 2 * mesh.dz );       //    pre-computed, outside the loop.
	double pos_2_div_3dxdx =  2. / ( 3 * mesh.dx * mesh.dx );
	double pos_1_div_2dydy =  1. / ( 2 * mesh.dy * mesh.dy );
	double pos_1_div_2dzdz =  1. / ( 2 * mesh.dz * mesh.dz );
	double neg_1_div_6dxdy = -1. / ( 6 * mesh.dx * mesh.dy );
	double pos_1_div_4dxdy =  1. / ( 4 * mesh.dx * mesh.dy );
	double neg_1_div_6dxdz = -1. / ( 6 * mesh.dx * mesh.dz );
	double pos_1_div_4dxdz =  1. / ( 4 * mesh.dx * mesh.dz );

	const Array3D_d& p{mesh.primitiveVariables.p}, &u{mesh.primitiveVariables.u}, &v{mesh.primitiveVariables.v}, &w{mesh.primitiveVariables.w},
		&mu{mesh.transportProperties.mu}; // For readability in math expression below.

	const Vector3_u nMeshNodes(mesh.NI, mesh.NJ, mesh.NK);

	for(size_t index1D : mesh.indexByType.fluidInterior)
	{
		Vector3_u indices3D = getIndices3D(index1D, nMeshNodes);
		size_t i = indices3D.i, j = indices3D.j, k = indices3D.k;
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

// Computes an RK4 step (slope: k1, ..., k4) for the y-momentum density. The previous rho_v can be
// supplied at time t (pass the member 'rho_v' as argument 'rho_v') or as an intermediate solution
// computed in the calling function etc. Result is stored in argument 'RK4_slope'.
void Solver::compute_RK4_step_yMomentum(const Array3D_d& rho_v, Array3D_d& RK4_slope)
{
	double neg_1_div_2dx   = -1. / ( 2 * mesh.dx );       // <- The parts of the flux/residual evaluation
	double neg_1_div_2dy   = -1. / ( 2 * mesh.dy );       //    that don't change from node to node are
	double neg_1_div_2dz   = -1. / ( 2 * mesh.dz );       //    pre-computed, outside the loop.
	double pos_1_div_2dxdx =  1. / ( 2 * mesh.dx * mesh.dx );
	double pos_2_div_3dydy =  2. / ( 3 * mesh.dy * mesh.dy );
	double pos_1_div_2dzdz =  1. / ( 2 * mesh.dz * mesh.dz );
	double neg_1_div_6dxdy = -1. / ( 6 * mesh.dx * mesh.dy );
	double pos_1_div_4dydx =  1. / ( 4 * mesh.dy * mesh.dx );
	double neg_1_div_6dzdy = -1. / ( 6 * mesh.dz * mesh.dy );
	double pos_1_div_4dydz =  1. / ( 4 * mesh.dy * mesh.dz );

	const Array3D_d& p{mesh.primitiveVariables.p}, &u{mesh.primitiveVariables.u}, &v{mesh.primitiveVariables.v}, &w{mesh.primitiveVariables.w},
		&mu{mesh.transportProperties.mu}; // For readability in math expression below.

	const Vector3_u nMeshNodes(mesh.NI, mesh.NJ, mesh.NK);

	for(size_t index1D : mesh.indexByType.fluidInterior)
	{
		Vector3_u indices3D = getIndices3D(index1D, nMeshNodes);
		size_t i = indices3D.i, j = indices3D.j, k = indices3D.k;
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

// Computes an RK4 step (slope: k1, ..., k4) for the z-momentum density. The previous rho_w can be
// supplied at time t (pass the member 'rho_w' as argument 'rho_w') or as an intermediate solution
// computed in the calling function etc. Result is stored in argument 'RK4_slope'.
void Solver::compute_RK4_step_zMomentum(const Array3D_d& rho_w, Array3D_d& RK4_slope)
{
	double neg_1_div_2dx   = -1. / ( 2 * mesh.dx );       // <- The parts of the flux/residual evaluation
	double neg_1_div_2dy   = -1. / ( 2 * mesh.dy );       //    that don't change from node to node are
	double neg_1_div_2dz   = -1. / ( 2 * mesh.dz );       //    pre-computed, outside the loop.
	double pos_1_div_2dxdx =  1. / ( 2 * mesh.dx * mesh.dx );
	double pos_1_div_2dydy =  1. / ( 2 * mesh.dy * mesh.dy );
	double pos_2_div_3dzdz =  2. / ( 3 * mesh.dz * mesh.dz );
	double pos_1_div_4dzdx =  1. / ( 4 * mesh.dz * mesh.dx );
	double neg_1_div_6dxdz = -1. / ( 6 * mesh.dx * mesh.dz );
	double pos_1_div_4dzdy =  1. / ( 4 * mesh.dz * mesh.dy );
	double neg_1_div_6dydz = -1. / ( 6 * mesh.dy * mesh.dz );

	const Array3D_d& p{mesh.primitiveVariables.p}, &u{mesh.primitiveVariables.u}, &v{mesh.primitiveVariables.v}, &w{mesh.primitiveVariables.w},
		&mu{mesh.transportProperties.mu}; // For readability in math expression below.

	const Vector3_u nMeshNodes(mesh.NI, mesh.NJ, mesh.NK);

	for(size_t index1D : mesh.indexByType.fluidInterior)
	{
		Vector3_u indices3D = getIndices3D(index1D, nMeshNodes);
		size_t i = indices3D.i, j = indices3D.j, k = indices3D.k;

		double convFluxX = neg_1_div_2dx   * ( rho_w(i+1,j,k) * u(i+1,j,k) - rho_w(i-1,j,k) * u(i-1,j,k) );
		double convFluxY = neg_1_div_2dy   * ( rho_w(i,j+1,k) * v(i,j+1,k) - rho_w(i,j-1,k) * v(i,j-1,k) );
		double convFluxZ = neg_1_div_2dz   * ( rho_w(i,j,k+1) * w(i,j,k+1) + p(i,j,k+1)
											 - rho_w(i,j,k-1) * w(i,j,k-1) - p(i,j,k-1) );

		double viscFluxXX = pos_1_div_2dxdx * ( (mu(i+1,j,k) + mu(i,j,k)) * (w(i+1,j,k) - w(i  ,j,k))
											  - (mu(i-1,j,k) + mu(i,j,k)) * (w(i  ,j,k) - w(i-1,j,k)) );
		double viscFluxYY = pos_1_div_2dydy * ( (mu(i,j+1,k) + mu(i,j,k)) * (w(i,j+1,k) - w(i,j  ,k))
											  - (mu(i,j-1,k) + mu(i,j,k)) * (w(i,j  ,k) - w(i,j-1,k)) );
		double viscFluxZZ = pos_2_div_3dzdz * ( (mu(i,j,k+1) + mu(i,j,k)) * (w(i,j,k+1) - w(i,j,k  ))
											  - (mu(i,j,k-1) + mu(i,j,k)) * (w(i,j,k  ) - w(i,j,k-1)) );
		double viscFluxZX = pos_1_div_4dzdx * ( mu(i+1,j,k) * (u(i+1,j,k+1) - u(i+1,j,k-1)) - mu(i-1,j,k) * (u(i-1,j,k+1) - u(i-1,j,k-1)) );
		double viscFluxXZ = neg_1_div_6dxdz * ( mu(i,j,k+1) * (u(i+1,j,k+1) - u(i-1,j,k+1)) - mu(i,j,k-1) * (u(i+1,j,k-1) - u(i-1,j,k-1)) );
		double viscFluxZY = pos_1_div_4dzdy * ( mu(i,j+1,k) * (v(i,j+1,k+1) - v(i,j+1,k-1)) - mu(i,j-1,k) * (v(i,j-1,k+1) - v(i,j-1,k-1)) );
		double viscFluxYZ = neg_1_div_6dydz * ( mu(i,j,k+1) * (v(i,j+1,k+1) - v(i,j-1,k+1)) - mu(i,j,k-1) * (v(i,j+1,k-1) - v(i,j-1,k-1)) );
		RK4_slope(i,j,k) = convFluxX + convFluxY + convFluxZ
						 + viscFluxXX + viscFluxYY + viscFluxZZ
						 + viscFluxZX + viscFluxXZ + viscFluxZY + viscFluxYZ ;
	}
}

// Computes an RK4 step (slope: k1, ..., k4) for the total specific energy. The previous energy can
// be supplied at time t (pass the member 'E' as argument 'E') or as an intermediate solution
// computed in the calling function etc. Result is stored in argument 'RK4_slope'.
void Solver::compute_RK4_step_energy(const Array3D_d& E, Array3D_d& RK4_slope)
{
	double neg_1_div_2dx   = -1. / ( 2 * mesh.dx );       // <- The parts of the flux/residual evaluation
	double neg_1_div_2dy   = -1. / ( 2 * mesh.dy );       //    that don't change from node to node are
	double neg_1_div_2dz   = -1. / ( 2 * mesh.dz );       //    pre-computed, outside the loop.
	double pos_2_div_3dxdx =  2. / ( 3 * mesh.dx * mesh.dx );
	double pos_1_div_2dxdx =  1. / ( 2 * mesh.dx * mesh.dx );
	double neg_1_div_6dydx = -1. / ( 6 * mesh.dy * mesh.dx );
	double neg_1_div_6dzdx = -1. / ( 6 * mesh.dz * mesh.dx );
	double pos_1_div_4dydx =  1. / ( 4 * mesh.dy * mesh.dx );
	double pos_1_div_4dzdx =  1. / ( 4 * mesh.dz * mesh.dx );
	double pos_1_div_2dydy =  1. / ( 2 * mesh.dy * mesh.dy );
	double pos_2_div_3dydy =  2. / ( 3 * mesh.dy * mesh.dy );
	double neg_1_div_6dzdy = -1. / ( 6 * mesh.dz * mesh.dy );
	double pos_1_div_4dzdy =  1. / ( 4 * mesh.dz * mesh.dy );
	double pos_1_div_2dzdz =  1. / ( 2 * mesh.dz * mesh.dz );
	double pos_2_div_3dzdz =  2. / ( 3 * mesh.dz * mesh.dz );
	double rho_H_0 = 1. / ( params.Gamma - 1 );     // <- Dimensionless reference enthalpy

	const Array3D_d& p{mesh.primitiveVariables.p}, &u{mesh.primitiveVariables.u}, &v{mesh.primitiveVariables.v}, &w{mesh.primitiveVariables.w},
		&T{mesh.primitiveVariables.T}, &mu{mesh.transportProperties.mu}, &kappa{mesh.transportProperties.kappa}; // For readability in math expression below.

	const Vector3_u nMeshNodes(mesh.NI, mesh.NJ, mesh.NK);

	for(size_t index1D : mesh.indexByType.fluidInterior)
	{
		Vector3_u indices3D = getIndices3D(index1D, nMeshNodes);
		size_t i = indices3D.i, j = indices3D.j, k = indices3D.k;
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

void Solver::computeAllIntermediateSolutions(const ConservedVariablesArrayGroup& RK4slopes, double timeIncrement)
{
	computeIntermediateSolution(mesh.conservedVariablesOld.rho,   RK4slopes.rho,   mesh.conservedVariables.rho,   timeIncrement);
	computeIntermediateSolution(mesh.conservedVariablesOld.rho_u, RK4slopes.rho_u, mesh.conservedVariables.rho_u, timeIncrement);
	computeIntermediateSolution(mesh.conservedVariablesOld.rho_v, RK4slopes.rho_v, mesh.conservedVariables.rho_v, timeIncrement);
	computeIntermediateSolution(mesh.conservedVariablesOld.rho_w, RK4slopes.rho_w, mesh.conservedVariables.rho_w, timeIncrement);
	computeIntermediateSolution(mesh.conservedVariablesOld.rho_E, RK4slopes.rho_E, mesh.conservedVariables.rho_E, timeIncrement);
}

// Evaluates an intermediate solution of the given conserved variable, using the given RK4 step (slope)
// and the appropriate time increment (dt/2 if k1 or k2 is given, dt if k3 is given).
// Result is stored in 'intermSolution'. Only writes interior nodes.
void Solver::computeIntermediateSolution(const Array3D_d& conservedVar, const Array3D_d& RK4_slope,
		                                       Array3D_d& intermSolution, double timeIncrement)
{
	for(size_t index1D : mesh.indexByType.fluidInterior)
		intermSolution(index1D) = conservedVar(index1D) + timeIncrement * RK4_slope(index1D);
}

// Updates the primitive variables and transport properties, using the specified conserved variables.
// The conserved variables can be given at time t (just pass the data-members rho, rho_u, ..., E),
// or as intermediate solutions computed in the calling function etc. Only writes interior nodes!
void Solver::updatePrimitiveVariables()
{
	double gammaMinusOne = params.Gamma - 1;
	double sutherlands_Sc = params.T_0 / params.sutherlands_C2;
	double ScPlusOne = 1 + sutherlands_Sc;
	double prandtlFactor = 1 / ( gammaMinusOne * params.Pr );
	const Array3D_d& rho_u{mesh.conservedVariables.rho_u},
					&rho_v{mesh.conservedVariables.rho_v},
					&rho_w{mesh.conservedVariables.rho_w},
					&rho  {mesh.conservedVariables.rho},
					&rho_E{mesh.conservedVariables.rho_E};
	Array3D_d& u{mesh.primitiveVariables.u},
			  &v{mesh.primitiveVariables.v},
			  &w{mesh.primitiveVariables.w},
			  &p{mesh.primitiveVariables.p},
			  &T{mesh.primitiveVariables.T},
			  &mu{mesh.transportProperties.mu},
			  &kappa{mesh.transportProperties.kappa};

	for(size_t i : mesh.indexByType.fluidInterior)
		u(i) = rho_u(i) / ( rho(i) + 1 );
	for(size_t i : mesh.indexByType.fluidInterior)
		v(i) = rho_v(i) / ( rho(i) + 1 );
	for(size_t i : mesh.indexByType.fluidInterior)
		w(i) = rho_w(i) / ( rho(i) + 1 );
	for(size_t i : mesh.indexByType.fluidInterior)
		p(i) = gammaMinusOne * ( rho_E(i) - (rho(i)+1)/2 * (u(i)*u(i) + v(i)*v(i) + w(i)*w(i)) );
	for(size_t i : mesh.indexByType.fluidInterior)
		T(i) = ( params.Gamma * p(i) - rho(i) ) / ( 1+rho(i) );
	for(size_t i : mesh.indexByType.fluidInterior)
		mu(i) = pow( 1+T(i), 1.5 ) * ScPlusOne / ( params.Re*( T(i) + ScPlusOne ) );
	for(size_t i : mesh.indexByType.fluidInterior)
		kappa(i) = mu(i) * prandtlFactor;
}

void Solver::computeRK4finalStepAllVariables()
{
	const ConservedVariablesArrayGroup& k1{mesh.RK4slopes.k1}, &k2{mesh.RK4slopes.k2}, &k3{mesh.RK4slopes.k3}, &k4{mesh.RK4slopes.k4};
	compute_RK4_final_step(k1.rho,   k2.rho,   k3.rho,   k4.rho,   mesh.conservedVariablesOld.rho  , mesh.conservedVariables.rho  );  // Use a weighted average of the four slopes
	compute_RK4_final_step(k1.rho_u, k2.rho_u, k3.rho_u, k4.rho_u, mesh.conservedVariablesOld.rho_u, mesh.conservedVariables.rho_u);  // to advance the solution from time t, to
	compute_RK4_final_step(k1.rho_v, k2.rho_v, k3.rho_v, k4.rho_v, mesh.conservedVariablesOld.rho_v, mesh.conservedVariables.rho_v);  // t + dt. The solutions at the new time level
	compute_RK4_final_step(k1.rho_w, k2.rho_w, k3.rho_w, k4.rho_w, mesh.conservedVariablesOld.rho_w, mesh.conservedVariables.rho_w);  // are stored in the intermediate arrays, so the
	compute_RK4_final_step(k1.rho_E, k2.rho_E, k3.rho_E, k4.rho_E, mesh.conservedVariablesOld.rho_E, mesh.conservedVariables.rho_E);  // norm of the change can be computed.
}

// Uses all the intermediate solutions for a variable, k1,2,3,4, to advance the variable from t to t+dt
void Solver::compute_RK4_final_step(const Array3D_d& k1, const Array3D_d& k2,
									const Array3D_d& k3, const Array3D_d& k4,
									const Array3D_d& conservedVar_old, Array3D_d& conservedVar_new)
{
	double dtFactor{dt / 6};
	for(size_t i : mesh.indexByType.fluidInterior)
	{
		double residualValue = dtFactor*( k1(i) + 2*k2(i) + 2*k3(i) + k4(i) );
		conservedVar_new(i) = conservedVar_old(i) + residualValue;
	}
}

void Solver::storeNormsOfChange()
{
	if(params.saveConvergenceHistory != ConvHistoryEnum::none)
		normHistory.rho.push_back( mesh.getNormOfChange(mesh.conservedVariablesOld.rho,
														mesh.conservedVariables.rho) );

	if(params.saveConvergenceHistory == ConvHistoryEnum::all)
	{
		normHistory.rho_u.push_back( mesh.getNormOfChange(mesh.conservedVariablesOld.rho_u,
														  mesh.conservedVariables   .rho_u) );

		normHistory.rho_v.push_back( mesh.getNormOfChange(mesh.conservedVariablesOld.rho_v,
														  mesh.conservedVariables   .rho_v) );

		normHistory.rho_w.push_back( mesh.getNormOfChange(mesh.conservedVariablesOld.rho_w,
														  mesh.conservedVariables   .rho_w) );

		normHistory.rho_E.push_back( mesh.getNormOfChange(mesh.conservedVariablesOld.rho_E,
														  mesh.conservedVariables   .rho_E) );
	}
}






