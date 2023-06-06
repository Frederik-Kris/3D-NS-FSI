/*
 * Solver.cpp
 *
 *  Created on: Apr 28, 2021
 *      Author: frederik
 */


#include "Solver.h"
#include "ConstantLiteralExpressions.h"
#include "SolutionImporter.h"

// Constructor. Takes ConfigSettings 'params', with necessary settings from config file
// Initializes mesh with 3D arrays (but not their elements)
Solver::Solver(const ConfigSettings& params) :
params(params),
mesh(params),
dt{0}
{
}

// Preparing solver before starting to march solution.
// Setup boundaries, apply initial condition (IC) and set the time-step size based on the IC.
void Solver::initialize(long unsigned timeLevel)
{
	mesh.categorizeNodes(params);
	applyInitialConditions();
	if(params.continueSimulation)
	{
		SolutionImporter reader(params);
		reader.importNormHistories(normHistory, timeLevel);
		reader.importLiftDrag(liftHistory, dragHistory, timeLevel);
	}
	updateTimeStepSize(0);
}

// Applies stagnation initial condition (IC) by setting the values of flow variables at every node
void Solver::applyStagnation_IC()
{
	cout << "Applying stagnation initial condition... ";
	double u_IC{0}, v_IC{0}, w_IC{0}, p_IC{0}, T_IC{0};                  	// ← Decide primitive variables
	PrimitiveVariablesScalars  decidedPrimitiveVariables(u_IC, v_IC, w_IC, p_IC, T_IC); // and derive the others
	ConservedVariablesScalars  derivedConservedVariables  = deriveConservedVariables(decidedPrimitiveVariables, params);
	TransportPropertiesScalars derivedTransportProperties = deriveTransportProperties(decidedPrimitiveVariables, params);
	for(SubMesh& subMesh : mesh.subMeshes)
	{
		subMesh.primitiveVariables.u.setAll(decidedPrimitiveVariables.u);
		subMesh.primitiveVariables.v.setAll(decidedPrimitiveVariables.v);
		subMesh.primitiveVariables.w.setAll(decidedPrimitiveVariables.w);
		subMesh.primitiveVariables.p.setAll(decidedPrimitiveVariables.p);
		subMesh.primitiveVariables.T.setAll(decidedPrimitiveVariables.T);
		subMesh.conservedVariables.rho  .setAll(derivedConservedVariables.rho);
		subMesh.conservedVariables.rho_u.setAll(derivedConservedVariables.rho_u);
		subMesh.conservedVariables.rho_v.setAll(derivedConservedVariables.rho_v);
		subMesh.conservedVariables.rho_w.setAll(derivedConservedVariables.rho_w);
		subMesh.conservedVariables.rho_E.setAll(derivedConservedVariables.rho_E);
		subMesh.transportProperties.mu   .setAll(derivedTransportProperties.mu);
		subMesh.transportProperties.kappa.setAll(derivedTransportProperties.kappa);
	}
	cout << "Done" << endl << endl;
}

// Applies uniform flow initial condition (IC) by setting the values of flow variables at every node
void Solver::applyUniformFlow_IC()
{
	cout << "Applying oblique uniform flow initial condition... ";
	double velocity = params.M_0 / sqrt(3);
	double u_IC{velocity}, v_IC{velocity}, w_IC{velocity}, p_IC{0}, T_IC{0}; // ← Decide primitive variables
	PrimitiveVariablesScalars  decidedPrimitiveVariables(u_IC, v_IC, w_IC, p_IC, T_IC);// and derive the others
	ConservedVariablesScalars  derivedConservedVariables  = deriveConservedVariables(decidedPrimitiveVariables, params);
	TransportPropertiesScalars derivedTransportProperties = deriveTransportProperties(decidedPrimitiveVariables, params);
	for(SubMesh& subMesh : mesh.subMeshes)
	{
		subMesh.primitiveVariables.u.setAll(decidedPrimitiveVariables.u);
		subMesh.primitiveVariables.v.setAll(decidedPrimitiveVariables.v);
		subMesh.primitiveVariables.w.setAll(decidedPrimitiveVariables.w);
		subMesh.primitiveVariables.p.setAll(decidedPrimitiveVariables.p);
		subMesh.primitiveVariables.T.setAll(decidedPrimitiveVariables.T);
		subMesh.conservedVariables.rho  .setAll(derivedConservedVariables.rho);
		subMesh.conservedVariables.rho_u.setAll(derivedConservedVariables.rho_u);
		subMesh.conservedVariables.rho_v.setAll(derivedConservedVariables.rho_v);
		subMesh.conservedVariables.rho_w.setAll(derivedConservedVariables.rho_w);
		subMesh.conservedVariables.rho_E.setAll(derivedConservedVariables.rho_E);
		subMesh.transportProperties.mu   .setAll(derivedTransportProperties.mu);
		subMesh.transportProperties.kappa.setAll(derivedTransportProperties.kappa);
	}
	cout << "Done" << endl << endl;
}

// Apply initial conditions (IC) for the solution. In case of continuation from previous simulation, the last
// output is loaded and used as IC.
void Solver::applyInitialConditions()
{
	if(!params.continueSimulation)
		applyStagnation_IC();
/*	else
	{
		SolutionImporter reader(params);
		reader.loadVtkFile();
		if( !reader.allFlowVarsDerivable() )
		{
			std::cerr << "Not possible to reconstruct all flow variables using the output in the .vtk file." << endl;
			applyStagnation_IC();
			return;
		}
		int index1D = 0; // NB! Note loop order. Don't use index1D here to access arrays in mesh.
		for(int k{0}; k<mesh.NK; ++k)
			for(int j{0}; j<mesh.NJ; ++j)
				for(int i{0}; i<mesh.NI; ++i)
				{
					ConservedVariablesScalars conservedVars;
					PrimitiveVariablesScalars primitiveVars;
					TransportPropertiesScalars transportProps;
					reader.computeAllFlowVars(index1D, conservedVars, primitiveVars, transportProps);
					setFlowVariablesAtNode( Vector3_i(i,j,k), conservedVars, primitiveVars, transportProps, mesh.flowVariableReferences);
					++index1D;
				}
	}*/
}

// Set the time step size (dt) as large as possible, within the stability criterion, without surpassing the end-time.
// Computes two time step sizes, using the inviscid CFL condition, and the viscous von Neumann condition
// and assigns the smallest one to dt (strictest criterion).
void Solver::updateTimeStepSize(double t)
{
	double dtConvective = getInviscidTimeStepLimit();
	double dtViscous    = getViscousTimeStepLimit();
	// Choose the strictest criterion. If that will take us past the end-time, then adapt dt to hit t_end exactly:
	dt = min( dtConvective, dtViscous );
	if ( params.stopCriterion == StopCriterionEnum::end_time && t + dt > params.t_end )
		dt = params.t_end - t;
}

// Find the inviscid time step limit (CFL condition):
double Solver::getInviscidTimeStepLimit()
{
	double maxSpectralRadiusX{0}, maxSpectralRadiusY{0}, maxSpectralRadiusZ{0};
	// Loop through all mesh nodes, to find largest spectral radii:
	for(SubMesh& subMesh : mesh.subMeshes)
	{
		const Array3D_d& rho{subMesh.conservedVariables.rho};
		const Array3D_d& p  {subMesh.primitiveVariables.p};
		const Array3D_d& u  {subMesh.primitiveVariables.u};
		const Array3D_d& v  {subMesh.primitiveVariables.v};
		const Array3D_d& w  {subMesh.primitiveVariables.w};
		for (int i{0}; i<subMesh.nNodesTotal; ++i)
		{
			double c_i = sqrt( (1 + params.Gamma * p(i)) / (1 + rho(i)) );    // ← Speed of sound at node i
			maxSpectralRadiusX = max( maxSpectralRadiusX, c_i + fabs(u(i)) );
			maxSpectralRadiusY = max( maxSpectralRadiusY, c_i + fabs(v(i)) );
			maxSpectralRadiusZ = max( maxSpectralRadiusZ, c_i + fabs(w(i)) );
		}
	}
	return fabs(params.convStabilityLimit) / ( maxSpectralRadiusX / mesh.smallestGridSpacings.x
			                                 + maxSpectralRadiusY / mesh.smallestGridSpacings.y
											 + maxSpectralRadiusZ / mesh.smallestGridSpacings.z );
}

// Find the viscous time step size limit (von Neumann condition):
double Solver::getViscousTimeStepLimit()
{
	double viscosityModifier = max( 4./3, params.Gamma / params.Pr );
	double maxViscosityFactor{0};   //← Modified viscosity 'nu' used in the stability criterion
	for(SubMesh& subMesh : mesh.subMeshes)
	{
		const Array3D_d& mu {subMesh.transportProperties.mu};
		const Array3D_d& rho{subMesh.conservedVariables.rho};
		for (int i{0}; i<subMesh.nNodesTotal; ++i)
		{
			double nu = mu(i) / (rho(i)+1) * viscosityModifier;
			maxViscosityFactor = max( maxViscosityFactor, nu );
		}
	}
	double dx = mesh.smallestGridSpacings.x;
	double dy = mesh.smallestGridSpacings.y;
	double dz = mesh.smallestGridSpacings.z;
	return fabs(params.viscStabilityLimit) / ( maxViscosityFactor * (1/(dx*dx) + 1/(dy*dy) + 1/(dz*dz) ) );
}

// Advances the conserved variables, primitive variables and transport properties from time t to t + dt, using RK4.
// Updates time step size per the stability criterion after advancing the solution.
void Solver::marchTimeStep(double& t, 			// ← IN-/OUTPUT, time
						   long unsigned& timeLevel)	// ← OUTPUT
{
	mesh.applyFilter_ifAppropriate(params, timeLevel, t);
	mesh.swapConservedVariables();	// Swap conserved variables and "old" arrays by move-semantics.
	// Now, the previously calculated variables (current time level) are stored in conservedVariablesOld.
	// We will use the "main" conserved variable arrays (mesh.conservedVariables) for the intermediate solutions
	// between RK4 steps, and eventually for the new solution at the end of this function.

	computeRK4slopes(ChooseSolutionStage::old, ChooseSlopeStage::k1);	// Compute step 1 (k1), i.e. the slopes at time t, using Euler's method
	computeAllIntermediateSolutions(ChooseSlopeStage::k1, dt/2);		// Compute intermediate solutions at time t + dt/2, using the slopes k1
	updatePrimitiveVariables();
	mesh.applyAllBoundaryConditions(t+dt/2, params);

	computeRK4slopes(ChooseSolutionStage::intermediate, ChooseSlopeStage::k2);  // k2: The slopes at time t + dt/2
	computeAllIntermediateSolutions(ChooseSlopeStage::k2, dt/2);	// Compute intermediate solutions at time t + dt/2, using the slopes k2.
	updatePrimitiveVariables();
	mesh.applyAllBoundaryConditions(t+dt/2, params);

	computeRK4slopes(ChooseSolutionStage::intermediate, ChooseSlopeStage::k3);	// k3: Again, the slopes at time t + dt/2, but now using k2
	computeAllIntermediateSolutions(ChooseSlopeStage::k3, dt);	// Compute intermediate solutions at time t + dt, using the slopes k3.
	updatePrimitiveVariables();
	mesh.applyAllBoundaryConditions(t+dt, params);

	computeRK4slopes(ChooseSolutionStage::intermediate, ChooseSlopeStage::k4);	// k4: The slopes at time t + dt
	computeRK4finalStepAllVariables();	// Use slopes k1, ..., k4 to compute solution at t+dt
	updatePrimitiveVariables();
	mesh.applyAllBoundaryConditions(t+dt, params);

	SubMesh& subMeshWithIB = mesh.subMeshes(2,2,0);
	if( subMeshWithIB.getImmersedBoundaries().empty() )
		throw std::runtime_error("'subMmeshWithIB' does not have any immersed boundary.");
	IntegralProperties integralProperties = subMeshWithIB.getImmersedBoundaries().front()->getIntegralProperties(params);
	liftHistory.push_back(integralProperties.lift);
	dragHistory.push_back(integralProperties.drag);
	separationAngles.clear();
	separationAngles = integralProperties.separationAngles;

	// At this stage, the conserved variables at the next time level are stored in mesh.conservedVariables .
	storeNormsOfChange();			// Compute norm of the change between old and new time levels, and add to history.

	t += dt;
	++timeLevel;
	updateTimeStepSize(t);
}

// Compute RK4 slopes for all the conserved variables, using specified solution.
void Solver::computeRK4slopes(ChooseSolutionStage solutionStage, 	// ← Current(old) or intermediate solution
							  ChooseSlopeStage slopeStage)			// ← Slopes to overwrite
{
	for(SubMesh& subMesh : mesh.subMeshes)
	{
		ConservedVariablesArrayGroup* conservedVars = nullptr;
		ConservedVariablesArrayGroup* RK4Slopes 	= nullptr;
		if	   (solutionStage == ChooseSolutionStage::old)
			conservedVars = &subMesh.conservedVariablesOld;
		else if(solutionStage == ChooseSolutionStage::intermediate)
			conservedVars = &subMesh.conservedVariables;
		else throw std::logic_error("Unexpected enum value for 'ChooseSolutionStage' when computing RK4 slopes");
		switch(slopeStage)
		{
		case ChooseSlopeStage::k1:
			RK4Slopes = &subMesh.RK4slopes.k1;	break;
		case ChooseSlopeStage::k2:
			RK4Slopes = &subMesh.RK4slopes.k2;	break;
		case ChooseSlopeStage::k3:
			RK4Slopes = &subMesh.RK4slopes.k3;	break;
		case ChooseSlopeStage::k4:
			RK4Slopes = &subMesh.RK4slopes.k4;	break;
		default:
			throw std::logic_error("Unexpected enum value for 'ChooseSlopeStage' when computing RK4 slopes");
		}
		compute_RK4_step_continuity(conservedVars->rho_u, conservedVars->rho_v, conservedVars->rho_w, RK4Slopes->rho  );
		compute_RK4_step_xMomentum (conservedVars->rho_u,               							  RK4Slopes->rho_u);
		compute_RK4_step_yMomentum (conservedVars->rho_v,               							  RK4Slopes->rho_v);
		compute_RK4_step_zMomentum (conservedVars->rho_w,               							  RK4Slopes->rho_w);
		compute_RK4_step_energy    (conservedVars->rho_E,                   						  RK4Slopes->rho_E);
	}
}

// Computes one RK4 step (slope: k1, ..., k4) for the mass density. Result is stored in argument 'RK4_slope'.
void Solver::compute_RK4_step_continuity(const Array3D_d& rho_u, // ← Momentum density
										 const Array3D_d& rho_v, // ← Momentum density
										 const Array3D_d& rho_w, // ← Momentum density
										 Array3D_d& RK4_slope)	 // ← OUTPUT, slope to overwrite
{
	for(SubMesh& subMesh : mesh.subMeshes)
	{
		double dx = subMesh.gridSpacings.x;
		double dy = subMesh.gridSpacings.y;
		double dz = subMesh.gridSpacings.z;
		for(int index1D : subMesh.indexByType.fluidInterior) // Loop through interior fluid nodes
		{
			Vector3_i indices3D = getIndices3D(index1D, subMesh.arrayLimits);
			int i = indices3D.i, j = indices3D.j, k = indices3D.k;
			double massFluxX = - ( rho_u(i+1, j  , k  )-rho_u(i-1, j  , k  ) ) / (2*dx);
			double massFluxY = - ( rho_v(i  , j+1, k  )-rho_v(i  , j-1, k  ) ) / (2*dy);
			double massFluxZ = - ( rho_w(i  , j  , k+1)-rho_w(i  , j  , k-1) ) / (2*dz);
			RK4_slope(i,j,k) = massFluxX + massFluxY + massFluxZ;
		}
	}
}

// Computes an RK4 step (slope: k1, ..., k4) for the x-momentum density. Result is stored in argument 'RK4_slope'.
void Solver::compute_RK4_step_xMomentum(const Array3D_d& rho_u, // ← Momentum density
										Array3D_d& RK4_slope)	// ← OUTPUT, slope to overwrite
{
	for(SubMesh& subMesh : mesh.subMeshes)
	{
		double dx = subMesh.gridSpacings.x;
		double dy = subMesh.gridSpacings.y;
		double dz = subMesh.gridSpacings.z;
		double neg_1_div_2dx   = -1. / ( 2 * dx );       // ← The parts of the flux/residual evaluation
		double neg_1_div_2dy   = -1. / ( 2 * dy );       //    that don't change from node to node.
		double neg_1_div_2dz   = -1. / ( 2 * dz );
		double pos_2_div_3dxdx =  2. / ( 3 * dx * dx );
		double pos_1_div_2dydy =  1. / ( 2 * dy * dy );
		double pos_1_div_2dzdz =  1. / ( 2 * dz * dz );
		double neg_1_div_6dydx = -1. / ( 6 * dy * dx );
		double pos_1_div_4dxdy =  1. / ( 4 * dx * dy );
		double neg_1_div_6dzdx = -1. / ( 6 * dz * dx );
		double pos_1_div_4dxdz =  1. / ( 4 * dx * dz );

		const Array3D_d& p {subMesh.primitiveVariables.p}; // For readability in math expression below.
		const Array3D_d& u {subMesh.primitiveVariables.u};
		const Array3D_d& v {subMesh.primitiveVariables.v};
		const Array3D_d& w {subMesh.primitiveVariables.w};
		const Array3D_d& mu{subMesh.transportProperties.mu};

		for(int index1D : subMesh.indexByType.fluidInterior) // Loop through interior fluid nodes
		{
			Vector3_i indices3D = getIndices3D(index1D, subMesh.arrayLimits);
			int i = indices3D.i, j = indices3D.j, k = indices3D.k;
			double convFluxX = neg_1_div_2dx   * ( rho_u(i+1,j,k) * u(i+1,j,k) + p(i+1,j,k) - rho_u(i-1,j,k) * u(i-1,j,k) - p(i-1,j,k) );
			double convFluxY = neg_1_div_2dy   * ( rho_u(i,j+1,k) * v(i,j+1,k) - rho_u(i,j-1,k) * v(i,j-1,k) );
			double convFluxZ = neg_1_div_2dz   * ( rho_u(i,j,k+1) * w(i,j,k+1) - rho_u(i,j,k-1) * w(i,j,k-1) );

			double viscFluxXX = pos_2_div_3dxdx * ( (mu(i+1,j,k) + mu(i,j,k)) * (u(i+1,j,k) - u(i  ,j,k))
												 - (mu(i-1,j,k) + mu(i,j,k)) * (u(i  ,j,k) - u(i-1,j,k)) );
			double viscFluxYY = pos_1_div_2dydy * ( (mu(i,j+1,k) + mu(i,j,k)) * (u(i,j+1,k) - u(i,j  ,k))
												 - (mu(i,j-1,k) + mu(i,j,k)) * (u(i,j  ,k) - u(i,j-1,k)) );
			double viscFluxZZ = pos_1_div_2dzdz * ( (mu(i,j,k+1) + mu(i,j,k)) * (u(i,j,k+1) - u(i,j,k  ))
												 - (mu(i,j,k-1) + mu(i,j,k)) * (u(i,j,k  ) - u(i,j,k-1)) );
			double viscFluxXY = pos_1_div_4dxdy * ( mu(i,j+1,k) * (v(i+1,j+1,k) - v(i-1,j+1,k)) - mu(i,j-1,k) * (v(i+1,j-1,k) - v(i-1,j-1,k)) );
			double viscFluxYX = neg_1_div_6dydx * ( mu(i+1,j,k) * (v(i+1,j+1,k) - v(i+1,j-1,k)) - mu(i-1,j,k) * (v(i-1,j+1,k) - v(i-1,j-1,k)) );
			double viscFluxXZ = pos_1_div_4dxdz * ( mu(i,j,k+1) * (w(i+1,j,k+1) - w(i-1,j,k+1)) - mu(i,j,k-1) * (w(i+1,j,k-1) - w(i-1,j,k-1)) );
			double viscFluxZX = neg_1_div_6dzdx * ( mu(i+1,j,k) * (w(i+1,j,k+1) - w(i+1,j,k-1)) - mu(i-1,j,k) * (w(i-1,j,k+1) - w(i-1,j,k-1)) );
			RK4_slope(i,j,k) = convFluxX + convFluxY + convFluxZ
							 + viscFluxXX + viscFluxYY + viscFluxZZ
							 + viscFluxXY + viscFluxYX + viscFluxXZ + viscFluxZX ;
		}
	}
}

// Computes an RK4 step (slope: k1, ..., k4) for the y-momentum density. Result is stored in argument 'RK4_slope'.
void Solver::compute_RK4_step_yMomentum(const Array3D_d& rho_v, // ← Momentum density
										Array3D_d& RK4_slope)	// ← OUTPUT, slope to overwrite
{
	for(SubMesh& subMesh : mesh.subMeshes)
	{
		double dx = subMesh.gridSpacings.x;
		double dy = subMesh.gridSpacings.y;
		double dz = subMesh.gridSpacings.z;
		double neg_1_div_2dx   = -1. / ( 2 * dx );       // ← The parts of the flux/residual evaluation
		double neg_1_div_2dy   = -1. / ( 2 * dy );       //    that don't change from node to node
		double neg_1_div_2dz   = -1. / ( 2 * dz );
		double pos_1_div_2dxdx =  1. / ( 2 * dx * dx );
		double pos_2_div_3dydy =  2. / ( 3 * dy * dy );
		double pos_1_div_2dzdz =  1. / ( 2 * dz * dz );
		double neg_1_div_6dxdy = -1. / ( 6 * dx * dy );
		double pos_1_div_4dydx =  1. / ( 4 * dy * dx );
		double neg_1_div_6dzdy = -1. / ( 6 * dz * dy );
		double pos_1_div_4dydz =  1. / ( 4 * dy * dz );

		const Array3D_d& p {subMesh.primitiveVariables.p}; // For readability in math expression below.
		const Array3D_d& u {subMesh.primitiveVariables.u};
		const Array3D_d& v {subMesh.primitiveVariables.v};
		const Array3D_d& w {subMesh.primitiveVariables.w};
		const Array3D_d& mu{subMesh.transportProperties.mu};

		for(int index1D : subMesh.indexByType.fluidInterior) // Loop through interior fluid nodes
		{
			Vector3_i indices3D = getIndices3D(index1D, subMesh.arrayLimits);
			int i = indices3D.i, j = indices3D.j, k = indices3D.k;
			double convFluxX = neg_1_div_2dx   * ( rho_v(i+1,j,k) * u(i+1,j,k) - rho_v(i-1,j,k) * u(i-1,j,k) );
			double convFluxY = neg_1_div_2dy   * ( rho_v(i,j+1,k) * v(i,j+1,k) + p(i,j+1,k) - rho_v(i,j-1,k) * v(i,j-1,k) - p(i,j-1,k) );
			double convFluxZ = neg_1_div_2dz   * ( rho_v(i,j,k+1) * w(i,j,k+1) - rho_v(i,j,k-1) * w(i,j,k-1) );

			double viscFluxXX = pos_1_div_2dxdx * ( (mu(i+1,j,k) + mu(i,j,k)) * (v(i+1,j,k) - v(i  ,j,k))
												 - (mu(i-1,j,k) + mu(i,j,k)) * (v(i  ,j,k) - v(i-1,j,k)) );
			double viscFluxYY = pos_2_div_3dydy * ( (mu(i,j+1,k) + mu(i,j,k)) * (v(i,j+1,k) - v(i,j  ,k))
												 - (mu(i,j-1,k) + mu(i,j,k)) * (v(i,j  ,k) - v(i,j-1,k)) );
			double viscFluxZZ = pos_1_div_2dzdz * ( (mu(i,j,k+1) + mu(i,j,k)) * (v(i,j,k+1) - v(i,j,k  ))
												 - (mu(i,j,k-1) + mu(i,j,k)) * (v(i,j,k  ) - v(i,j,k-1)) );
			double viscFluxYX = pos_1_div_4dydx * ( mu(i+1,j,k) * (u(i+1,j+1,k) - u(i+1,j-1,k)) - mu(i-1,j,k) * (u(i-1,j+1,k) - u(i-1,j-1,k)) );
			double viscFluxXY = neg_1_div_6dxdy * ( mu(i,j+1,k) * (u(i+1,j+1,k) - u(i-1,j+1,k)) - mu(i,j-1,k) * (u(i+1,j-1,k) - u(i-1,j-1,k)) );
			double viscFluxYZ = pos_1_div_4dydz * ( mu(i,j,k+1) * (w(i,j+1,k+1) - w(i,j-1,k+1)) - mu(i,j,k-1) * (w(i,j+1,k-1) - w(i,j-1,k-1)) ) ;
			double viscFluxZY = neg_1_div_6dzdy * ( mu(i,j+1,k) * (w(i,j+1,k+1) - w(i,j+1,k-1)) - mu(i,j-1,k) * (w(i,j-1,k+1) - w(i,j-1,k-1)) );
			RK4_slope(i,j,k) = convFluxX + convFluxY + convFluxZ
							 + viscFluxXX + viscFluxYY + viscFluxZZ
							 + viscFluxYX + viscFluxXY + viscFluxYZ + viscFluxZY ;
		}
	}
}

// Computes an RK4 step (slope: k1, ..., k4) for the z-momentum density. Result is stored in argument 'RK4_slope'.
void Solver::compute_RK4_step_zMomentum(const Array3D_d& rho_w, // ← Momentum density
										Array3D_d& RK4_slope)	// ← OUTPUT, slope to overwrite
{
	for(SubMesh& subMesh : mesh.subMeshes)
	{
		double dx = subMesh.gridSpacings.x;
		double dy = subMesh.gridSpacings.y;
		double dz = subMesh.gridSpacings.z;
		double neg_1_div_2dx   = -1. / ( 2 * dx );       // ← The parts of the flux/residual evaluation
		double neg_1_div_2dy   = -1. / ( 2 * dy );       //    that don't change from node to node
		double neg_1_div_2dz   = -1. / ( 2 * dz );
		double pos_1_div_2dxdx =  1. / ( 2 * dx * dx );
		double pos_1_div_2dydy =  1. / ( 2 * dy * dy );
		double pos_2_div_3dzdz =  2. / ( 3 * dz * dz );
		double pos_1_div_4dzdx =  1. / ( 4 * dz * dx );
		double neg_1_div_6dxdz = -1. / ( 6 * dx * dz );
		double pos_1_div_4dzdy =  1. / ( 4 * dz * dy );
		double neg_1_div_6dydz = -1. / ( 6 * dy * dz );

		const Array3D_d& p {subMesh.primitiveVariables.p}; // For readability in math expression below.
		const Array3D_d& u {subMesh.primitiveVariables.u};
		const Array3D_d& v {subMesh.primitiveVariables.v};
		const Array3D_d& w {subMesh.primitiveVariables.w};
		const Array3D_d& mu{subMesh.transportProperties.mu};

		for(int index1D : subMesh.indexByType.fluidInterior) // Loop through interior fluid nodes
		{
			Vector3_i indices3D = getIndices3D(index1D, subMesh.arrayLimits);
			int i = indices3D.i, j = indices3D.j, k = indices3D.k;

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
}

// Computes an RK4 step (slope: k1, ..., k4) for the total specific energy. Result is stored in argument 'RK4_slope'.
void Solver::compute_RK4_step_energy(const Array3D_d& E,	// ← Total specific energy per volume
									 Array3D_d& RK4_slope)	// ← OUTPUT, slope to overwrite
{
	for(SubMesh& subMesh : mesh.subMeshes)
	{
		double dx = subMesh.gridSpacings.x;
		double dy = subMesh.gridSpacings.y;
		double dz = subMesh.gridSpacings.z;
		double neg_1_div_2dx   = -1. / ( 2 * dx );       // ← The parts of the flux/residual evaluation
		double neg_1_div_2dy   = -1. / ( 2 * dy );       //    that don't change from node to node
		double neg_1_div_2dz   = -1. / ( 2 * dz );
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
		double rho_H_0 = 1. / ( params.Gamma - 1 );     // ← Dimensionless reference enthalpy

		const Array3D_d& p{subMesh.primitiveVariables.p}; // For readability in math expression below.
		const Array3D_d& u{subMesh.primitiveVariables.u};
		const Array3D_d& v{subMesh.primitiveVariables.v};
		const Array3D_d& w{subMesh.primitiveVariables.w};
		const Array3D_d& T{subMesh.primitiveVariables.T};
		const Array3D_d& mu   {subMesh.transportProperties.mu};
		const Array3D_d& kappa{subMesh.transportProperties.kappa};

		for(int index1D : subMesh.indexByType.fluidInterior) // Loop through interior fluid nodes
		{
			Vector3_i indices3D = getIndices3D(index1D, subMesh.arrayLimits);
			int i = indices3D.i, j = indices3D.j, k = indices3D.k;
			// To make the colossal expression below a bit more readable I define some subindices for neighbor nodes:
			// P(actual point), W(West, i-1), E(East, i+1), S(South, j-1), N(North, j+1), D(Down, k-1), U(Up, k+1)
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

// Compute intermediate solutions for all conserved variables, at given time increment after current time,
// using the given set of slopes.
void Solver::computeAllIntermediateSolutions(ChooseSlopeStage slopeStage, // ← Slopes to use
											 double timeIncrement) // ← How long after current time to evaluate
{
	for(SubMesh& subMesh : mesh.subMeshes)
	{
		const ConservedVariablesArrayGroup* RK4Slopes = nullptr;
		switch(slopeStage)
		{
		case ChooseSlopeStage::k1:
			RK4Slopes = &subMesh.RK4slopes.k1;	break;
		case ChooseSlopeStage::k2:
			RK4Slopes = &subMesh.RK4slopes.k2;	break;
		case ChooseSlopeStage::k3:
			RK4Slopes = &subMesh.RK4slopes.k3;	break;
		case ChooseSlopeStage::k4:
			RK4Slopes = &subMesh.RK4slopes.k4;	break;
		default:
			throw std::logic_error("Unexpected enum value for 'ChooseSlopeStage' when computing intermediate solutions.");
		}
		const vector<int>& interiorNodes = subMesh.indexByType.fluidInterior;
		const ConservedVariablesArrayGroup& previousLevel = subMesh.conservedVariablesOld;
		ConservedVariablesArrayGroup& intermediateSolution = subMesh.conservedVariables;
		computeIntermediateSolution(previousLevel.rho,   RK4Slopes->rho,   timeIncrement, interiorNodes, intermediateSolution.rho);
		computeIntermediateSolution(previousLevel.rho_u, RK4Slopes->rho_u, timeIncrement, interiorNodes, intermediateSolution.rho_u);
		computeIntermediateSolution(previousLevel.rho_v, RK4Slopes->rho_v, timeIncrement, interiorNodes, intermediateSolution.rho_v);
		computeIntermediateSolution(previousLevel.rho_w, RK4Slopes->rho_w, timeIncrement, interiorNodes, intermediateSolution.rho_w);
		computeIntermediateSolution(previousLevel.rho_E, RK4Slopes->rho_E, timeIncrement, interiorNodes, intermediateSolution.rho_E);
	}
}

// Evaluates an intermediate solution of the given conserved variable, using the given RK4 step (slope)
// and the appropriate time increment (dt/2 if k1 or k2 is given, dt if k3 is given).
// Result is stored in 'intermSolution'. Only writes interior fluid nodes.
void Solver::computeIntermediateSolution(const Array3D_d& conservedVar, // ← Old values
										 const Array3D_d& RK4_slope,
										 double timeIncrement,
										 const vector<int>& interiorNodes,
										 Array3D_d& intermSolution) // ← OUTPUT, new values
{
	for(int index1D : interiorNodes)
		intermSolution(index1D) = conservedVar(index1D) + timeIncrement * RK4_slope(index1D);
}

// Updates the primitive variables and transport properties, from the conserved variables.
// Only writes interior fluid nodes!
void Solver::updatePrimitiveVariables()
{
	double gammaMinusOne = params.Gamma - 1;
	double sutherlands_Sc = params.sutherlands_C2 / params.T_0;
	double ScPlusOne = 1 + sutherlands_Sc;
	double prandtlFactor = 1 / ( gammaMinusOne * params.Pr );
	for(SubMesh& subMesh : mesh.subMeshes)
	{
		const Array3D_d& rho_u{subMesh.conservedVariables.rho_u};
		const Array3D_d& rho_v{subMesh.conservedVariables.rho_v};
		const Array3D_d& rho_w{subMesh.conservedVariables.rho_w};
		const Array3D_d& rho  {subMesh.conservedVariables.rho  };
		const Array3D_d& rho_E{subMesh.conservedVariables.rho_E};
		Array3D_d& u{subMesh.primitiveVariables.u};
		Array3D_d& v{subMesh.primitiveVariables.v};
		Array3D_d& w{subMesh.primitiveVariables.w};
		Array3D_d& p{subMesh.primitiveVariables.p};
		Array3D_d& T{subMesh.primitiveVariables.T};
		Array3D_d& mu   {subMesh.transportProperties.mu};
		Array3D_d& kappa{subMesh.transportProperties.kappa};

		for(int i : subMesh.indexByType.fluidInterior)
			u(i) = rho_u(i) / ( rho(i) + 1 );
		for(int i : subMesh.indexByType.fluidInterior)
			v(i) = rho_v(i) / ( rho(i) + 1 );
		for(int i : subMesh.indexByType.fluidInterior)
			w(i) = rho_w(i) / ( rho(i) + 1 );
		for(int i : subMesh.indexByType.fluidInterior)
			p(i) = gammaMinusOne * ( rho_E(i) - (rho(i)+1)/2 * (u(i)*u(i) + v(i)*v(i) + w(i)*w(i)) );
		for(int i : subMesh.indexByType.fluidInterior)
			T(i) = ( params.Gamma * p(i) - rho(i) ) / ( 1+rho(i) );
		for(int i : subMesh.indexByType.fluidInterior)
			mu(i) = pow( 1+T(i), 1.5 ) * ScPlusOne / ( params.Re_0*( T(i) + ScPlusOne ) );
		for(int i : subMesh.indexByType.fluidInterior)
			kappa(i) = mu(i) * prandtlFactor;
	}
}

// Do the final step of the RK4 method for all conserved variables.
// Use all the 4 slopes to bring the solution to the next time level
void Solver::computeRK4finalStepAllVariables()
{
	for(SubMesh& subMesh : mesh.subMeshes)
	{
		const ConservedVariablesArrayGroup& k1{subMesh.RK4slopes.k1};
		const ConservedVariablesArrayGroup& k2{subMesh.RK4slopes.k2};
		const ConservedVariablesArrayGroup& k3{subMesh.RK4slopes.k3};
		const ConservedVariablesArrayGroup& k4{subMesh.RK4slopes.k4};
		const ConservedVariablesArrayGroup& conservedVarsOld{subMesh.conservedVariablesOld};
		ConservedVariablesArrayGroup& conservedVarsNew{subMesh.conservedVariables};
		const vector<int>& interiorNodes{subMesh.indexByType.fluidInterior};
		compute_RK4_final_step(k1.rho,   k2.rho,   k3.rho,   k4.rho,   interiorNodes, conservedVarsOld.rho  , conservedVarsNew.rho  );  // Use a weighted average of the four slopes
		compute_RK4_final_step(k1.rho_u, k2.rho_u, k3.rho_u, k4.rho_u, interiorNodes, conservedVarsOld.rho_u, conservedVarsNew.rho_u);  // to advance the solution from time t, to
		compute_RK4_final_step(k1.rho_v, k2.rho_v, k3.rho_v, k4.rho_v, interiorNodes, conservedVarsOld.rho_v, conservedVarsNew.rho_v);  // t + dt. The solutions at the new time level
		compute_RK4_final_step(k1.rho_w, k2.rho_w, k3.rho_w, k4.rho_w, interiorNodes, conservedVarsOld.rho_w, conservedVarsNew.rho_w);  // are stored in the intermediate arrays, so the
		compute_RK4_final_step(k1.rho_E, k2.rho_E, k3.rho_E, k4.rho_E, interiorNodes, conservedVarsOld.rho_E, conservedVarsNew.rho_E);  // norm of the change can be computed.
	}
}

// Uses all the intermediate RK4 slopes k1,2,3,4, for a variable, to advance the variable from t to t+dt
void Solver::compute_RK4_final_step(const Array3D_d& k1, // ← Slopes
									const Array3D_d& k2,
									const Array3D_d& k3,
									const Array3D_d& k4,
									const vector<int>& interiorNodes,
									const Array3D_d& conservedVar_old, // ← Values at current time
									Array3D_d& conservedVar_new) // ← OUTPUT, values at next time level
{
	for(int i : interiorNodes)
	{
		double residualValue = (dt/6) * ( k1(i) + 2*k2(i) + 2*k3(i) + k4(i) );
		conservedVar_new(i) = conservedVar_old(i) + residualValue;
	}
}

// Compute and store norms of change for conserved variables, if specified in config file
void Solver::storeNormsOfChange()
{

	if(params.saveConvergenceHistory != ConvHistoryEnum::none)
	{
		double normRho = 0;
		for(const SubMesh& subMesh : mesh.subMeshes)
			normRho += subMesh.getNormOfChange(subMesh.conservedVariablesOld.rho, subMesh.conservedVariables.rho);
		normHistory.rho.push_back(normRho);
	}

	if(params.saveConvergenceHistory == ConvHistoryEnum::all)
	{
		double normRhoU = 0;
		double normRhoV = 0;
		double normRhoW = 0;
		double normRhoE = 0;
		for(const SubMesh& subMesh : mesh.subMeshes)
		{
			const ConservedVariablesArrayGroup& conservedVarsOld = subMesh.conservedVariablesOld;
			const ConservedVariablesArrayGroup& conservedVarsNew = subMesh.conservedVariables;
			normRhoU += subMesh.getNormOfChange(conservedVarsOld.rho_u, conservedVarsNew.rho_u);
			normRhoV += subMesh.getNormOfChange(conservedVarsOld.rho_v, conservedVarsNew.rho_v);
			normRhoW += subMesh.getNormOfChange(conservedVarsOld.rho_w, conservedVarsNew.rho_w);
			normRhoE += subMesh.getNormOfChange(conservedVarsOld.rho_E, conservedVarsNew.rho_E);
		}
		normHistory.rho_u.push_back(normRhoU);
		normHistory.rho_v.push_back(normRhoV);
		normHistory.rho_w.push_back(normRhoW);
		normHistory.rho_E.push_back(normRhoE);
	}
}








