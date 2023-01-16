/*
 * solver.h
 *
 *  Created on: Apr 28, 2021
 *      Author: frederik
 */

#ifndef SRC_SOLVER_H_
#define SRC_SOLVER_H_

#include "Array3D.h"
#include "includes_and_names.h"
#include "Mesh.h"
#include "ConfigSettings.h"

// Numerical solver class, which contains mesh with flow variables and methods needed to march in time.
class Solver
{
public:

	Solver(const ConfigSettings& params);

	void initialize();

	void marchTimeStep(double& t, ulong& timeLevel);

	const ConservedVariablesVectorGroup& getConvergenceHistory() const {return normHistory;}

	const ConfigSettings params;	// Parameters and settings, imported from ConfigFile
	Mesh mesh;						// Computational mesh, containing flow variables
	double dt;						// Time-step size

private:
	void applyStagnation_IC();

	void applyUniformFlow_IC();

	void updateTimeStepSize(double t);

	double getInviscidTimeStepLimit();

	double getViscousTimeStepLimit();

	void computeRK4slopes(const ConservedVariablesArrayGroup& conservedVariables,
						  ConservedVariablesArrayGroup& RK4slopes);

	void compute_RK4_step_continuity(const Array3D_d& rho_u,
									 const Array3D_d& rho_v,
									 const Array3D_d& rho_w,
									 Array3D_d& RK4_slope);

	void compute_RK4_step_xMomentum(const Array3D_d& rho_u, Array3D_d& RK4_slope);

	void compute_RK4_step_yMomentum(const Array3D_d& rho_v, Array3D_d& RK4_slope);

	void compute_RK4_step_zMomentum(const Array3D_d& rho_w, Array3D_d& RK4_slope);

	void compute_RK4_step_energy(const Array3D_d& E, Array3D_d& RK4_slope);

	void computeAllIntermediateSolutions(const ConservedVariablesArrayGroup& RK4slopes, double timeIncrement);

	void computeIntermediateSolution(const Array3D_d& conservedVar,
									 const Array3D_d& RK4_slope,
									 double timeIncrement,
									 Array3D_d& intermSolution);

	void updatePrimitiveVariables();

	void computeRK4finalStepAllVariables();

	void compute_RK4_final_step(const Array3D_d& k1,
								const Array3D_d& k2,
								const Array3D_d& k3,
								const Array3D_d& k4,
								const Array3D_d& conservedVar_old,
								Array3D_d& conservedVar_new);

	void storeNormsOfChange();

	ConservedVariablesVectorGroup normHistory; // Vectors of the developments of the change-norms.
};


#endif /* SRC_SOLVER_H_ */
