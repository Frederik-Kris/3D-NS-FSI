/*
 * solver.h
 *
 *  Created on: Apr 28, 2021
 *      Author: frederik
 */

#ifndef SRC_SOLVER_H_
#define SRC_SOLVER_H_

#include "includes_and_names.h"
#include "Mesh.h"
#include "Array3D_d.h"
#include "ConfigSettings.h"

// Numerical solver class, which contains mesh with flow variables and methods needed to march in time.
class Solver
{
public:
	Solver(const ConfigSettings& params);
private:
	void applyStagnation_IC();
	void applyUniformFlow_IC();
	void updateTimeStepSize(const Array3D_d& p, const Array3D_d& rho, const Array3D_d& u, const Array3D_d& v, const Array3D_d& w);
	void getDerivedVariables_atPoint(double u_0,double v_0, double w_0, double p_0, double T_0,
			double& rho_0, double& rho_u_0, double& rho_v_0, double& rho_w_0, double& E_0, double& mu_0, double& kappa_0);
	void marchTimeStep();
	void compute_RK4_step_continuity(const Array3D_d& rho_u, const Array3D_d& rho_v, const Array3D_d& rho_w,
                                           Array3D_d& RK4_slope);
	void compute_RK4_step_xMomentum(const Array3D_d& rho_u, Array3D_d& RK4_slope);
	void compute_RK4_step_yMomentum(const Array3D_d& rho_v, Array3D_d& RK4_slope);
	void compute_RK4_step_zMomentum(const Array3D_d& rho_w, Array3D_d& RK4_slope);
	void compute_RK4_step_energy(const Array3D_d& E, Array3D_d& RK4_slope);
	void computeIntermediateSolution(const Array3D_d& conservedVar, const Array3D_d& RK4_slope,
			                               Array3D_d& intermSolution, double timeIncrement);
	void updatePrimitiveVariables(const Array3D_d& rho, const Array3D_d& rho_u, const Array3D_d& rho_v,
                                  const Array3D_d& rho_w, const Array3D_d& E);
	void applyInjectionBC_risingInletVelocityChannelFlow(Array3D_d& rho, Array3D_d& rho_u, Array3D_d& rho_v, Array3D_d& rho_w, Array3D_d& E, double time);
	void compute_RK4_final_step(const Array3D_d& k1, const Array3D_d& k2, const Array3D_d& k3, const Array3D_d& k4,
								const Array3D_d& conservedVar_old, Array3D_d& conservedVar_new);

	ConfigSettings params;								// Parameters and settings, imported from ConfigFile
	Mesh mesh;											// Computational mesh, containing flow variables
	double dt;											// Time-step size
};


#endif /* SRC_SOLVER_H_ */
