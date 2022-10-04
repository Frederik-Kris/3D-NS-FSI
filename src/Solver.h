/*
 * solver.h
 *
 *  Created on: Apr 28, 2021
 *      Author: frederik
 */

#ifndef SRC_SOLVER_H_
#define SRC_SOLVER_H_



#include "includes_and_names.h"
#include "Array3D_d.h"
#include "ConfigSettings.h"

// Top level class, which contains all flow variables and methods needed to run a simulation
class Solver
{
public:
	Solver();
private:
	void setGridSpacings();
	void applyStagnation_IC();
	void applyUniformFlow_IC();
	void getDerivedVariables_atPoint(double u_0,double v_0, double w_0, double p_0, double T_0,
			double& rho_0, double& rho_u_0, double& rho_v_0, double& rho_w_0, double& E_0, double& mu_0, double& kappa_0);
	void updateTimeStepSize();
	bool checkStoppingCriterion();
	void marchTimeStep();
	void applyFilter_ifAppropriate(Array3D_d& variable_old, Array3D_d& variable_new);
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
	void computeNorms_conservedVariables();
	double getNormOfChange(const Array3D_d& oldValue, const Array3D_d& newValue);
	void swapConservedVariables();
	void processOutput(Clock& statusReportTimer);
	void storeCurrentSolution_csv();
	void storeCurrentSolution_csv_paraview();
	void storeCurrentSolution_csv_matlab();
	vector<Array3D_d*> getPlotVariables();
	string get_csvHeaderString();
	vector<string> getVariableFileNames();
	void writePlaneTo_csv(ofstream& outputFile, Array3D_d* flowVariable);
	void writeStatusReport_toScreen();
	void writeOutputTimes();
	void writeNormHistoryFiles();
	void writeChannelFlowResults();
	void checkMassConservation(double& inFluxSum, double& outFluxSum);

	ConfigSettings params;								// Parameters and settings, imported from ConfigFile
	Array3D_d rho, rho_u, rho_v, rho_w, E;				// Conserved variables
	Array3D_d u, v, w, p, T;							// Primitive variables
	Array3D_d mu, kappa;								// Transport properties
	Array3D_d k1_rho  , k2_rho  , k3_rho  , k4_rho  ;	// RK4 stages, continuity
	Array3D_d k1_rho_u, k2_rho_u, k3_rho_u, k4_rho_u;	// RK4 stages, x-momentum
	Array3D_d k1_rho_v, k2_rho_v, k3_rho_v, k4_rho_v;	// RK4 stages, y-momentum
	Array3D_d k1_rho_w, k2_rho_w, k3_rho_w, k4_rho_w;	// RK4 stages, z-momentum
	Array3D_d k1_E    , k2_E    , k3_E    , k4_E    ;	// RK4 stages, specific total energy
	Array3D_d interm_rho, interm_rho_u, interm_rho_v,	// Intermediate solutions.
			  interm_rho_w, interm_E;
	vector<double> normHistory_rho,  normHistory_rho_u, // Vectors of the development of the change-norms.
		normHistory_rho_v, normHistory_rho_w, normHistory_E;
	uint savedSolutions;                                // No. of times saved to disk
	vector<double> outputTimes;                         // The exact times when solution was saved
	uint timeLevel;                                     // No. of timesteps computed. Zero is IC.
	double t, dt;										// Current simulated time and time step size
	double dx, dy, dz;									// Grid spacing in x-, y- and z-direction
	Clock wallClockTimer;								// Timer that starts when simulation starts
};




#endif /* SRC_SOLVER_H_ */
