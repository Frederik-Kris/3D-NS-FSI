/*
 * Simulation.h
 *
 *  Created on: Oct 6, 2022
 *      Author: frederk
 */

#ifndef SRC_SIMULATION_H_
#define SRC_SIMULATION_H_

#include "includes_and_names.h"
#include "Solver.h"
#include "OutputManager.h"
#include "Array3D_d.h"
#include "ConfigSettings.h"

// Top level class, containing all variables and methods needed to run a simulation.
class Simulation
{
public:
	Simulation();
private:

	bool checkStoppingCriterion();

	ConfigSettings params;		// Parameters and settings, imported from ConfigFile
	Solver solver;				// Solver class, for handling numerical method.
	OutputManager output;		// Class for managing output to screen and files.
	vector<double> normHistory_rho,  normHistory_rho_u, // Vectors of the development of the change-norms.
			normHistory_rho_v, normHistory_rho_w, normHistory_E;
	double t;		// Simulated time
	uint timeLevel;                                     // No. of timesteps computed. Zero is IC.
	Clock wallClockTimer;								// Timer that starts when simulation starts
};



#endif /* SRC_SIMULATION_H_ */
