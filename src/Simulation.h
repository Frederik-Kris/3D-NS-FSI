/*
 * Simulation.h
 *
 *  Created on: Oct 6, 2022
 *      Author: frederk
 */

#ifndef SRC_SIMULATION_H_
#define SRC_SIMULATION_H_

#include "Array3D.h"
#include "includes_and_names.h"
#include "Solver.h"
#include "OutputManager.h"
#include "ConfigSettings.h"

// Top level class, containing all variables and methods needed to run a simulation.
class Simulation
{
public:

	Simulation(const ConfigSettings& params);

	void run();

private:

	bool checkStoppingCriterion();

	const ConfigSettings params;// Parameters and settings, imported from config file
	Solver solver;				// Solver class, for handling numerical method.
	OutputManager output;		// Class for managing output to screen and files.
	double t;					// Simulated time
	ulong timeLevel;			// No. of timesteps computed. Zero is IC.
	ulong timeLevelStart;		// Time level to start from, if simulation is continued from previous result.
};



#endif /* SRC_SIMULATION_H_ */
