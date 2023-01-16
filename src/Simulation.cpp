/*
 * Simulation.cpp
 *
 *  Created on: Oct 6, 2022
 *      Author: frederk
 */

#include "Simulation.h"

// Constructor taking parameters from config file
Simulation::Simulation(const ConfigSettings& params) :
params(params),
solver(params),
output(params),
t{0}, timeLevel{0}
{}

// Initialize and run simulation, and handle output as specified in config file
void Simulation::run()
{
	solver.initialize();
	output.initialize();
	output.processInitialOutput(solver.mesh, t);
	while ( !checkStoppingCriterion() )
	{
		solver.marchTimeStep(t, timeLevel);
		output.processIntermediateOutput(solver.mesh, t, timeLevel, solver.dt);
	}
	output.processFinalOutput(solver.mesh, t, timeLevel, solver.dt, solver.getConvergenceHistory());
}

// Check if computing another timestep will violate the stopping criterion. Returns true if simulation should stop.
bool Simulation::checkStoppingCriterion()
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


