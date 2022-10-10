/*
 * Simulation.cpp
 *
 *  Created on: Oct 6, 2022
 *      Author: frederk
 */

#include "Simulation.h"

Simulation::Simulation() :
params(),
solver(params),
output(),
t{0}, timeLevel{0}
{
	output.processInitialOutput(params);
	Clock statusReportTimer;
	while ( !checkStoppingCriterion() )
	{
		solver.marchTimeStep();
		output.processIntermediateOutput(params, statusReportTimer, t, solver.dt);
	}
	const vector<const vector<double>*> convergenceHistoryPointers = {&normHistory_rho, &normHistory_rho_u, &normHistory_rho_v, &normHistory_rho_w, &normHistory_E};
	output.processFinalOutput(params);
}

// Check if computing another timestep will violate the stopping criterion. Returns true if simulation should stop.
// If the the criterion is to stop after a given number of timesteps, no more timesteps are required after this.
// To stop at an exact time, one last timestep with an adapted dt is required to hit the end-time.
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


