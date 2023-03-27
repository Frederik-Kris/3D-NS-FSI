/*
 * main.cpp
 *
 *  Created on: Apr 8, 2021
 *      Author: frederik
 */

#include "includes_and_names.h"
#include "Simulation.h"
#include "tests.h"

int main()
{
	const ConfigSettings params("ConfigFile");	// Read parameters from config file.
	if(params.errorOccurred)
		return -1;
	Simulation sim(params);						// Initialize simulation with given parameters.
	sim.run();									// Run simulation. Produce output as specified.
	return 0;
}
