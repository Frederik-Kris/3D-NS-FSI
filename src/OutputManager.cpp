/*
 * OutputManager.cpp
 *
 *  Created on: Oct 7, 2022
 *      Author: frederk
 */

#include "OutputManager.h"

OutputManager::OutputManager()
{
}

// Checks whether it is time to store the solution to disk, or write out a status report to screen and does it if appropriate.
// It will save solution if the next timestep would take the solution past the save-time.
// Thus, in general it saves too early, but the deviation is dt at most.
void Solver::processOutput(Clock& statusReportTimer)
{
	double timeSinceSave = fmod(t, params.save_period);
	if ( timeSinceSave + dt >= params.save_period && params.save_intervals )
		storeCurrentSolution_csv();
	double timeSinceStatusReport = statusReportTimer.getElapsedTime().asSeconds();
	if ( timeSinceStatusReport >= params.statusReportInterval )
	{
		writeStatusReport_toScreen();
		statusReportTimer.restart();
	}
}

// Check whether the influx and outflux of mass are equal, by integrating momentum over inlet and outlet boundaries.
// Integration is done by double trapezoidal rule, i.e., a double loop summing the nodes in the inlet and outlet planes, wighting the outer nodes with 1/2.
void Solver::checkMassConservation(double& inFluxSum, double& outFluxSum)
{
	inFluxSum = outFluxSum = 0;
	for (uint j=0; j<params.NJ; ++j)
	{
		double inletColumnSum = 0, outletColumnSum = 0;
		inletColumnSum  += rho_u(0          ,j,0)/2 + rho_u(0          ,j,params.NK-1)/2;	// Add the boundary nodes with weigth 1/2
		outletColumnSum += rho_u(params.NI-1,j,0)/2 + rho_u(params.NI-1,j,params.NK-1)/2;
		for(uint k=1; k<params.NK-1; ++k)
		{
			inletColumnSum  += rho_u(0,j,k);
			outletColumnSum += rho_u(params.NI-1,j,k);
		}
		inletColumnSum  *= dz;
		outletColumnSum *= dz;
		double weightMultiplier = (j==0 || j==params.NJ-1)? 0.5 : 1.;
		inFluxSum  += weightMultiplier * inletColumnSum;
		outFluxSum += weightMultiplier * outletColumnSum;
	}
}

// Store selected variables from the solution at current time level, using the format(s) specified
void Solver::storeCurrentSolution_csv()
{
	if(params.saveForParaview)
		storeCurrentSolution_csv_paraview();
	if (params.saveForMatlab)
		storeCurrentSolution_csv_matlab();
	++savedSolutions;
	outputTimes.push_back(t);
}

// Writes a .csv file with the specified flow variables at the current time level.
// The file is formatted to fit how ParaView wants to import it, with coordinates in the leftmost column.
// Only the flow variables selected in the ConfigFile are saved.
void Solver::storeCurrentSolution_csv_paraview()
{
	vector<Array3D_d*> flowVariables = getPlotVariables();
	if ( flowVariables.empty() )
		return;

	ofstream outputFile;
	string filename = "output/out.csv." + to_string(savedSolutions);
	outputFile.open( filename );
	if ( !outputFile )
	{
		cout << "Could not open file:  " + filename << endl
		     << "Solution was not saved. You may have to move the 'output' folder to the location of the source files or to the executable itself, depending on whether you run the executable through an IDE or by itself." << endl;
		return;
	}
	outputFile << get_csvHeaderString();

	for (uint i{0}; i<params.NI; ++i)
		for (uint j{0}; j<params.NJ; ++j)
			for (uint k{0}; k<params.NK; ++k)
			{
				outputFile << endl;
				double x{ i*dx }, y{ j*dy }, z{ k*dz};
				outputFile << x << ", " << y << ", " << z;
				for (Array3D_d* flowVar : flowVariables)
					outputFile << ", " << (*flowVar)(i,j,k);
			}
	outputFile.close();
}

// Writes multiple .csv files with the specified flow variables at the current time level.
// The file is formatted as a table or matrix, with the values from a specified 2D plane. Matlab
// can read this file to a matrix. Only the flow variables selected in the ConfigFile are saved.
void Solver::storeCurrentSolution_csv_matlab()
{
	vector<Array3D_d*> flowVariables = getPlotVariables();		// Get pointers to the arrays with data to save
	vector<string> variableFileNames = getVariableFileNames();	// Get a vector with the names of the variables

	for(uint i=0; i<flowVariables.size(); ++i)	// Let i go from zero to no. of variables to save.
	{
		ofstream outputFile;
		string filename = "output/" + variableFileNames.at(i) + "_" + to_string(savedSolutions) + ".csv";
		outputFile.open( filename );
		if ( !outputFile )
		{
			cout << "Could not open file:  " + filename << endl
			     << "Solution was not saved. You may have to move the 'output' folder to the location of the source files or to the executable itself, depending on whether you run the executable through an IDE or by itself." << endl;
			return;
		}
		writePlaneTo_csv(outputFile, flowVariables.at(i));
	}
}

// Get a vector with pointers to the flow variables that should be saved.
vector<Array3D_d*> Solver::getPlotVariables()
{
	vector<Array3D_d*> flowVariables;
	if ( params.save_rho )
		flowVariables.push_back(&rho);
	if ( params.save_rho_u )
		flowVariables.push_back(&rho_u);
	if ( params.save_rho_v )
		flowVariables.push_back(&rho_v);
	if ( params.save_rho_w )
		flowVariables.push_back(&rho_w);
	if ( params.save_E )
		flowVariables.push_back(&E);
	if ( params.save_u )
		flowVariables.push_back(&u);
	if ( params.save_v )
		flowVariables.push_back(&v);
	if ( params.save_w )
		flowVariables.push_back(&w);
	if ( params.save_p )
		flowVariables.push_back(&p);
	if ( params.save_T )
		flowVariables.push_back(&T);
	if ( params.save_mu )
		flowVariables.push_back(&mu);
	if ( params.save_kappa )
		flowVariables.push_back(&kappa);
	return flowVariables;
}

// Returns a camma-separated string with the names of the flow variables to save.
// For ParaView to read a .csv file it needs the headers on the first line.
string Solver::get_csvHeaderString()
{
	string headers = "x, y, z";
	if ( params.save_rho )
		headers += ", density";
	if ( params.save_rho_u )
		headers += ", x-momentum";
	if ( params.save_rho_v )
		headers += ", y-momentum";
	if ( params.save_rho_w )
		headers += ", z-momentum";
	if ( params.save_E )
		headers += ", total energy";
	if ( params.save_u )
		headers += ", velocity comp u";
	if ( params.save_v )
		headers += ", velocity comp v";
	if ( params.save_w )
		headers += ", velocity comp w";
	if ( params.save_p )
		headers += ", pressure";
	if ( params.save_T )
		headers += ", temperature";
	if ( params.save_mu )
		headers += ", dynamic viscosity";
	if ( params.save_kappa )
		headers += ", thermal conductivity";
	return headers;
}

// Returns a vector of strings, which are the names of the variables to save.
vector<string> Solver::getVariableFileNames()
{
	vector<string> variableNames;
	if ( params.save_rho )
		variableNames.push_back("rho");
	if ( params.save_rho_u )
		variableNames.push_back("rho_u");
	if ( params.save_rho_v )
		variableNames.push_back("rho_v");
	if ( params.save_rho_w )
		variableNames.push_back("rho_w");
	if ( params.save_E )
		variableNames.push_back("E");
	if ( params.save_u )
		variableNames.push_back("u");
	if ( params.save_v )
		variableNames.push_back("v");
	if ( params.save_w )
		variableNames.push_back("w");
	if ( params.save_p )
		variableNames.push_back("p");
	if ( params.save_T )
		variableNames.push_back("T");
	if ( params.save_mu )
		variableNames.push_back("mu");
	if ( params.save_kappa )
		variableNames.push_back("kappa");
	return variableNames;
}

// Write the values from a plane, defined in ConfigFile, of one flow variable. It's written like a comma separated table into 'outputFile'.
void Solver::writePlaneTo_csv(ofstream& outputFile, Array3D_d* flowVariable)
{
	switch(params.saveNormalAxis)
	{
	case saveNormalAxisEnum::x:
		for(uint j=0; j<params.NJ; ++j)
		{
			outputFile << (*flowVariable)(params.saveConstantIndex, j, 0);
			for(uint k=1; k<params.NK; ++k)
				outputFile << ", " << (*flowVariable)(params.saveConstantIndex, j, k);
			outputFile << endl;
		}
		break;
	case saveNormalAxisEnum::y:
		for(uint i=0; i<params.NI; ++i)
		{
			outputFile << (*flowVariable)(i, params.saveConstantIndex, 0);
			for(uint k=1; k<params.NK; ++k)
				outputFile << ", " << (*flowVariable)(i, params.saveConstantIndex, k);
			outputFile << endl;
		}
		break;
	case saveNormalAxisEnum::z:
		for(uint i=0; i<params.NI; ++i)
		{
			outputFile << (*flowVariable)(i, 0, params.saveConstantIndex);
			for(uint j=1; j<params.NJ; ++j)
				outputFile << ", " << (*flowVariable)(i, j, params.saveConstantIndex);
			outputFile << endl;
		}
		break;
	}
}

// Writes a brief report with progression, timestep size and wall clock time elapsed, etc.
// 'setprecision' is used to control no. of significant digits, default is 6.
void Solver::writeStatusReport_toScreen()
{
	cout << "Simulated time: t = " << t;
	if (params.stopCriterion == StopCriterionEnum::end_time)
	{
		double progressPercentage = t / params.t_end * 100.;
		cout << " , t_end = " << params.t_end << " ( " << setprecision(3) << progressPercentage << setprecision(6) << " % )" << endl;
	}
	else
		cout << endl;

	cout << "Time level: n = " << timeLevel;
	if (params.stopCriterion == StopCriterionEnum::timesteps)
	{
		double progressPercentage = (double)timeLevel / (double)params.stopTimeLevel * 100.;  // Cast to double to avoid integer division.
		cout << " , n_max = " << params.stopTimeLevel << " ( " << setprecision(3) << progressPercentage << setprecision(6) << " % )" << endl;
	}
	else
		cout << endl;
	cout << "Timestep size: dt = " << dt << endl;
	double inletMassFlux=0, outletMassFlux=0;
	checkMassConservation(inletMassFlux, outletMassFlux);
	cout << "Inlet mass flux: " << inletMassFlux << " , Outlet mass flux: " << outletMassFlux << " , Difference: " << inletMassFlux-outletMassFlux << endl;
	cout << "Wall clock time: " << setprecision(3) << wallClockTimer.getElapsedTime().asSeconds() << setprecision(6) << " sec" << endl << endl;
}

// Write a file with a list of the actual times when the solution was saved. For debugging.
void Solver::writeOutputTimes()
{
	ofstream timeFile;
	string filename = "output/times.dat";
	timeFile.open( filename );
	if ( !timeFile )
	{
		cout << "Could not open file:  " + filename << endl
		     << "Output times were not written to .dat file. You may have to move the 'output' folder to the location of the source files or to the executable itself, depending on whether you run the executable through an IDE or by itself." << endl;
		return;
	}
	for (double time : outputTimes)
		timeFile << time << endl;
	timeFile.close();
}

// Write files with lists of the norm of change for the conserved variables.
// TODO: This function should be customizable to choose in ConfigFile whether to log norms of particular variables. OR, find out which variable is the best indicator for convergence, and only log that one.
void Solver::writeNormHistoryFiles()
{
	vector<string> filenames = {"output/norm_rho.dat", "output/norm_rho_u.dat", "output/norm_rho_v.dat", "output/norm_rho_w.dat", "output/norm_E.dat"};
	vector<vector<double>*> normVectorPointers = {&normHistory_rho, &normHistory_rho_u, &normHistory_rho_v, &normHistory_rho_w, &normHistory_E};
	for(uint i=0; i<filenames.size(); ++i)
	{
		ofstream normFile;
		normFile.open( filenames.at(i) );
		if ( !normFile )
		{
			cout << "Could not open file:  " + filenames.at(i) << endl
				 << "Norm history was not written to .dat file. You may have to move the 'output' folder to the location of the source files or to the executable itself, depending on whether you run the executable through an IDE or by itself." << endl;
		}
		else
			for (double norm : *normVectorPointers.at(i))
				normFile << norm << endl;
		normFile.close();
	}
}






