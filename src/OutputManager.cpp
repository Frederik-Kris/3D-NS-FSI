/*
 * OutputManager.cpp
 *
 *  Created on: Oct 7, 2022
 *      Author: frederk
 */

#include "OutputManager.h"

OutputManager::OutputManager(const ConfigSettings& params) :
params(params),
savedSolutions{0}
{
}

void OutputManager::processInitialOutput(const Mesh& mesh, double t)
{
	if ( params.saveIC )
		storeCurrentSolution_csv(mesh, t);
}

// Checks whether it is time to store the solution to disk, or write out a status report to screen and does it if appropriate.
// It will save solution if the next timestep would take the solution past the save-time.
// Thus, in general it saves too early, but the deviation is dt at most.
void OutputManager::processIntermediateOutput(const Mesh& mesh, Clock& statusReportTimer, double t, ulong timeLevel, double dt)
{
	// Save to disk:
	if(params.saveIntervals)
	{
		bool withinRelevantTimeRegime = t >= params.saveIntervalsStartTime
								   && ( t <= params.saveIntervalsEndTime || params.saveIntervalsEndTime <= 0);
		bool enoughTimeSinceSave = true;
		if(params.savePeriod > 0)
		{
			double timeSinceSave = fmod(t, params.savePeriod);
			enoughTimeSinceSave = timeSinceSave + dt >= params.savePeriod;
		}
		if ( withinRelevantTimeRegime && enoughTimeSinceSave )
			storeCurrentSolution_csv(mesh, t);
	}
	// Write to screen:
	double timeSinceStatusReport = statusReportTimer.getElapsedTime().asSeconds();
	if ( timeSinceStatusReport >= params.statusReportInterval )
	{
		writeStatusReport_toScreen(t, timeLevel, dt);
		statusReportTimer.restart();
	}
}

void OutputManager::processFinalOutput(const Mesh& mesh, double t, ulong timeLevel, double dt,
								const vector<ConservedVariablesScalars> convergenceHistory)
{
	if ( params.saveFinal )
		storeCurrentSolution_csv(mesh, t);
	writeStatusReport_toScreen(t, timeLevel, dt);
	writeOutputTimes();
	writeConvergenceHistoryFiles(convergenceHistory);
}

// Store selected variables from the solution at current time level, using the format(s) specified
void OutputManager::storeCurrentSolution_csv(const Mesh& mesh, double t)
{
	if(params.saveForParaview)
		storeCurrentSolution_csv_paraview(mesh);
	if (params.saveForMatlab)
		storeCurrentSolution_csv_matlab(mesh);
	++savedSolutions;
	outputTimes.push_back(t);
}

// Writes a .csv file with the specified flow variables at the current time level.
// The file is formatted to fit how ParaView wants to import it, with coordinates in the leftmost column.
// Only the flow variables selected in the ConfigFile are saved.
void OutputManager::storeCurrentSolution_csv_paraview(const Mesh& mesh)
{
	vector<const Array3D_d*> flowVariables = getPlotVariables(mesh);
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

	for (size_t i{0}; i<params.NI; ++i)
		for (size_t j{0}; j<params.NJ; ++j)
			for (size_t k{0}; k<params.NK; ++k)
			{
				outputFile << endl;
				double x{ i*mesh.dx }, y{ j*mesh.dy }, z{ k*mesh.dz};
				outputFile << x << ", " << y << ", " << z;
				for (const Array3D_d* flowVar : flowVariables)
					outputFile << ", " << (*flowVar)(i,j,k);
			}
	outputFile.close();
}

// Writes multiple .csv files with the specified flow variables at the current time level.
// The file is formatted as a table or matrix, with the values from a specified 2D plane. Matlab
// can read this file to a matrix. Only the flow variables selected in the ConfigFile are saved.
void OutputManager::storeCurrentSolution_csv_matlab(const Mesh& mesh)
{
	vector<const Array3D_d*> flowVariables = getPlotVariables(mesh);		// Get pointers to the arrays with data to save
	vector<string> variableFileNames = getVariableFileNames();	// Get a vector with the names of the variables

	for(size_t i=0; i<flowVariables.size(); ++i)	// Let i go from zero to no. of variables to save.
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
vector<const Array3D_d*> OutputManager::getPlotVariables(const Mesh& mesh)
{
	vector<const Array3D_d*> flowVariables;
	if ( params.save_rho )
		flowVariables.push_back(&mesh.conservedVariables.rho);
	if ( params.save_rho_u )
		flowVariables.push_back(&mesh.conservedVariables.rho_u);
	if ( params.save_rho_v )
		flowVariables.push_back(&mesh.conservedVariables.rho_v);
	if ( params.save_rho_w )
		flowVariables.push_back(&mesh.conservedVariables.rho_w);
	if ( params.save_E )
		flowVariables.push_back(&mesh.conservedVariables.rho_E);
	if ( params.save_u )
		flowVariables.push_back(&mesh.primitiveVariables.u);
	if ( params.save_v )
		flowVariables.push_back(&mesh.primitiveVariables.v);
	if ( params.save_w )
		flowVariables.push_back(&mesh.primitiveVariables.w);
	if ( params.save_p )
		flowVariables.push_back(&mesh.primitiveVariables.p);
	if ( params.save_T )
		flowVariables.push_back(&mesh.primitiveVariables.T);
	if ( params.save_mu )
		flowVariables.push_back(&mesh.transportProperties.mu);
	if ( params.save_kappa )
		flowVariables.push_back(&mesh.transportProperties.kappa);
	return flowVariables;
}

// Returns a camma-separated string with the names of the flow variables to save.
// For ParaView to read a .csv file it needs the headers on the first line.
string OutputManager::get_csvHeaderString()
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
vector<string> OutputManager::getVariableFileNames()
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
void OutputManager::writePlaneTo_csv(ofstream& outputFile, const Array3D_d* flowVariable)
{
	switch(params.saveNormalAxis)
	{
	case saveNormalAxisEnum::x:
		for(size_t j=0; j<params.NJ; ++j)
		{
			outputFile << (*flowVariable)(params.saveConstantIndex, j, 0);
			for(size_t k=1; k<params.NK; ++k)
				outputFile << ", " << (*flowVariable)(params.saveConstantIndex, j, k);
			outputFile << endl;
		}
		break;
	case saveNormalAxisEnum::y:
		for(size_t i=0; i<params.NI; ++i)
		{
			outputFile << (*flowVariable)(i, params.saveConstantIndex, 0);
			for(size_t k=1; k<params.NK; ++k)
				outputFile << ", " << (*flowVariable)(i, params.saveConstantIndex, k);
			outputFile << endl;
		}
		break;
	case saveNormalAxisEnum::z:
		for(size_t i=0; i<params.NI; ++i)
		{
			outputFile << (*flowVariable)(i, 0, params.saveConstantIndex);
			for(size_t j=1; j<params.NJ; ++j)
				outputFile << ", " << (*flowVariable)(i, j, params.saveConstantIndex);
			outputFile << endl;
		}
		break;
	}
}

// Writes a brief report with progression, timestep size and wall clock time elapsed, etc.
// 'setprecision' is used to control no. of significant digits, default is 6.
void OutputManager::writeStatusReport_toScreen(double t, ulong timeLevel, double dt)
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
	cout << "Wall clock time: " << setprecision(3) << wallClockTimer.getElapsedTime().asSeconds() << setprecision(6) << " sec" << endl << endl;
}

// Write a file with a list of the actual times when the solution was saved. For debugging.
void OutputManager::writeOutputTimes()
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
void OutputManager::writeConvergenceHistoryFiles(const vector<ConservedVariablesScalars>& convergenceHistory)
{
	ofstream normFileRho, normFileRho_u, normFileRho_v, normFileRho_w, normFileRho_E;
	normFileRho.  open( "output/norm_rho.dat"   );
	normFileRho_u.open( "output/norm_rho_u.dat" );
	normFileRho_v.open( "output/norm_rho_v.dat" );
	normFileRho_w.open( "output/norm_rho_w.dat" );
	normFileRho_E.open( "output/norm_rho_E.dat" );
	if ( !normFileRho || !normFileRho_u || !normFileRho_v || !normFileRho_w || !normFileRho_E )
	{
		cout << "Could not open a convergence history file. " << endl
				<< "Norm history was not written to .dat file. You may have to move the 'output' folder to the location of the source files or to the executable itself, depending on whether you run the executable through an IDE or by itself." << endl;
	}
	else
		for (ConservedVariablesScalars norms : convergenceHistory)
		{
			normFileRho   << norms.rho   << endl;
			normFileRho_u << norms.rho_u << endl;
			normFileRho_v << norms.rho_v << endl;
			normFileRho_w << norms.rho_w << endl;
			normFileRho_E << norms.rho_E << endl;
		}
	normFileRho.  close();
	normFileRho_u.close();
	normFileRho_v.close();
	normFileRho_w.close();
	normFileRho_E.close();
}






