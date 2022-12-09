/*
 * OutputManager.cpp
 *
 *  Created on: Oct 7, 2022
 *      Author: frederk
 */

#include "OutputManager.h"

OutputManager::OutputManager(const ConfigSettings& params) :
params(params),
savedSolutions{0},
vorticity(params.NI, params.NJ, params.NK)
{
}

void OutputManager::initialize()
{
	bool outputFolderExists = false;
	std::filesystem::path outputPath("./output");
	if(std::filesystem::exists(outputPath))
		if(std::filesystem::is_directory(outputPath))
			outputFolderExists = true;

	if(!outputFolderExists)
	{
		bool folderWasCreated = std::filesystem::create_directory(outputPath);
		if(!folderWasCreated)
		{
			std::cerr << "Output folder could not be created.\n";
			return;
		}
	}
	for(std::filesystem::directory_entry entry : std::filesystem::recursive_directory_iterator(outputPath))
		if( entry.path().generic_string().find(".pvsm") == string::npos ) // if path does not contain ".."
			try
				{ std::filesystem::remove(entry.path()); }
			catch (std::exception& ex)
			{
				std::cerr << "Could not remove " << entry.path() << endl;
				std::cerr << ex.what() << endl;
			}
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
								const ConservedVariablesVectorGroup& convergenceHistory)
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

void OutputManager::writeValuesFromIndices_csv_paraview(const Mesh& mesh, ofstream& outputFile, const vector<size_t>& nodesToWrite, const vector<Array3D_d const*>& flowVariables)
{
	Vector3_u nMeshNodes(mesh.NI, mesh.NJ, mesh.NK);
	Vector3_d gridSpacing(mesh.dx, mesh.dy, mesh.dz);
	for (size_t index1D : nodesToWrite)
	{
		outputFile << endl;
		Vector3_d position = getNodePosition( getIndices3D(index1D, nMeshNodes), gridSpacing, mesh.positionOffset );
		outputFile << position.x << ", " << position.y << ", " << position.z;
		for (const Array3D_d* flowVar : flowVariables)
			outputFile << ", " << flowVar->at(index1D);
	}
}

// Writes a .csv file with the specified flow variables at the current time level.
// The file is formatted to fit how ParaView wants to import it, with coordinates in the leftmost column.
// Only the flow variables selected in the ConfigFile are saved.
void OutputManager::storeCurrentSolution_csv_paraview(const Mesh& mesh)
{
	vector<Array3D_d const*> flowVariables = getPlotVariables(mesh);
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
	writeValuesFromIndices_csv_paraview(mesh, outputFile, mesh.indexByType.fluidInterior, flowVariables);
	writeValuesFromIndices_csv_paraview(mesh, outputFile, mesh.indexByType.fluidEdge,	  flowVariables);
//	writeValuesFromIndices_csv_paraview(mesh, outputFile, mesh.indexByType.ghost,		  flowVariables);

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

void OutputManager::computeVorticity(const Mesh& mesh, const AxisOrientationEnum axis)
{
	for (size_t index1D : mesh.indexByType.fluidInterior)
	{
		Vector3_u indices = getIndices3D(index1D, Vector3_u(mesh.NI, mesh.NJ, mesh.NK));
		if (axis == AxisOrientationEnum::x)
		{
			vorticity.x(index1D) = ( mesh.primitiveVariables.w(indices.i, indices.j+1, indices.k)
								   - mesh.primitiveVariables.w(indices.i, indices.j-1, indices.k) ) / (2*mesh.dy)
								 - ( mesh.primitiveVariables.v(indices.i, indices.j, indices.k+1)
								   - mesh.primitiveVariables.v(indices.i, indices.j, indices.k-1) ) / (2*mesh.dz);
		}
		if (axis == AxisOrientationEnum::y)
		{
			vorticity.y(index1D) = ( mesh.primitiveVariables.u(indices.i, indices.j, indices.k+1)
								   - mesh.primitiveVariables.u(indices.i, indices.j, indices.k-1) ) / (2*mesh.dz)
								 - ( mesh.primitiveVariables.w(indices.i+1, indices.j, indices.k)
								   - mesh.primitiveVariables.w(indices.i-1, indices.j, indices.k) ) / (2*mesh.dx);
		}
		if (axis == AxisOrientationEnum::z)
		{
			vorticity.z(index1D) = ( mesh.primitiveVariables.v(indices.i+1, indices.j, indices.k)
								   - mesh.primitiveVariables.v(indices.i-1, indices.j, indices.k) ) / (2*mesh.dx)
								 - ( mesh.primitiveVariables.u(indices.i, indices.j+1, indices.k)
								   - mesh.primitiveVariables.u(indices.i, indices.j-1, indices.k) ) / (2*mesh.dy);
		}
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
	if ( params.save_vorticity_x )
	{
		computeVorticity(mesh, AxisOrientationEnum::x);
		flowVariables.push_back(&vorticity.x);
	}
	if ( params.save_vorticity_y )
	{
		computeVorticity(mesh, AxisOrientationEnum::y);
		flowVariables.push_back(&vorticity.y);
	}
	if ( params.save_vorticity_z )
	{
		computeVorticity(mesh, AxisOrientationEnum::z);
		flowVariables.push_back(&vorticity.z);
	}
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
	if ( params.save_vorticity_x )
		headers += ", x-vorticity";
	if ( params.save_vorticity_y )
		headers += ", y-vorticity";
	if ( params.save_vorticity_z )
		headers += ", z-vorticity";
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
	if ( params.save_vorticity_x )
		variableNames.push_back("vorticity_x");
	if ( params.save_vorticity_y )
		variableNames.push_back("vorticity_y");
	if ( params.save_vorticity_z )
		variableNames.push_back("vorticity_z");
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
void OutputManager::writeConvergenceHistoryFiles(const ConservedVariablesVectorGroup& convergenceHistory)
{
	vector< vector<double> const* > normHistoriesToSave;
	vector<string> variableNames;
	if(params.saveConvergenceHistory != ConvHistoryEnum::none)
	{
		normHistoriesToSave.push_back(&convergenceHistory.rho);
		variableNames.push_back("rho");
	}
	if(params.saveConvergenceHistory == ConvHistoryEnum::all)
	{
		normHistoriesToSave.push_back(&convergenceHistory.rho_u);
		variableNames.push_back("rho_u");
		normHistoriesToSave.push_back(&convergenceHistory.rho_v);
		variableNames.push_back("rho_v");
		normHistoriesToSave.push_back(&convergenceHistory.rho_w);
		variableNames.push_back("rho_w");
		normHistoriesToSave.push_back(&convergenceHistory.rho_E);
		variableNames.push_back("rho_E");
	}
	for(size_t varIndex{0}; varIndex<variableNames.size(); ++varIndex)
	{
		ofstream normFile;
		normFile.open( "output/norm_" + variableNames.at(varIndex) + ".dat"   );
		if ( !normFile )
		{
			cout << "Could not open convergence history file for " + variableNames.at(varIndex) + ". " << endl
					<< "Norm history was not written to .dat file. You may have to move the 'output' folder to the location of the source files or to the executable itself, depending on whether you run the executable through an IDE or by itself." << endl;
		}
		else
			for (double norm : *normHistoriesToSave.at(varIndex))
				normFile << norm << endl;

		normFile.close();
	}
}






