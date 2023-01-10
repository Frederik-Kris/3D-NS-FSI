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
		storeCurrentSolution(mesh, t);
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
			storeCurrentSolution(mesh, t);
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
		storeCurrentSolution(mesh, t);
	writeStatusReport_toScreen(t, timeLevel, dt);
	writeOutputTimes();
	writeConvergenceHistoryFiles(convergenceHistory);
}

// Store selected variables from the solution at current time level, using the format(s) specified
void OutputManager::storeCurrentSolution(const Mesh& mesh, double t)
{
	storeCurrentSolution_vtk(mesh, t);
	++savedSolutions;
	outputTimes.push_back(t);
}

string OutputManager::getVtkHeader(const Mesh& mesh, double t)
{
	string vtkHeader;
	vtkHeader += "# vtk DataFile Version 5.1\n";
	vtkHeader += "Output at time t = " + to_string(t) + "\n";
	vtkHeader += "ASCII\n";
	vtkHeader += "DATASET STRUCTURED_POINTS\n";
	vtkHeader += "DIMENSIONS " + to_string(mesh.NI) + " " + to_string(mesh.NJ) + " " + to_string(mesh.NK) + "\n";
	vtkHeader += "ORIGIN " + to_string(-mesh.positionOffset.x * mesh.dx) + " "
						   + to_string(-mesh.positionOffset.y * mesh.dy) + " "
						   + to_string(-mesh.positionOffset.z * mesh.dz) + "\n";
	vtkHeader += "SPACING " + to_string(mesh.dx) + " " + to_string(mesh.dy) + " " + to_string(mesh.dz) + "\n";
	vtkHeader += "POINT_DATA " + to_string(mesh.nNodesTotal) + "\n";
	return vtkHeader;
}

void OutputManager::writeVtkNodeFlags(const Mesh& mesh, ofstream& outputFile)
{
	outputFile << "SCALARS Node_Flag int 1\n";
	outputFile << "LOOKUP_TABLE default\n";
	for (size_t k{0}; k<mesh.NK; ++k)
		for (size_t j{0}; j<mesh.NJ; ++j)
			for (size_t i{0}; i<mesh.NI; ++i)
			{
				int flagValue = 0;
				if(mesh.nodeType(i,j,k) == NodeTypeEnum::FluidInterior)
					flagValue = 0;
				else if(mesh.nodeType(i,j,k) == NodeTypeEnum::FluidEdge)
					flagValue = 1;
				else if(mesh.nodeType(i,j,k) == NodeTypeEnum::FluidGhost)
					flagValue = 2;
				else if(mesh.nodeType(i,j,k) == NodeTypeEnum::SolidGhost)
					flagValue = 3;
				else if(mesh.nodeType(i,j,k) == NodeTypeEnum::SolidInactive)
					flagValue = 4;
				else
					throw std::logic_error("Unexpected enum value.");
				outputFile << flagValue << "\n";
			}
}

vector<const Array3D_d*> OutputManager::getScalarVariablePointers(const Mesh& mesh)
{
	vector<const Array3D_d*> scalarFlowVariables;
	if ( params.saveDensity )
		scalarFlowVariables.push_back(&mesh.conservedVariables.rho);
	if ( params.saveEnergy )
		scalarFlowVariables.push_back(&mesh.conservedVariables.rho_E);
	if ( params.savePressure )
		scalarFlowVariables.push_back(&mesh.primitiveVariables.p);
	if ( params.saveTemperature )
		scalarFlowVariables.push_back(&mesh.primitiveVariables.T);
	if ( params.saveViscosity )
		scalarFlowVariables.push_back(&mesh.transportProperties.mu);
	if ( params.saveThermalCond )
		scalarFlowVariables.push_back(&mesh.transportProperties.kappa);
	return scalarFlowVariables;
}

vector<string> OutputManager::getScalarVariableNames()
{
	vector<string> variableNames;
	if ( params.saveDensity )
		variableNames.push_back("Density");
	if ( params.saveEnergy )
		variableNames.push_back("Total_Specific_Energy_per_Volume");
	if ( params.savePressure )
		variableNames.push_back("Pressure");
	if ( params.saveTemperature )
		variableNames.push_back("Temperature");
	if ( params.saveViscosity )
		variableNames.push_back("Dynamic_Viscosity");
	if ( params.saveThermalCond )
		variableNames.push_back("Thermal_Conductivity");
	return variableNames;
}

void OutputManager::writeVtkScalarData(const Mesh& mesh, ofstream& outputFile)
{
	vector<const Array3D_d*> scalarFlowVariables = getScalarVariablePointers(mesh);
	vector<string> variableNames = getScalarVariableNames();
	uint counter = 0;
	for (const Array3D_d* flowVariable : scalarFlowVariables)
	{
		outputFile << "SCALARS " << variableNames.at(counter) << " double 1\n";
		outputFile << "LOOKUP_TABLE default\n";
		for (size_t k{0}; k<mesh.NK; ++k)
			for (size_t j{0}; j<mesh.NJ; ++j)
				for (size_t i{0}; i<mesh.NI; ++i)
					outputFile << flowVariable->at(i,j,k) << "\n";
		++counter;
	}
}

vector<std::array<const Array3D_d*, 3>> OutputManager::getVectorVariablePointers(const Mesh& mesh)
{
	vector<std::array<const Array3D_d*, 3>> vectorFlowVariables;
	if ( params.saveMomentum )
	{
		vectorFlowVariables.push_back( std::array<const Array3D_d*, 3>({&mesh.conservedVariables.rho_u,
																		&mesh.conservedVariables.rho_v,
																		&mesh.conservedVariables.rho_w}) );
	}
	if ( params.saveVelocity )
	{
		vectorFlowVariables.push_back( std::array<const Array3D_d*, 3>({&mesh.primitiveVariables.u,
																		&mesh.primitiveVariables.v,
																		&mesh.primitiveVariables.w}) );
	}
	return vectorFlowVariables;
}

vector<string> OutputManager::getVectorVariableNames()
{
	vector<string> variableNames;
	if ( params.saveMomentum )
		variableNames.push_back("Momentum_Density");
	if ( params.saveVelocity )
		variableNames.push_back("Velocity");
	return variableNames;
}

void OutputManager::writeVtkVectorData(const Mesh& mesh, ofstream& outputFile)
{
	vector<std::array<const Array3D_d*, 3>> vectorFlowVariables = getVectorVariablePointers(mesh);
	vector<string> variableNames = getVectorVariableNames();
	uint counter = 0;
	for (std::array<const Array3D_d*, 3> flowVariableVector : vectorFlowVariables)
	{
		outputFile << "VECTORS " << variableNames.at(counter) << " double\n";
		for (size_t k{0}; k<mesh.NK; ++k)
			for (size_t j{0}; j<mesh.NJ; ++j)
				for (size_t i{0}; i<mesh.NI; ++i)
					outputFile << flowVariableVector[0]->at(i,j,k) << " "
							   << flowVariableVector[1]->at(i,j,k) << " "
							   << flowVariableVector[2]->at(i,j,k) << "\n";
		++counter;
	}
}

// Writes a .csv file with the specified flow variables at the current time level.
// The file is formatted to fit how ParaView wants to import it, with coordinates in the leftmost column.
// Only the flow variables selected in the ConfigFile are saved.
void OutputManager::storeCurrentSolution_vtk(const Mesh& mesh, double t)
{
	ofstream outputFile;
	string filename = "output/out.vtk." + to_string(savedSolutions);
	outputFile.open( filename );
	if ( !outputFile )
	{
		cout << "Could not open file:  " + filename << endl
		     << "Solution was not saved. You may have to move the 'output' folder to the location of the source files or to the executable itself, depending on whether you run the executable through an IDE or by itself." << endl;
		return;
	}
	outputFile << getVtkHeader(mesh, t);
	writeVtkNodeFlags (mesh, outputFile);
	writeVtkScalarData(mesh, outputFile);
	writeVtkVectorData(mesh, outputFile);
	outputFile.close();
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






