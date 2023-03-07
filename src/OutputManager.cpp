/*
 * OutputManager.cpp
 *
 *  Created on: Oct 7, 2022
 *      Author: frederk
 */

#include "OutputManager.h"
#include "ConstantLiteralExpressions.h"
#include "SolutionImporter.h"

// Constructor, taking parameters from config file.
OutputManager::OutputManager(const ConfigSettings& params) :
params(params),
savedSolutions{0},
timeLevelStart{0},
previousProgression{0}
{
}

// Preparation before simulation start. Checking for output folder and deleting old results.
void OutputManager::initialize(long _timeLevelStart)
{
	timeLevelStart = _timeLevelStart;
	bool outputFolderExists = false;
	std::filesystem::path outputPath(StringLiterals::outputFolder);
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
	if(params.continueSimulation)
	{
		SolutionImporter reader(params);
		outputTimes = reader.getSolutionTimes();
		savedSolutions = outputTimes.size();
	}
	else // Start from scratch
	{
		// Delete all files except ParaView state files:
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
}

// Check if anything should be written before sim start.
void OutputManager::processInitialOutput(const Mesh& mesh, double t, long timeLevel)
{
	if ( params.saveIC )
		storeCurrentSolution(mesh, t, timeLevel);
}

// Checks whether it is time to store the solution to disk, or write out a status report to screen and does it if appropriate.
// It will save solution if the next timestep would take the solution past the save-time.
// Thus, in general it saves too early, but the deviation is dt at most.
void OutputManager::processIntermediateOutput(const Mesh& mesh,
											  double t, // <- time
											  long timeLevel,
											  double dt) // <- timestep size
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
			storeCurrentSolution(mesh, t, timeLevel);
	}
	// Write to screen:
	double timeSinceStatusReport = statusReportTimer.getElapsedTime().asSeconds();
	if ( timeSinceStatusReport >= params.statusReportInterval )
	{
		writeStatusReport_toScreen(t, timeLevel, dt);
		statusReportTimer.restart();
	}
}

// If appropriate, process output after simulation has reached its end.
void OutputManager::processFinalOutput(const Mesh& mesh,
									   double t, // <- time
									   long timeLevel,
									   double dt, // <- timestep size
									   const ConservedVariablesVectorGroup& convergenceHistory,
									   const vector<double>& lift,
									   const vector<double>& drag,
									   const vector<double>& separationAngles)
{
	if ( params.saveFinal )
		storeCurrentSolution(mesh, t, timeLevel);
	writeStatusReport_toScreen(t, timeLevel, dt);
	writeOutputTimes();
	writeConvergenceHistoryFiles(convergenceHistory);
	writeIntegralProperties(lift, drag, separationAngles);
}

// Store selected variables from the solution at current time level, to disk, and remember the time.
void OutputManager::storeCurrentSolution(const Mesh& mesh, double t, long timeLevel)
{
	storeCurrentSolution_vtk(mesh, t);
	++savedSolutions;
	outputTimes.push_back(t);
	writeTimeLevel(timeLevel);
}

// Get a string with the first lines we need in the .vtk file.
string OutputManager::getVtkHeader(const Mesh& mesh, double t)
{
	string vtkHeader;
	vtkHeader += "# vtk DataFile Version 4.1\n";
	vtkHeader += "Output at time t = " + to_string(t) + "\n";
	vtkHeader += "ASCII\n";
	vtkHeader += "DATASET STRUCTURED_POINTS\n";
	vtkHeader += "DIMENSIONS " + to_string(mesh.NI) + " " + to_string(mesh.NJ) + " " + to_string(mesh.NK) + "\n";
	vtkHeader += "ORIGIN " + to_string(mesh.positionOffset.x) + " "
						   + to_string(mesh.positionOffset.y) + " "
						   + to_string(mesh.positionOffset.z) + "\n";
	vtkHeader += "SPACING " + to_string(mesh.dx) + " " + to_string(mesh.dy) + " " + to_string(mesh.dz) + "\n";
	vtkHeader += "POINT_DATA " + to_string(mesh.nNodesTotal) + "\n";
	return vtkHeader;
}

// Define a scalar integer flag value to describe the type of each node, and write that to the .vtk file.
void OutputManager::writeVtkNodeFlags(const Mesh& mesh,
									  ofstream& outputFile) // <- OUTPUT, vtk file
{
	outputFile << "SCALARS Node_Flag int 1\n";
	outputFile << "LOOKUP_TABLE default\n";
	for (int k{0}; k<mesh.NK; ++k)
		for (int j{0}; j<mesh.NJ; ++j)
			for (int i{0}; i<mesh.NI; ++i)
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

// Get pointers to the arrays containing scalar flow variables (non-vectors).
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

// Get the names of the scalar (non-vector) flow variables as strings without spaces.
vector<string> OutputManager::getScalarVariableNames()
{
	vector<string> variableNames;
	if ( params.saveDensity )
		variableNames.push_back(StringLiterals::densityString);
	if ( params.saveEnergy )
		variableNames.push_back(StringLiterals::energyString);
	if ( params.savePressure )
		variableNames.push_back(StringLiterals::pressureString);
	if ( params.saveTemperature )
		variableNames.push_back(StringLiterals::temperatureString);
	if ( params.saveViscosity )
		variableNames.push_back(StringLiterals::viscosityString);
	if ( params.saveThermalCond )
		variableNames.push_back(StringLiterals::thermalCondString);
	return variableNames;
}

// Write the current values of the scalar (non-vector) flow variables to specified vtk file.
void OutputManager::writeVtkScalarData(const Mesh& mesh,
									   ofstream& outputFile) // <- OUTPUT, vtk file
{
	vector<const Array3D_d*> scalarFlowVariables = getScalarVariablePointers(mesh);
	vector<string> variableNames = getScalarVariableNames();
	int counter = 0;
	for (const Array3D_d* flowVariable : scalarFlowVariables)
	{
		outputFile << "SCALARS " << variableNames.at(counter) << " double 1\n";
		outputFile << "LOOKUP_TABLE default\n";
		for (int k{0}; k<mesh.NK; ++k)
			for (int j{0}; j<mesh.NJ; ++j)
				for (int i{0}; i<mesh.NI; ++i)
					outputFile << flowVariable->at(i,j,k) << "\n";
		++counter;
	}
}

// Get pointers to the arrays containing vector flow variables. The std::array with fixed size 3
// represents the 3D vectors' three components.
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

// Get the names of the vector flow variables as strings without spaces.
vector<string> OutputManager::getVectorVariableNames()
{
	vector<string> variableNames;
	if ( params.saveMomentum )
		variableNames.push_back(StringLiterals::momentumString);
	if ( params.saveVelocity )
		variableNames.push_back(StringLiterals::velocityString);
	return variableNames;
}

// Write the current values of the vector flow variables to specified vtk file.
void OutputManager::writeVtkVectorData(const Mesh& mesh,
									   ofstream& outputFile) // <- OUTPUT, vtk file
{
	// The std::array with fixed size 3 represents the 3D vectors' three components:
	vector<std::array<const Array3D_d*, 3>> vectorFlowVariables = getVectorVariablePointers(mesh);
	vector<string> variableNames = getVectorVariableNames();
	int counter = 0;
	for (std::array<const Array3D_d*, 3> flowVariableVector : vectorFlowVariables)
	{
		outputFile << "VECTORS " << variableNames.at(counter) << " double\n";
		for (int k{0}; k<mesh.NK; ++k)
			for (int j{0}; j<mesh.NJ; ++j)
				for (int i{0}; i<mesh.NI; ++i)
					outputFile << flowVariableVector[0]->at(i,j,k) << " "	// <- x-component
							   << flowVariableVector[1]->at(i,j,k) << " "	// <- y-component
							   << flowVariableVector[2]->at(i,j,k) << "\n";	// <- z-component
		++counter;
	}
}

// Writes a .vtk file with the flow variables specified in the config file at the current time level.
// The file is formatted to fit how ParaView wants to import it. If there are multiple files, they
// are numbered 0, 1, ... so ParaView will import it as a time series.
void OutputManager::storeCurrentSolution_vtk(const Mesh& mesh, double t)
{
	ofstream outputFile;
	string filename = string(StringLiterals::outputFolder) + StringLiterals::solutionFileNameBase + to_string(savedSolutions);
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
// Progression is computed based on the stopping criterion (end time or time level).
// 'setprecision' is used to control no. of significant digits, default is 6.
void OutputManager::writeStatusReport_toScreen(double t,	// <- time
											   long timeLevel,
											   double dt)	// <- timestep size
{
	cout << "Simulated time: t = " << t;
	double progressPercentage{0};
	if (params.stopCriterion == StopCriterionEnum::end_time)
	{
		progressPercentage = t / params.t_end * 100.;
		cout << " , t_end = " << params.t_end << " ( " << setprecision(4) << progressPercentage << setprecision(6) << " % )" << endl;
	}
	else
		cout << endl;

	cout << "Time level: n = " << timeLevel;
	if (params.stopCriterion == StopCriterionEnum::timesteps)
	{
		progressPercentage = (double)(timeLevel-timeLevelStart) / (double)params.stopTimeLevel * 100.;  // Cast to double to avoid integer division.
		cout << " , n_max = " << timeLevelStart+params.stopTimeLevel << " ( " << setprecision(4) << progressPercentage << setprecision(6) << " % )" << endl;
	}
	else
		cout << endl;
	cout << "Timestep size: dt = " << dt << endl;
	cout << "Wall clock time: " << getTimeString( wallClockTimer.getElapsedTime().asSeconds() ) << endl;
	double newProgress = progressPercentage - previousProgression;
	double timeSinceReport = statusReportTimer.getElapsedTime().asSeconds();
	double estimatedTimeLeft = timeSinceReport / newProgress * (100-progressPercentage);
	cout << "Estimated time to finish: " + getTimeString(estimatedTimeLeft) << endl << endl;
	previousProgression = progressPercentage;
}

// Write a file with a list of the actual times when the solution was saved.
void OutputManager::writeOutputTimes()
{
	ofstream timeFile;
	string filename = string(StringLiterals::outputFolder) + StringLiterals::timesFile;
	timeFile.open( filename );
	if ( !timeFile )
	{
		cout << "Could not open file:  " + filename << endl
		     << "Output times were not written to .dat file. You may have to move the 'output' folder to the location of the source files or to the executable itself, depending on whether you run the executable through an IDE or by itself." << endl;
		return;
	}
	for (double time : outputTimes)
		timeFile << setprecision(8) << time << endl;
	timeFile.close();
}

// Write files with history of lift and drag. Also separation angles in final solution.
void OutputManager::writeIntegralProperties(const vector<double>& liftHistory,
											const vector<double>& dragHistory,
											const vector<double>& separationAngles)
{
	vector<string> names = {"lift", "drag", "separation_angles"};
	vector<vector<double>> dataSets = {liftHistory, dragHistory, separationAngles};
	for (int i : {0,1,2})
	{
		ofstream outputFile;
		string filename = StringLiterals::outputFolder + names.at(i) + ".dat";
		outputFile.open( filename );
		if ( !outputFile)
		{
			cout << "Could not open file:  " + filename << endl
				 << names.at(i) + " was not written to .dat file. You may have to move the 'output' folder to the location of the source files or to the executable itself, depending on whether you run the executable through an IDE or by itself." << endl;
			return;
		}
		for (double scalar : dataSets.at(i))
			outputFile << scalar << endl;
		outputFile.close();
	}
}

void OutputManager::writeTimeLevel(int timeLevel)
{
	ofstream outputFile( string(StringLiterals::outputFolder) + StringLiterals::timeLevelFile );
	if(!outputFile)
		std::cerr << "Could not create file with number of time steps." << endl;
	else
		outputFile << timeLevel;
	outputFile.close();
}

string OutputManager::getTimeString(double secondsInput)
{
	double hours = floor(secondsInput / 3600);
	double minutes = floor( fmod(secondsInput, 3600) / 60 );
	double seconds = fmod( fmod(secondsInput, 3600), 60 );
	std::stringstream timeSS;
	if(secondsInput > 60)
		timeSS << hours << ":" << setfill('0') << setw(2) <<  minutes << ":" << setfill('0') << setw(2) << round(seconds);
	else
		timeSS << setprecision(3) << secondsInput << " sec";
	return timeSS.str();
}

// Write files with lists of the norm of change for the conserved variables.
void OutputManager::writeConvergenceHistoryFiles(const ConservedVariablesVectorGroup& convergenceHistory)
{
	vector< vector<double> const* > normHistoriesToSave;
	vector<string> fileNames;
	if(params.saveConvergenceHistory != ConvHistoryEnum::none)
	{
		normHistoriesToSave.push_back(&convergenceHistory.rho);
		fileNames.push_back(StringLiterals::normRhoFile);
	}
	if(params.saveConvergenceHistory == ConvHistoryEnum::all)
	{
		normHistoriesToSave.push_back(&convergenceHistory.rho_u);
		fileNames.push_back(StringLiterals::normRhoUFile);
		normHistoriesToSave.push_back(&convergenceHistory.rho_v);
		fileNames.push_back(StringLiterals::normRhoVFile);
		normHistoriesToSave.push_back(&convergenceHistory.rho_w);
		fileNames.push_back(StringLiterals::normRhoWFile);
		normHistoriesToSave.push_back(&convergenceHistory.rho_E);
		fileNames.push_back(StringLiterals::normRhoEFile);
	}
	for(int varIndex{0}; varIndex<fileNames.size(); ++varIndex)
	{
		ofstream normFile;
		normFile.open( StringLiterals::outputFolder + fileNames.at(varIndex) );
		if ( !normFile )
		{
			cout << "Could not open convergence history file for " + fileNames.at(varIndex) + ". " << endl
					<< "Norm history was not written to .dat file. You may have to move the 'output' folder to the location of the source files or to the executable itself, depending on whether you run the executable through an IDE or by itself." << endl;
		}
		else
			for (double norm : *normHistoriesToSave.at(varIndex))
				normFile << norm << endl;

		normFile.close();
	}
}






