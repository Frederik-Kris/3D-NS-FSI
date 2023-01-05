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
	if(params.saveForParaview)
		storeCurrentSolution_vtk(mesh, t);
	if (params.saveForMatlab)
		storeCurrentSolutionPlane_csv(mesh);
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

// Writes multiple .csv files with the specified flow variables at the current time level.
// The file is formatted as a table or matrix, with the values from a specified 2D plane. Matlab
// can read this file to a matrix. Only the flow variables selected in the ConfigFile are saved.
void OutputManager::storeCurrentSolutionPlane_csv(const Mesh& mesh)
{
	vector<const Array3D_d*> flowVariables = getFlowVariablePointers_csv(mesh);		// Get pointers to the arrays with data to save
	vector<string> variableFileNames = getVariableCsvFileNames();	// Get a vector with the names of the variables

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
		writePlaneToCsv(outputFile, flowVariables.at(i));
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
vector<const Array3D_d*> OutputManager::getFlowVariablePointers_csv(const Mesh& mesh)
{
	vector<const Array3D_d*> flowVariables;
	if ( params.saveDensity )
		flowVariables.push_back(&mesh.conservedVariables.rho);
	if ( params.saveMomentum )
	{
		flowVariables.push_back(&mesh.conservedVariables.rho_u);
		flowVariables.push_back(&mesh.conservedVariables.rho_v);
		flowVariables.push_back(&mesh.conservedVariables.rho_w);
	}
	if ( params.saveEnergy )
		flowVariables.push_back(&mesh.conservedVariables.rho_E);
	if ( params.saveVelocity )
	{
		flowVariables.push_back(&mesh.primitiveVariables.u);
		flowVariables.push_back(&mesh.primitiveVariables.v);
		flowVariables.push_back(&mesh.primitiveVariables.w);
	}
	if ( params.savePressure )
		flowVariables.push_back(&mesh.primitiveVariables.p);
	if ( params.saveTemperature )
		flowVariables.push_back(&mesh.primitiveVariables.T);
	if ( params.saveViscosity )
		flowVariables.push_back(&mesh.transportProperties.mu);
	if ( params.saveThermalCond )
		flowVariables.push_back(&mesh.transportProperties.kappa);
	if ( params.saveVorticity )
	{
		computeVorticity(mesh, AxisOrientationEnum::x);
		flowVariables.push_back(&vorticity.x);
		computeVorticity(mesh, AxisOrientationEnum::y);
		flowVariables.push_back(&vorticity.y);
		computeVorticity(mesh, AxisOrientationEnum::z);
		flowVariables.push_back(&vorticity.z);
	}
	return flowVariables;
}

// Returns a camma-separated string with the names of the flow variables to save.
// For ParaView to read a .csv file it needs the headers on the first line.
string OutputManager::getCsvHeaderString()
{
	string headers = "x, y, z";
	if ( params.saveDensity )
		headers += ", Density";
	if ( params.saveMomentum )
	{
		headers += ", x-momentum";
		headers += ", y-momentum";
		headers += ", z-momentum";
	}
	if ( params.saveEnergy )
		headers += ", Total energy";
	if ( params.saveVelocity )
	{
		headers += ", Velocity comp u";
		headers += ", Velocity comp v";
		headers += ", Velocity comp w";
	}
	if ( params.savePressure )
		headers += ", Pressure";
	if ( params.saveTemperature )
		headers += ", Temperature";
	if ( params.saveViscosity )
		headers += ", Dynamic viscosity";
	if ( params.saveThermalCond )
		headers += ", Thermal conductivity";
	if ( params.saveVorticity )
	{
		headers += ", x-vorticity";
		headers += ", y-vorticity";
		headers += ", z-vorticity";
	}
	return headers;
}

// Returns a vector of strings, which are the names of the variables to save.
vector<string> OutputManager::getVariableCsvFileNames()
{
	vector<string> variableNames;
	if ( params.saveDensity )
		variableNames.push_back("rho");
	if ( params.saveMomentum )
	{
		variableNames.push_back("rho_u");
		variableNames.push_back("rho_v");
		variableNames.push_back("rho_w");
	}
	if ( params.saveEnergy )
		variableNames.push_back("E");
	if ( params.saveVelocity )
	{
		variableNames.push_back("u");
		variableNames.push_back("v");
		variableNames.push_back("w");
	}
	if ( params.savePressure )
		variableNames.push_back("p");
	if ( params.saveTemperature )
		variableNames.push_back("T");
	if ( params.saveViscosity )
		variableNames.push_back("mu");
	if ( params.saveThermalCond )
		variableNames.push_back("kappa");
	if ( params.saveVorticity )
	{
		variableNames.push_back("vorticity_x");
		variableNames.push_back("vorticity_y");
		variableNames.push_back("vorticity_z");
	}
	return variableNames;
}

// Write the values from a plane, defined in ConfigFile, of one flow variable. It's written like a comma separated table into 'outputFile'.
void OutputManager::writePlaneToCsv(ofstream& outputFile, const Array3D_d* flowVariable)
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






