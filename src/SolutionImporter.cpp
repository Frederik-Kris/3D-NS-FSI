/*
 * SolutionImporter.cpp
 *
 *  Created on: Feb 8, 2023
 *      Author: frederk
 */

#include "SolutionImporter.h"
#include "ConstantLiteralExpressions.h"
#include <filesystem>
#include <vtk-7.1/vtkStructuredPoints.h>
#include <vtk-7.1/vtkStructuredPointsReader.h>
#include <vtk-7.1/vtkPointData.h>
#include <vtk-7.1/vtkDataArray.h>
#include <vtk-7.1/vtkSmartPointer.h>
#include <vtk-7.1/vtkImageAlgorithm.h>

// Constructor. Finds correct vtk-file and reads it. Keeps pointer to the data in the pointData member.
// Also checks which of the flow variables that are available in the file.
SolutionImporter::SolutionImporter(const ConfigSettings& params)
: params{params},
  pointData{NULL}
{}

// Search the output folder for solution files. Return the last in the time series.
string SolutionImporter::findLatestSolutionFile()
{
	vector<int> solutionFileNumbers;
	std::filesystem::path outputPath(StringLiterals::outputFolder);
	// Loop through entries in folder "./output"
	for(std::filesystem::directory_entry entry : std::filesystem::recursive_directory_iterator(outputPath))
		if( entry.path().generic_string().find(StringLiterals::solutionFileNameBase) != string::npos ) // if entry path contains ".."
			{
				int lastDot = entry.path().generic_string().find_last_of('.');	// Get position of last dot
				std::stringstream ss;
				ss << entry.path().generic_string().substr(lastDot+1);	// Read the part after last dot (the no.) into stream
				int fileNumber;
				ss >> fileNumber;
				solutionFileNumbers.push_back(fileNumber);
			}
	if( !solutionFileNumbers.empty() )
	{
		std::stringstream ss;
		ss << StringLiterals::outputFolder
		   << StringLiterals::solutionFileNameBase
		   << *std::max_element( solutionFileNumbers.begin(), solutionFileNumbers.end() );
		return ss.str();
	}
	else
		throw std::runtime_error("Did not find any output vtk files.");
}

void SolutionImporter::loadVtkFile()
{
	string filename = findLatestSolutionFile();
	reader = vtkSmartPointer<vtkStructuredPointsReader>::New();
	reader->SetFileName( filename.c_str() );
	reader->ReadAllScalarsOn();
	reader->ReadAllVectorsOn();
	reader->Update();
	vtkStructuredPoints* readerOutput = reader->GetOutput();
	pointData = readerOutput->GetPointData();
	vtkFileLoaded = true;
	checkAvailableFlowVars();
}

// Check which flow variables that were saved in the vtk-file.
void SolutionImporter::checkAvailableFlowVars()
{
	if( !vtkFileLoaded )
		throw std::logic_error("Cannot check for flow variables before vtk-file is loaded.");

	for(int i{0}; i<pointData->GetNumberOfArrays(); ++i)
	{
		string arrayName = pointData->GetArrayName(i);
		if(arrayName == StringLiterals::densityString)
			densityExists = true;
		else if(arrayName == StringLiterals::momentumString)
			momentumExists = true;
		else if(arrayName == StringLiterals::energyString)
			energyExists = true;
		else if(arrayName == StringLiterals::velocityString)
			velocityExists = true;
		else if(arrayName == StringLiterals::pressureString)
			pressureExists = true;
		else if(arrayName == StringLiterals::temperatureString)
			temperatureExists = true;
	}
}

// Check whether the vtk-file contains enough flow variable to derive the rest.
// Actually it is possible with many more combinations than this, but I just enable the most likely ones.
bool SolutionImporter::allFlowVarsDerivable()
{
	if( !vtkFileLoaded )
		throw std::logic_error("Cannot check for flow variables before vtk-file is loaded.");

	bool densityDerivable = densityExists || (pressureExists && temperatureExists);
	bool momentumDerivable = momentumExists || (velocityExists && densityDerivable);
	bool energyDerivable = energyExists || (pressureExists && densityDerivable && momentumDerivable);
	return densityDerivable && momentumDerivable && energyDerivable;
}

// Compute all conserved and primitive variables and transport properties in one node.
// NB! 1D index in vtk, "i" must be ordered oppositely from the mesh arrays (first index increases most frequent, Fortran style).
void SolutionImporter::computeAllFlowVars(int i,
										  ConservedVariablesScalars& conservedVars,
										  PrimitiveVariablesScalars& primitiveVars,
										  TransportPropertiesScalars& transportProps)
{
	if( !vtkFileLoaded )
		throw std::logic_error("Cannot import flow variables before vtk-file is loaded.");

	double rho=0, rhoU=0, rhoV=0, rhoW=0, rhoE=0;
	double u=0, v=0, w=0, p=0, T=0;
	if(densityExists)
	{
		rho = pointData->GetArray(StringLiterals::densityString)->GetComponent(i, 0);
		if(pressureExists)
			p = pointData->GetArray(StringLiterals::pressureString)->GetComponent(i, 0);
		if(temperatureExists)
			T = pointData->GetArray(StringLiterals::temperatureString)->GetComponent(i, 0);
	}
	else // p and T exist
	{
		p = pointData->GetArray(StringLiterals::pressureString)->GetComponent(i, 0);
		T = pointData->GetArray(StringLiterals::temperatureString)->GetComponent(i, 0);
		rho = ( params.Gamma * p - T ) / ( 1 + T );
	}
	if(!temperatureExists)
		T = ( params.Gamma * p - rho ) / ( 1 + rho );
	if(momentumExists)
	{
		rhoU = pointData->GetArray(StringLiterals::momentumString)->GetComponent(i, 0);
		rhoV = pointData->GetArray(StringLiterals::momentumString)->GetComponent(i, 1);
		rhoW = pointData->GetArray(StringLiterals::momentumString)->GetComponent(i, 2);
		if(!velocityExists)
		{
			u = rhoU / (1+rho);
			v = rhoV / (1+rho);
			w = rhoW / (1+rho);
		}
	}
	else // velocity exists
	{
		u = pointData->GetArray(StringLiterals::velocityString)->GetComponent(i, 0);
		v = pointData->GetArray(StringLiterals::velocityString)->GetComponent(i, 1);
		w = pointData->GetArray(StringLiterals::velocityString)->GetComponent(i, 2);
		rhoU = (1+rho) * u;
		rhoV = (1+rho) * v;
		rhoW = (1+rho) * w;
	}
	if(energyExists)
	{
		rhoE = pointData->GetArray(StringLiterals::energyString)->GetComponent(i, 0);
		if(!pressureExists)
			p = ( params.Gamma - 1 )*( rhoE - (1 + rho)/2 * ( u*u + v*v + w*w ));
	}
	else // rho and p exist
		rhoE = p / ( params.Gamma - 1 ) + (1 + rho)/2 * ( u*u + v*v + w*w );
	conservedVars.rho = rho;
	conservedVars.rho_u = rhoU;
	conservedVars.rho_v = rhoV;
	conservedVars.rho_w = rhoW;
	conservedVars.rho_E = rhoE;
	primitiveVars.u = u;
	primitiveVars.v = v;
	primitiveVars.w = w;
	primitiveVars.p = p;
	primitiveVars.T = T;
	transportProps = deriveTransportProperties(primitiveVars, params);
}

void SolutionImporter::importNormHistories(ConservedVariablesVectorGroup &normHistory, long unsigned timeLevel)
{
	vector<string> filenames = {StringLiterals::normRhoFile,
								StringLiterals::normRhoUFile,
								StringLiterals::normRhoVFile,
								StringLiterals::normRhoWFile,
								StringLiterals::normRhoEFile };
	vector<vector<double>*> normVectors = { &normHistory.rho,
											&normHistory.rho_u,
											&normHistory.rho_v,
											&normHistory.rho_w,
											&normHistory.rho_E };
	for(size_t i{0}; i<filenames.size(); ++i)
	{
		std::filesystem::path filePath( StringLiterals::outputFolder + filenames.at(i) );
		if( std::filesystem::exists(filePath) )
		{
			std::ifstream normFile(filePath);
			if( !normFile.is_open() )
				throw std::runtime_error( "Could not open file: " + filenames.at(i) );
			double norm;
			while( normFile >> norm )
			{
				normVectors.at(i)->push_back(norm);
				if( normVectors.at(i)->size() == timeLevel )
					break;	// If one or more time steps were computed after the result
			}				// was saved, we don't want the norms from those time steps.
		}
	}
}

vector<double> SolutionImporter::getSolutionTimes()
{
	std::filesystem::path timeFilePath( string(StringLiterals::outputFolder) + StringLiterals::timesFile );
	std::ifstream timeFile(timeFilePath);
	if( !timeFile.is_open() )
		throw std::runtime_error( "Could not open file: " + timeFilePath.generic_string() );
	vector<double> solutionTimes;
	double t;
	while( timeFile >> t )
	{
		solutionTimes.push_back(t);
	}
	return solutionTimes;
}

long SolutionImporter::getStartTimeLevel()
{
	ifstream timeLevelFile( string(StringLiterals::outputFolder) + StringLiterals::timeLevelFile );
	if(!timeLevelFile)
		throw std::runtime_error("Could not open time level file.");
	long timeLevel;
	if(! (timeLevelFile >> timeLevel) )
		throw std::runtime_error("Could not read time level from file.");
	return timeLevel;
}

void SolutionImporter::importLiftDrag(vector<double> &lift, vector<double> &drag, long unsigned timeLevel)
{
	std::ifstream liftFile("./output/lift.dat");
	if(!liftFile)
		throw std::runtime_error("Didn't find lift file");
	double liftValue;
	while(liftFile >> liftValue)
	{
		lift.push_back(liftValue);
		if(lift.size() == timeLevel)
			break;
	}
	std::ifstream dragFile("./output/drag.dat");
	if(!dragFile)
		throw std::runtime_error("Didn't find drag file");
	double dragValue;
	while(dragFile >> dragValue)
	{
		drag.push_back(dragValue);
		if(drag.size() == timeLevel)
			break;
	}
}






