/*
 * SolutionImporter.h
 *
 *  Created on: Feb 8, 2023
 *      Author: frederk
 */

#ifndef SRC_SOLUTIONIMPORTER_H_
#define SRC_SOLUTIONIMPORTER_H_

#include "includes_and_names.h"
#include <vtk-7.1/vtkSmartPointer.h>
#include <vtk-7.1/vtkStructuredPointsReader.h>
#include <vtk-7.1/vtkPointData.h>
#include "FlowVariableGroupStructs.h"

// A class to handle the finding and reading of a previous solution from a .vtk file:
class SolutionImporter
{
public:

	SolutionImporter(const ConfigSettings& params);

	void loadVtkFile();

	bool allFlowVarsDerivable();

	void computeAllFlowVars(size_t i,
							ConservedVariablesScalars& conservedVars,
							PrimitiveVariablesScalars& primitiveVars,
							TransportPropertiesScalars& transportProps);

	void importNormHistories(ConservedVariablesVectorGroup& normHistory, ulong timeLevel);

	vector<double> getSolutionTimes();

	ulong getStartTimeLevel();

private:

	const ConfigSettings& params;

	vtkSmartPointer<vtkStructuredPointsReader> reader;
	vtkPointData* pointData;

	// Which flow variables are found in the vtk file:
	bool densityExists{false}, momentumExists{false}, energyExists{false};
	bool velocityExists{false}, pressureExists{false}, temperatureExists{false};

	bool vtkFileLoaded{false}; // False until we load the file

	string findLatestSolutionFile();

	void checkAvailableFlowVars();
};

#endif /* SRC_SOLUTIONIMPORTER_H_ */
