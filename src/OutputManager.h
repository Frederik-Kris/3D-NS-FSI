/*
 * OutputManager.h
 *
 *  Created on: Oct 7, 2022
 *      Author: frederk
 */

#ifndef SRC_OUTPUTMANAGER_H_
#define SRC_OUTPUTMANAGER_H_

#include <filesystem>
#include "Array3D.h"
#include "includes_and_names.h"
#include "ConfigSettings.h"
#include "Mesh.h"

class OutputManager
{
public:
	OutputManager(const ConfigSettings& params);

	void initialize();

	void processInitialOutput(const Mesh& mesh, double t);

	void processIntermediateOutput(const Mesh& mesh, Clock& statusReportTimer,
									double t, ulong timeLevel, double dt);

	void processFinalOutput(const Mesh& mesh, double t, ulong timeLevel, double dt,
							const ConservedVariablesVectorGroup& convergenceHistory);
private:
	void storeCurrentSolution(const Mesh& mesh, double t);

	string getVtkHeader(const Mesh&, double t);

	void writeVtkNodeFlags(const Mesh&, ofstream& outputFile);

	vector<const Array3D_d*> getScalarVariablePointers(const Mesh& mesh);

	vector<string> getScalarVariableNames();

	void writeVtkScalarData(const Mesh&, ofstream& outputFile);

	vector<std::array<const Array3D_d*, 3>> getVectorVariablePointers(const Mesh& mesh);

	vector<string> getVectorVariableNames();

	void writeVtkVectorData(const Mesh&, ofstream& outputFile);

	void storeCurrentSolution_vtk(const Mesh& mesh, double t);

	void writeStatusReport_toScreen(double t, ulong timeLevel, double dt);

	void writeOutputTimes();

	void writeConvergenceHistoryFiles(const ConservedVariablesVectorGroup& convergenceHistory);

	const ConfigSettings params;	// Parameters and settings, imported from ConfigFile
	uint savedSolutions;			// No. of times saved to disk
	vector<double> outputTimes;     // The exact times when solution was saved
	Clock wallClockTimer;			// Timer that starts when simulation starts
};

#endif /* SRC_OUTPUTMANAGER_H_ */
