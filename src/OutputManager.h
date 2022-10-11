/*
 * OutputManager.h
 *
 *  Created on: Oct 7, 2022
 *      Author: frederk
 */

#ifndef SRC_OUTPUTMANAGER_H_
#define SRC_OUTPUTMANAGER_H_

#include "includes_and_names.h"
#include "Array3D_d.h"
#include "ConfigSettings.h"
#include "Mesh.h"

class OutputManager
{
public:
	OutputManager(const ConfigSettings& params);
	void processInitialOutput(const Mesh& mesh, double t);
	void processIntermediateOutput(const Mesh& mesh, Clock& statusReportTimer, double t, uint timeLevel, double dt);
	void processFinalOutput(const Mesh& mesh, double t, uint timeLevel, double dt,
			const vector<ConservedVariablesScalars> convergenceHistory);
private:
	void storeCurrentSolution_csv(const Mesh& mesh, double t);
	void storeCurrentSolution_csv_paraview(const Mesh& mesh);
	void storeCurrentSolution_csv_matlab(const Mesh& mesh);
	vector<const Array3D_d*> getPlotVariables(const Mesh& mesh);
	string get_csvHeaderString();
	vector<string> getVariableFileNames();
	void writePlaneTo_csv(ofstream& outputFile, const Array3D_d* flowVariable);
	void writeStatusReport_toScreen(double t, uint timeLevel, double dt);
	void writeOutputTimes();
	void writeConvergenceHistoryFiles(const vector<ConservedVariablesScalars>& convergenceHistory);

	const ConfigSettings params;	// Parameters and settings, imported from ConfigFile
	uint savedSolutions;			// No. of times saved to disk
	vector<double> outputTimes;     // The exact times when solution was saved
	Clock wallClockTimer;			// Timer that starts when simulation starts
};

#endif /* SRC_OUTPUTMANAGER_H_ */
