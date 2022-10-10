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

class OutputManager
{
public:
	OutputManager();
	void processInitialOutput(const ConfigSettings& params);
	void processIntermediateOutput(const ConfigSettings& params, Clock& statusReportTimer, double t, double dt);
	void processFinalOutput(const ConfigSettings& params);
	void checkMassConservation(double& inFluxSum, double& outFluxSum);
private:
	void storeCurrentSolution_csv();
	void storeCurrentSolution_csv_paraview();
	void storeCurrentSolution_csv_matlab();
	vector<Array3D_d*> getPlotVariables();
	string get_csvHeaderString();
	vector<string> getVariableFileNames();
	void writePlaneTo_csv(ofstream& outputFile, Array3D_d* flowVariable);
	void writeStatusReport_toScreen();
	void writeOutputTimes();
	void writeNormHistoryFiles();

	uint savedSolutions;                                // No. of times saved to disk
	vector<double> outputTimes;                         // The exact times when solution was saved
};

#endif /* SRC_OUTPUTMANAGER_H_ */
