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

// Subclass of Simulation, meant to take care of program output, e.g. writing to screen
// and storing results to disk.
class OutputManager
{
public:

	OutputManager(const ConfigSettings& params);

	void initialize(long unsigned timeLevelStart);

	void processInitialOutput(const Mesh& mesh, double t, long unsigned timeLevel);

	void processIntermediateOutput(const Mesh& mesh,
								   double t,
								   long unsigned timeLevel,
								   double dt);

	void processFinalOutput(const Mesh& mesh,
							double t,
							long unsigned timeLevel,
							double dt,
							const ConservedVariablesVectorGroup& convergenceHistory,
							const vector<double>& lift,
							const vector<double>& drag,
							const vector<double>& separationAngles);

private:

	void storeCurrentSolution(const Mesh& mesh, double t, long timeLevel);

	string getVtkHeader(const Mesh&, double t);

	void writeVtkNodeFlags(const Mesh&, ofstream& outputFile);

	vector<const Array3D_d*> getScalarVariablePointers(const Mesh& mesh);

	vector<string> getScalarVariableNames();

	void writeVtkScalarData(const Mesh&, ofstream& outputFile);

	vector<std::array<const Array3D_d*, 3>> getVectorVariablePointers(const Mesh& mesh);

	vector<string> getVectorVariableNames();

	void writeVtkVectorData(const Mesh&, ofstream& outputFile);

	void storeCurrentSolution_vtk(const Mesh& mesh, double t);

	void writeStatusReport_toScreen(double t, long timeLevel, double dt);

	string getTimeString(double seconds);

	void writeOutputTimes();

	void writeTimeLevel(int timeLevel);

	void writeIntegralProperties(const vector<double>& liftHistory,
								 const vector<double>& dragHistory,
								 const vector<double>& separationAngles);

	void writeConvergenceHistoryFiles(const ConservedVariablesVectorGroup& convergenceHistory);

	const ConfigSettings params;	// Parameters and settings, imported from config file
	int savedSolutions;			// No. of times saved to disk
	vector<double> outputTimes;     // The exact times when solution was saved
	long timeLevelStart;			// Time level when the simulation started. Zero unless simulation continued from earlier result.
	Clock wallClockTimer;			// Wall clock time for the entire simulation
	Clock statusReportTimer;		// Wall clock time since last status report to screen
	double previousProgression;		// Store progression between status reports, to estimate runtime
};

#endif /* SRC_OUTPUTMANAGER_H_ */



