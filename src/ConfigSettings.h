/*
 * ConfigSettings.h
 *
 *  Created on: Mar 15, 2021
 *      Author: frederik
 */

#ifndef SRC_CONFIGSETTINGS_H_
#define SRC_CONFIGSETTINGS_H_



#include "includes_and_names.h"
#include "SmallVectors.h"

enum class StopCriterionEnum
{
	timesteps, end_time
};

enum class saveNormalAxisEnum
{
	x, y, z
};

// What variables to store change-norm for:
enum class ConvHistoryEnum
{
	none, density, all
};

struct RefinementSpecification
{
	Vector3_i region;
	int level;
};

// This is a class to store ALL the settings imported from the config file.
// This is more convenient than constantly calling lookup("string") from a 'libconfig::Config' variable.
class ConfigSettings
{
public:
	ConfigSettings(string filename);
	void loadSettings(string filename);
	// MORE ELABORATE DESCRIPTIONS OF THESE PARAMETERS ARE IN THE CONFIG FILE ITSELF.
	bool continueSimulation;		 // Use newest output as initial condition
	int NI, NJ, NK;			     	 // Number of grid points in x-, y- and z-direction, respectively
	double L_x, L_y, L_z;			 // Dimensionless size of domain
	vector<double> refineBoundX;	 // Where to split domain into local grid refinement regions
	vector<double> refineBoundY;
	vector<double> refineBoundZ;
	vector<RefinementSpecification> specifiedRefinementLevels; // Regions where we specified refine level in config file


	StopCriterionEnum stopCriterion; // How the stopping criterion is defined
	long stopTimeLevel;  			 // Time level to stop simulation, if stopCriterion is 'timesteps'. Zero is IC.
	double t_end;                    // Time to stop simulation, if stopCriterion is 'end_time'.
	double convStabilityLimit;       // Constant specifying the inviscid stability criterion
	double viscStabilityLimit;       // Constant specifying the viscous stability criterion
	vector<int> filterIntervals;				// Number of timesteps between each time second order filter is applied
	vector<double> filterIntervalChangeTimes;	// Times to change the filter interval

	double statusReportInterval;	// How much wall clock time between each status/progression report to screen

	bool saveDensity, saveMomentum, saveEnergy; 		// Specifies what conserved variables to save to disk
	bool saveVelocity, savePressure, saveTemperature;	// Specifies what primitive variables to save to disk
	bool saveViscosity, saveThermalCond;				// Specifies what transport properties to save to disk

	bool saveIC, saveFinal; 		// Specifies whether to save IC and final solution
	bool saveIntervals; 			// Specifies whether to save solutions at given interval
	double savePeriod;				// How much time between each save, if save_intervals=true
	double saveIntervalsStartTime;	// When to start periodic saving, if save_intervals=true
	double saveIntervalsEndTime;	// When to stop periodic saving, if save_intervals=true
	ConvHistoryEnum saveConvergenceHistory; 	// Save 2-norm of change for each time step?

	double Gamma;				// Ratio of specific heats
	double Pr;					// Prandtl number, assumed constant
	double R;					// Specific gas constant
	double Re;					// Reynolds number. Scaled by ref.velocity.
	double Re_0;				// Stagnation Reynolds number. Scaled by ref. speed of sound.
	double sutherlands_C2;		// Second constant in Sutherlands law
	double M_0;					// Reference Mach number
	double T_0;					// Reference temperature, NOT DIMENSIONLESS

	bool errorOccurred;			// Error flag, set true if reading the 'ConfigFile' fails

	// Derived (not given in config file):

	double machinePrecisionBuffer;
private:
	void tryOpenConfigFile(Config& cfg, string filename);
	void tryReadSettingValues(Config& cfg);
	void readSettingValues(Config& cfg);
	void setDerivedParameters();
};



#endif /* SRC_CONFIGSETTINGS_H_ */
