/*
 * ConfigSettings.cpp
 *
 *  Created on: Mar 15, 2021
 *      Author: frederik
 */

#include "ConfigSettings.h"

// Default constructor. Simply loads settings from ConfigFile, if possible.
ConfigSettings::ConfigSettings(string filename)
{
	loadSettings(filename);
}

// This function tries to load the settings from the 'ConfigFile', using the 'libconfig' library.
// If something goes wrong, the 'error' member is assigned true. Caller should check this.
void ConfigSettings::loadSettings(string filename)
{
	errorOccurred = false;	// <- until the opposite is proven
	Config cfg;
	tryOpenConfigFile(cfg, filename);

	if (!errorOccurred)
		tryReadSettingValues(cfg);
	if (!errorOccurred)
		setDerivedParameters();
}

// Tries to open and parse the 'ConfigFile'. If the file can't be opened or parsed, 'error' is assigned true.
void ConfigSettings::tryOpenConfigFile(Config& cfg, string filename)
{
	try
	{
		cfg.readFile(filename.c_str());
	}
	catch (libconfig::FileIOException& ex)
	{
		cout << "File I/O error: " << ex.what() << endl
			 << "Make sure the 'ConfigFile' exists in correct path." << endl;
		errorOccurred = true;
	}
	catch (libconfig::ParseException& ex)
	{
		cout << "Parse error in ConfigFile: " << ex.what() << endl
			 << "See documentation for 'libconfig' library and ensure that the 'ConfigFile' is structured correctly." << endl;
		errorOccurred = true;
	}
}

// Tries to call 'readSettingValues'. If a setting can't be found or has the wrong type 'error' is assigned true.
void ConfigSettings::tryReadSettingValues(Config& cfg)
{
	try
	{
		readSettingValues(cfg);
	}
	catch (libconfig::SettingNotFoundException& ex)
	{
		cout << "ConfigFile setting not found: " << ex.what() << endl
			 << "Did you specify the entire setting path? (group.setting)" << endl
			 << "Is the name in 'ConfigFile' the same as in the 'lookup' in the code?" << endl;
		errorOccurred = true;
	}
	catch (libconfig::SettingTypeException& ex)
	{
		cout << "ConfigFile setting type mismatch: " << ex.what() << endl
			 << "The return value of a 'lookup' was probably assigned to a variable of another type." << endl
			 << "Integers must NOT have a decimal point. Floating points must have a decimal point (1. or 1.0, not 1)."
			 << "Strings must be in \"double quotes\"." << endl
			 << "See the documentation of the 'libconfig' library for info on how to specify data types." << endl;
		errorOccurred = true;
	}
}

// Looks up settings in the Config variable 'cfg' and assign them to members in the ConfigSettings
void ConfigSettings::readSettingValues(Config& cfg)
{
	NI = (uint) cfg.lookup("NI");
	NJ = (uint) cfg.lookup("NJ");
	NK = (uint) cfg.lookup("NK");

	string stopCriterionString = cfg.lookup("stopCriterion").c_str();
	if (stopCriterionString == "timesteps")
	{
		stopCriterion = StopCriterionEnum::timesteps;
		stopTimeLevel = (uint) cfg.lookup("stopTimeLevel");
	}
	else if (stopCriterionString == "end_time")
	{
		stopCriterion = StopCriterionEnum::end_time;
		t_end = cfg.lookup("t_end");
	}
	else
	{
		cout << "Error when reading ConfigFile. 'stopCriterion' needs to be 'timesteps' or 'end_time'." << endl;
		errorOccurred = true;
	}

	convStabilityLimit = cfg.lookup("convStabilityLimit");
	viscStabilityLimit = cfg.lookup("viscStabilityLimit");
	filterInterval	   = cfg.lookup("filterInterval");

	statusReportInterval = cfg.lookup("statusReportInterval");

	save_rho       = cfg.lookup("save_rho");
	save_rho_u     = cfg.lookup("save_rho_u");
	save_rho_v     = cfg.lookup("save_rho_v");
	save_rho_w     = cfg.lookup("save_rho_w");
	save_E         = cfg.lookup("save_E");
	save_u         = cfg.lookup("save_u");
	save_v         = cfg.lookup("save_v");
	save_w         = cfg.lookup("save_w");
	save_p         = cfg.lookup("save_p");
	save_T         = cfg.lookup("save_T");
	save_mu        = cfg.lookup("save_mu");
	save_kappa     = cfg.lookup("save_kappa");
	saveIC        = cfg.lookup("saveIC");
	saveFinal     = cfg.lookup("saveFinal");
	saveIntervals = cfg.lookup("saveIntervals");
	savePeriod    = cfg.lookup("savePeriod");
	saveIntervalsStartTime = cfg.lookup("saveIntervalsStartTime");
	saveIntervalsEndTime   = cfg.lookup("saveIntervalsEndTime");
	saveForParaview= cfg.lookup("saveForParaview");
	saveForMatlab  = cfg.lookup("saveForMatlab");

	if(saveForMatlab)
	{
		string saveNormalAxisString = cfg.lookup("saveNormalAxis").c_str();
		if(saveNormalAxisString=="x")
			saveNormalAxis = saveNormalAxisEnum::x;
		else if(saveNormalAxisString=="y")
			saveNormalAxis = saveNormalAxisEnum::y;
		else if(saveNormalAxisString=="z")
			saveNormalAxis = saveNormalAxisEnum::z;
		else
		{
			errorOccurred = true;
			std::cerr << "saveNormalAxis must be x, y, or z.\n";
		}
		saveConstantIndex = (uint) cfg.lookup("saveConstantIndex");
	}

	string saveConvergenceHistoryString = cfg.lookup("saveConvergenceHistory").c_str();
	if(saveConvergenceHistoryString == "none")
		saveConvergenceHistory = ConvHistoryEnum::none;
	else if (saveConvergenceHistoryString == "density")
		saveConvergenceHistory = ConvHistoryEnum::density;
	else if (saveConvergenceHistoryString == "all")
		saveConvergenceHistory = ConvHistoryEnum::all;
	else
	{
		errorOccurred = true;
		std::cerr << "saveConvergenceHistory must be 'none', 'density' or 'all'. \n";
	}

	Gamma 			= cfg.lookup("Gamma");
	Pr				= cfg.lookup("Pr");
	R				= cfg.lookup("R");
	Re				= cfg.lookup("Re");
	sutherlands_C2	= cfg.lookup("sutherlands_C2");
	M_0				= cfg.lookup("M_0");
	T_0				= cfg.lookup("T_0");
	L_x				= cfg.lookup("L_x");
	L_y				= cfg.lookup("L_y");
	L_z				= cfg.lookup("L_z");
}

// Set parameters that are not given in the config file.
// They are derived from the given values etc.
void ConfigSettings::setDerivedParameters()
{
	machinePrecisionBuffer = std::numeric_limits<double>::epsilon()*(L_x+L_y+L_z);
}


