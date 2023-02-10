/*
 * ConstantLiteralExpressions.h
 *
 *  Created on: Feb 7, 2023
 *      Author: frederk
 */

#ifndef SRC_CONSTANTLITERALEXPRESSIONS_H_
#define SRC_CONSTANTLITERALEXPRESSIONS_H_

/*  Strings that describe variables or concepts in the program. By declaring them here we never have to
 *  worry about spelling it the same way. Then we can use these names when saving/loading results,
 *  or even as keys in associative containers.
 */
namespace StringLiterals
{
constexpr char densityString 	 [] = "Density";
constexpr char energyString 	 [] = "Total_Specific_Energy_per_Volume";
constexpr char pressureString 	 [] = "Pressure";
constexpr char temperatureString [] = "Temperature";
constexpr char viscosityString   [] = "Dynamic_Viscosity";
constexpr char thermalCondString [] = "Thermal_Conductivity";
constexpr char velocityString 	 [] = "Velocity";
constexpr char momentumString 	 [] = "Momentum_Density";

constexpr char outputFolder 		[] = "./output/";
constexpr char solutionFileNameBase [] = "out.vtk.";
constexpr char timesFile 			[] = "times.dat";
constexpr char normRhoFile   		[] = "norm_rho.dat";
constexpr char normRhoUFile 		[] = "norm_rho_u.dat";
constexpr char normRhoVFile 		[] = "norm_rho_v.dat";
constexpr char normRhoWFile 		[] = "norm_rho_w.dat";
constexpr char normRhoEFile 		[] = "norm_rho_E.dat";
}



#endif /* SRC_CONSTANTLITERALEXPRESSIONS_H_ */
