/*
 * FlowVariableScalarStructs.h
 *
 *  Created on: Oct 10, 2022
 *      Author: frederk
 */

#ifndef SRC_FLOWVARIABLESCALARSTRUCTS_H_
#define SRC_FLOWVARIABLESCALARSTRUCTS_H_


struct ConservedVariablesScalars
{
public:
	ConservedVariablesScalars(double rho, double rho_u, double rho_v, double rho_w, double E) :
		rho{rho}, rho_u{rho_u}, rho_v{rho_v}, rho_w{rho_w}, E{E}
	{}
	double rho;
	double rho_u;
	double rho_v;
	double rho_w;
	double E;
};


#endif /* SRC_FLOWVARIABLESCALARSTRUCTS_H_ */
