/*
 * Mesh.h
 *
 *  Created on: Oct 7, 2022
 *      Author: frederk
 */

#ifndef SRC_MESH_H_
#define SRC_MESH_H_

#include "includes_and_names.h"
#include "Array3D_d.h"
#include "FlowVariableScalarStructs.h"

class Mesh
{
public:
	Mesh(uint nMeshNodesX, uint nMeshNodesY, uint nMeshNodesZ,
		double domainLengthX, double domainLengthY, double domainLengthZ);
	void applyFilter_ifAppropriate(Array3D_d& variable_old, Array3D_d& variable_new,
									uint filterInterval, uint timeLevel);
	ConservedVariablesScalars computeNorms_conservedVariables();
	void swapConservedVariables();
	uint NI, NJ, NK;									// Mesh size. Number of nodes in x,y,z directions
	uint nNodesTotal;									// Total number of nodes in the mesh
	Array3D_d rho, rho_u, rho_v, rho_w, E;				// Conserved variables
	Array3D_d u, v, w, p, T;							// Primitive variables
	Array3D_d mu, kappa;								// Transport properties
	Array3D_d k1_rho  , k2_rho  , k3_rho  , k4_rho  ;	// RK4 stages, continuity
	Array3D_d k1_rho_u, k2_rho_u, k3_rho_u, k4_rho_u;	// RK4 stages, x-momentum
	Array3D_d k1_rho_v, k2_rho_v, k3_rho_v, k4_rho_v;	// RK4 stages, y-momentum
	Array3D_d k1_rho_w, k2_rho_w, k3_rho_w, k4_rho_w;	// RK4 stages, z-momentum
	Array3D_d k1_E    , k2_E    , k3_E    , k4_E    ;	// RK4 stages, specific total energy
	Array3D_d interm_rho, interm_rho_u, interm_rho_v,	// Intermediate solutions.
	interm_rho_w, interm_E;
	double dx, dy, dz;									// Grid spacing in x-, y- and z-direction
private:
	void setGridSpacings(double domainLengthX, double domainLengthY, double domainLengthZ );
	double getNormOfChange(const Array3D_d& oldValue, const Array3D_d& newValue);
};

#endif /* SRC_MESH_H_ */
