/*
 * Mesh.h
 *
 *  Created on: Oct 7, 2022
 *      Author: frederk
 */

#ifndef SRC_MESH_H_
#define SRC_MESH_H_

#include "includes_and_names.h"
#include "ConfigSettings.h"
#include "Array3D_d.h"
#include "Boundary.h"
#include "FlowVariableGroupStructs.h"

struct IndexVectorGroup
{
	vector<uint> fluidNodes;
	vector<uint> activeNodes;
	vector<uint> ghostNodes;
	vector<uint> surfaceNodes;
};

class Mesh
{
public:
	Mesh(const ConfigSettings& params);
	void setupBoundaries(const ConfigSettings& params);
	void assertBoundaryConditionCompliance();
	void categorizeNodes();
	void applyFilter_ifAppropriate(Array3D_d& variable_old, Array3D_d& variable_new,
									uint filterInterval, uint timeLevel);
	ConservedVariablesScalars computeNorms_conservedVariables();
	void swapConservedVariables();
	uint NI, NJ, NK;	// Mesh size. Number of nodes in x,y,z directions
	uint nNodesTotal;	// Total number of nodes in the mesh
	double dx, dy, dz;	// Grid spacing in x-, y- and z-direction
	ConservedVariablesArrayGroup conservedVariables;				// Mass density, momentum & total energy
	PrimitiveVariablesArrayGroup primitiveVariables;				// Velocity, pressure & Temperature
	TransportPropertiesArrayGroup transportProperties;				// Viscosity and thermal conductivity
	ConservedVariablesArrayGroup intermediateConservedVariables;	// Intermediate states of conserved variables.
	RK4slopesArrayGroup RK4slopes;									// 4 slopes for each conserved variable
	IndexVectorGroup nodeIndices;									// Indices to different types of nodes
private:
	vector<unique_ptr<MeshEdgeBoundary>> edgeBoundaries;		// Boundaries at the edges of the Cartesian mesh
	vector<unique_ptr<ImmersedBoundary>> immersedBoundaries;	// Boundaries at immersed bodies
	void setGridSpacings(double domainLengthX, double domainLengthY, double domainLengthZ );
	double getNormOfChange(const Array3D_d& oldValue, const Array3D_d& newValue);
};

#endif /* SRC_MESH_H_ */
