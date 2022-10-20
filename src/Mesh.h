/*
 * Mesh.h
 *
 *  Created on: Oct 7, 2022
 *      Author: frederk
 */

#ifndef SRC_MESH_H_
#define SRC_MESH_H_

class Mesh;

#include "Array3D.h"
#include "includes_and_names.h"
#include "ConfigSettings.h"
#include "Boundary.h"
#include "FlowVariableGroupStructs.h"

typedef vector<unique_ptr<MeshEdgeBoundary>> EdgeBoundaryCollection;
typedef vector<unique_ptr<ImmersedBoundary>> ImmersedBoundaryCollection;

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

	void categorizeNodes(const ConfigSettings& params);

	void applyFilter_ifAppropriate(Array3D_d& variable_old, Array3D_d& variable_new,
									uint filterInterval, uint timeLevel);

	ConservedVariablesScalars computeNorms_conservedVariables();

	void swapConservedVariables();

	Vector3_u getIndices3D(uint index1D) const;

	uint getIndex1D(uint i, uint j, uint k) const;

	uint getIndex1D(Vector3_u indices) const
	{ return getIndex1D(indices.i, indices.j, indices.k); }

	Vector3_d getNodePosition(uint i, uint j, uint k) const;

	Vector3_d getNodePosition(Vector3_u indices) const
	{ return getNodePosition(indices.i, indices.j, indices.k); }

	IndexBoundingBox getSurroundingNodesBox(Vector3_d point) const;

	uint NI, NJ, NK;	// Mesh size. Number of nodes in x,y,z directions
	uint nNodesTotal;	// Total number of nodes in the mesh
	double dx, dy, dz;	// Grid spacing in x-, y- and z-direction
	ConservedVariablesArrayGroup conservedVariables;				// Mass density, momentum & total energy
	PrimitiveVariablesArrayGroup primitiveVariables;				// Velocity, pressure & Temperature
	TransportPropertiesArrayGroup transportProperties;				// Viscosity and thermal conductivity
	ConservedVariablesArrayGroup intermediateConservedVariables;	// Intermediate states of conserved variables.
	RK4slopesArrayGroup RK4slopes;									// 4 slopes for each conserved variable
	IndexVectorGroup nodeIndices;									// Indices to different types of nodes
	Array3D_nodeType nodeTypes;										// Type/category of each node (ghost, fluid, etc.)
private:
	EdgeBoundaryCollection edgeBoundaries;			// Boundaries at the edges of the Cartesian mesh
	ImmersedBoundaryCollection immersedBoundaries;	// Boundaries at immersed bodies

	void setGridSpacings(double domainLengthX, double domainLengthY, double domainLengthZ );

	double getNormOfChange(const Array3D_d& oldValue, const Array3D_d& newValue);
};

#endif /* SRC_MESH_H_ */
