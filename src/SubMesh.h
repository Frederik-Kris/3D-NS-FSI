/*
 * SubMesh.h
 *
 *  Created on: Mar 6, 2023
 *      Author: frederk
 */

#ifndef SRC_SUBMESH_H_
#define SRC_SUBMESH_H_

#include "includes_and_names.h"
#include "FlowVariableGroupStructs.h"
#include "Boundary.h"

// The following typedefs are containers to store different derived classes (boundary types),
// to utilize polymorphism.
typedef vector<unique_ptr<MeshEdgeBoundary>> EdgeBoundaryCollection;
typedef vector<unique_ptr<ImmersedBoundary>> ImmersedBoundaryCollection;

// Class that keeps the actual arrays containing the flow variables, and the boundary conditions (BC)
class SubMesh
{

public:

	SubMesh(Vector3_i subMeshSize);

	SubMesh() : SubMesh( Vector3_i(0,0,0) ) {}

	void setSize(int NI, int NJ, int NK);

	void setBoundaries(EdgeBoundaryCollection& edgeBoundaries, ImmersedBoundaryCollection& immersedBoundaries);

	void categorizeNodes(const ConfigSettings& params);

	void swapConservedVariables();

	void applyAllBoundaryConditions(double t, const ConfigSettings& params);

	const ImmersedBoundaryCollection& getImmersedBoundaries() {return immersedBoundaries;}

	const int NI, NJ, NK;	// Mesh size. Number of nodes in x,y,z directions
	const int nNodesTotal;	// Total number of nodes in the mesh
	double dx, dy, dz;	// Grid spacing in x-, y- and z-direction
	ConservedVariablesArrayGroup conservedVariables;	// Mass density, momentum & total energy
	ConservedVariablesArrayGroup conservedVariablesOld;	// Previous time level. Or temporary storage, when needed.
	PrimitiveVariablesArrayGroup primitiveVariables;	// Velocity, pressure & Temperature
	TransportPropertiesArrayGroup transportProperties;	// Viscosity and thermal conductivity
	AllFlowVariablesArrayGroup flowVariableReferences;	// References to all flow variables
	RK4slopesArrayGroup RK4slopes;						// 4 slopes for each conserved variable
	IndexVectorGroup indexByType;	// 1D Indices to nodes of certain types.
	Array3D_nodeType nodeType;		// Type/category of each node (ghost, fluid, etc.)
	Vector3_d positionOffset;		// Offset of origin point, due to mesh edge boundary conditions.

private:

	EdgeBoundaryCollection edgeBoundaries;			// Boundaries at the edges of the Cartesian mesh
	ImmersedBoundaryCollection immersedBoundaries;	// Boundaries at immersed bodies

};


#endif /* SRC_SUBMESH_H_ */
