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

	SubMesh();

	void setSize(int NI, int NJ, int NK,
			     const IndexBoundingBox& indexDomain,
				 const SpaceBoundingBox& spaceDomain);

	void setBoundaries(EdgeBoundaryCollection& edgeBoundaries, ImmersedBoundaryCollection& immersedBoundaries);

	void categorizeNodes(const ConfigSettings& params, Array3D<AllFlowVariablesArrayGroup>& neighborSubMeshReferences);

	void swapConservedVariables();

	void applyAllBoundaryConditions(double t, const ConfigSettings& params);

	const ImmersedBoundaryCollection& getImmersedBoundaries() {return immersedBoundaries;}

	Vector3_i regionID;				// Indices to this sub-mesh in a 3D array of sub-meshes.
	int NI, NJ, NK;					// Number of nodes in x,y,z directions in the submesh, not including overlap nodes OUTSIDE of the region.
	IndexBoundingBox arrayLimits;	// Index limits for the submesh array
	int nNodesTotal;				// Total number of nodes in the submesh array
	SpaceBoundingBox boundingBox;	// Coordinates to the limits of the region. Boundary nodes may fall outside depending on BC.
	double dx, dy, dz;				// Grid spacing in x-, y- and z-direction
	ConservedVariablesArrayGroup conservedVariables;	// Mass density, momentum & total energy
	ConservedVariablesArrayGroup conservedVariablesOld;	// Previous time level. Or temporary storage, when needed.
	PrimitiveVariablesArrayGroup primitiveVariables;	// Velocity, pressure & Temperature
	TransportPropertiesArrayGroup transportProperties;	// Viscosity and thermal conductivity
	AllFlowVariablesArrayGroup flowVariableReferences;	// References to all flow variables
	RK4slopesArrayGroup RK4slopes;						// 4 slopes for each conserved variable
	IndexVectorGroup indexByType;	// 1D Indices to nodes of certain types.
	Array3D<NodeTypeEnum> nodeType;	// Type/category of each node (ghost, fluid, etc.)

private:

	EdgeBoundaryCollection edgeBoundaries;			// Boundaries at the edges of the Cartesian mesh
	ImmersedBoundaryCollection immersedBoundaries;	// Boundaries at immersed bodies

};


#endif /* SRC_SUBMESH_H_ */
