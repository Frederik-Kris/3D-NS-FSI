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

	void categorizeNodes(const ConfigSettings& params);

	void swapConservedVariables();

	void applyAllBoundaryConditions(double t, const ConfigSettings& params);

	const ImmersedBoundaryCollection& getImmersedBoundaries() {return immersedBoundaries;}

	ØØØ; // NEI! VI BYTTER TILBAKE NI ETC. TIL ANTALL NODER INNI REGIONEN. INDEX LIMITS BLIR HELLER ARRAY LIMITS.
	// DA MÅ KANSKJE CUSTOM ARRAY LIMITS INN I ARRAY CLASSEN OGSÅ.
	int NI, NJ, NK;					// Array size including ghost nodes. Number of nodes in x,y,z directions
	int nNodesTotal;				// Total number of nodes in the submesh
	IndexBoundingBox indexLimits;	// Index bounding box for the submesh region
	SpaceBoundingBox boundingBox;	// Coordinates to the limits of the region. Boundary nodes may fall outside depending on BC.
	double dx, dy, dz;				// Grid spacing in x-, y- and z-direction
	ConservedVariablesArrayGroup conservedVariables;	// Mass density, momentum & total energy
	ConservedVariablesArrayGroup conservedVariablesOld;	// Previous time level. Or temporary storage, when needed.
	PrimitiveVariablesArrayGroup primitiveVariables;	// Velocity, pressure & Temperature
	TransportPropertiesArrayGroup transportProperties;	// Viscosity and thermal conductivity
	AllFlowVariablesArrayGroup flowVariableReferences;	// References to all flow variables
	RK4slopesArrayGroup RK4slopes;						// 4 slopes for each conserved variable
	IndexVectorGroup indexByType;	// 1D Indices to nodes of certain types.
	Array3D_nodeType nodeType;		// Type/category of each node (ghost, fluid, etc.)

private:

	EdgeBoundaryCollection edgeBoundaries;			// Boundaries at the edges of the Cartesian mesh
	ImmersedBoundaryCollection immersedBoundaries;	// Boundaries at immersed bodies

};


#endif /* SRC_SUBMESH_H_ */
