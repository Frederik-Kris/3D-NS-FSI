/*
 * Mesh.h
 *
 *  Created on: Oct 7, 2022
 *      Author: frederk
 */

#ifndef SRC_MESH_H_
#define SRC_MESH_H_

#include "Array3D.h"
#include "FlowVariableGroupStructs.h"
#include "includes_and_names.h"
#include "ConfigSettings.h"
#include "Boundary.h"

// The following typedefs are containers to store different derived classes (boundary types),
// to utilize polymorphism.
typedef vector<unique_ptr<MeshEdgeBoundary>> EdgeBoundaryCollection;
typedef vector<unique_ptr<ImmersedBoundary>> ImmersedBoundaryCollection;

// Class that keeps the actual arrays containing the flow variables, and the boundary conditions (BC)
class Mesh
{
public:

	Mesh(const ConfigSettings& params);

	void setupBoundaries(const ConfigSettings& params);

	void categorizeNodes(const ConfigSettings& params);

	bool checkFilterCondition(const ConfigSettings& params, ulong timeLevel, double t);

	void applyFilter_ifAppropriate(Array3D_d& variable_old, // <- To filter
								   Array3D_d& variable_tmp, // <- Temporary storage
								   const ConfigSettings& params,
								   ulong timeLevel,
								   double t);

	void swapConservedVariables();

	double getNormOfChange(const Array3D_d& oldValue, const Array3D_d& newValue);

	void applyAllBoundaryConditions(double t, const ConfigSettings& params);

	// Check whether any value is NaN or infinite.
	void checkFinity(const ConservedVariablesArrayGroup& arrays)
	{
		string message("Non-finite value found in ");
		if( !arrays.rho.allFinite() )
			cout << message << "rho.\n";
		if( !arrays.rho_u.allFinite() )
			cout << message << "rho_u.\n";
		if( !arrays.rho_v.allFinite() )
			cout << message << "rho_v.\n";
		if( !arrays.rho_w.allFinite() )
			cout << message << "rho_w.\n";
		if( !arrays.rho_E.allFinite() )
			cout << message << "rho_E.\n";
	}

	// Check whether any value is NaN or infinite.
	void checkFinity(const PrimitiveVariablesArrayGroup& arrays)
	{
		string message("Non-finite value found in ");
		if( !arrays.u.allFinite() )
			cout << message << "u.\n";
		if( !arrays.v.allFinite() )
			cout << message << "v.\n";
		if( !arrays.w.allFinite() )
			cout << message << "w.\n";
		if( !arrays.p.allFinite() )
			cout << message << "p.\n";
		if( !arrays.T.allFinite() )
			cout << message << "T.\n";
	}

	const ImmersedBoundaryCollection& getImmersedBoundaries() {return immersedBoundaries;}

	const size_t NI, NJ, NK;	// Mesh size. Number of nodes in x,y,z directions
	const size_t nNodesTotal;	// Total number of nodes in the mesh
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

#endif /* SRC_MESH_H_ */
