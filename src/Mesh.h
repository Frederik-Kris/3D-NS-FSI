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
#include "SubMesh.h"

// Simple container of submeshes with indexing in three directions.
// Could have used an std::vector<SubMesh> instead, but 3D indexing is handy for neighbor-finding.
class SubMeshCollection
{
public:

	SubMeshCollection(int nRegionsX, int nRegionsY, int nRegionsZ)
	: nRegionsX{nRegionsX},
	  nRegionsY{nRegionsY},
	  nRegionsZ{nRegionsZ},
	  subMeshes(nRegionsX * nRegionsY * nRegionsZ)
	{}

	SubMeshCollection() : SubMeshCollection(0,0,0) {}

	void setNumberOfRegions(int _nRegionsX, int _nRegionsY, int _nRegionsZ)
	{
		nRegionsX = _nRegionsX;
		nRegionsY = _nRegionsY;
		nRegionsZ = _nRegionsZ;
		subMeshes = vector<SubMesh>(nRegionsX * nRegionsY * nRegionsZ);
	}

	// Access submesh (lvalue)
	SubMesh& operator()(int i, int j, int k)
	{ return subMeshes.at( i*nRegionsY*nRegionsZ + j*nRegionsZ + k ); }

	// Access submesh (rvalue)
	const SubMesh& operator()(int i, int j, int k) const
	{ return subMeshes.at( i*nRegionsY*nRegionsZ + j*nRegionsZ + k ); }

	const int nRegionsX, nRegionsY, nRegionsZ; // No. of regions

private:
	vector<SubMesh> subMeshes; // 1D array with the actual submeshes.
};

// Class that represents the computational mesh. Also handles boundary conditions.
class Mesh
{
public:

	Mesh(const ConfigSettings& params);

	void setupBoundaries(const ConfigSettings& params);

	Array3D<AllFlowVariablesArrayGroup> getNeighborSubMeshVariables(const IndexBoundingBox& neighborSubMeshIndices);

	void categorizeNodes();

	bool checkFilterCondition(const ConfigSettings& params, long timeLevel, double t);

	void applyFilter_ifAppropriate(Array3D_d& variable_old, // <- To filter
								   Array3D_d& variable_tmp, // <- Temporary storage
								   const ConfigSettings& params,
								   long timeLevel,
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

	const int NI, NJ, NK;			// Base mesh size. Number of nodes in x,y,z directions, if refinement level is zero.
	Vector3_i nRegions;				// Number of sub-meshes in each direction.
	Array3D<SubMesh> subMeshes;		// Regions in the mesh, containing the actual node data. Sub-meshes can have different refinement levels.
	Array3D<int> refinementLevels;	// Refinement levels of the submesh regions in the mesh
private:
	// Boundary conditions. These are just proxies to represent what BCs we want. The actual BCs are applied in the submeshes.
	// TODO: consider making a separate class for these proxies, for clarity, so we don't try to applyBC() from these by accident.
	EdgeBoundaryCollection edgeBoundaries;			// Boundaries at the edges of the Cartesian mesh
	ImmersedBoundaryCollection immersedBoundaries;	// Boundaries at immersed bodies
};

#endif /* SRC_MESH_H_ */


