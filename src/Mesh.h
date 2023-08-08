/*
 * Mesh.h
 *
 *  Created on: Oct 7, 2022
 *      Author: frederk
 */

#ifndef SRC_MESH_H_
#define SRC_MESH_H_

#include "Array3D.h"
#include "BoundaryProxy.h"
#include "FlowVariableGroupStructs.h"
#include "includes_and_names.h"
#include "ConfigSettings.h"
#include "SubMesh.h"


// Class that represents the computational mesh. Also handles boundary conditions.
class Mesh
{
public:

	Mesh(const ConfigSettings& params);

	void setupBoundaries(const ConfigSettings& params);

	Array3D<SubMeshDescriptor> getSubMeshDescriptors(const IndexBoundingBox& subMeshIndices);

	void categorizeNodes(const ConfigSettings& params);

	bool checkFilterCondition(const ConfigSettings& params, long timeLevel, double t);

	void applyFilter_ifAppropriate(const ConfigSettings& params,
								   long timeLevel,
								   double t);

	void swapConservedVariables();

	void applyAllBoundaryConditions(double t, const ConfigSettings& params);

	void checkFinityAll() const
	{
		for(const SubMesh& subMesh : subMeshes)
		{
			checkFinity(subMesh.conservedVariables);
			checkFinity(subMesh.conservedVariablesOld);
			checkFinity(subMesh.primitiveVariables);
			checkFinity(subMesh.RK4slopes.k1);
			checkFinity(subMesh.RK4slopes.k2);
			checkFinity(subMesh.RK4slopes.k3);
			checkFinity(subMesh.RK4slopes.k4);
		}
	}

	// Check whether any value is NaN or infinite.
	void checkFinity(const ConservedVariablesArrayGroup& arrays) const
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
	void checkFinity(const PrimitiveVariablesArrayGroup& arrays) const
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

	// Find a sub-mesh which contains an immersed boundary:
	const SubMesh& getSubMeshWithIB() const
	{
		for(const SubMesh& subMesh : subMeshes)
			if( !subMesh.getImmersedBoundaries().empty() )
				return subMesh;
		throw std::runtime_error("Did not encounter any submesh with IBs.");
	}

	const int NI, NJ, NK;			// Base mesh size. Number of nodes in x,y,z directions, if refinement level is zero.
	Vector3_d smallestGridSpacings;	// Smallest dx, dy, dz, in all the sub-mesh regions.
	Vector3_i nRegions;				// Number of sub-meshes in each direction.
	Array3D<SubMesh> subMeshes;		// Regions in the mesh, containing the actual node data. Sub-meshes can have different refinement levels.
	Array3D<int> refinementLevels;	// Refinement levels of the submesh regions in the mesh
private:
	// Boundary conditions. These are just proxies to represent what BCs we want. The actual BCs are applied in the submeshes.
	EdgeBoundaryProxyCollection edgeBoundaries;			// Boundaries at the edges of the Cartesian mesh
	ImmersedBoundaryProxyCollection immersedBoundaries;	// Boundaries at immersed bodies

    vector<double> xThresholds;	// ↰
    vector<double> yThresholds;	// ← Coordinates for thresholds dividing submeshes
    vector<double> zThresholds;	// ↲
    vector<int> iThresholds;	// ↰
    vector<int> jThresholds;	// ← indices for thresholds dividing submeshes
    vector<int> kThresholds;	// ↲
    double dxBase;	// ↰
    double dyBase;	// ← Grid spacings at refinement level 0
    double dzBase;	// ↲

    void initializeThresholds(const ConfigSettings& params);
    void setRegionIDs();
	void makeThresholdsCoincideWithGridLines();
    void setRefinementLevels(const ConfigSettings& params);
    EdgeBoundaryCollection setupSubMeshEdgeBoundaries(int i, int j, int k, IndexBoundingBox& subMeshArrayLimits);
    void setupSubMeshes(const ConfigSettings& params);
    void findSmallestGridSpacings();
	void getRequestedThresholdsAndBaseGridSpacings(const ConfigSettings &params);
	void createBoundariesFromProxies(int i, int j, int k, IndexBoundingBox& subMeshArrayLimits, EdgeBoundaryCollection& subMeshEdgeBCs);
	void createSubmeshInterfacBoundaries(int i, int j, int k, EdgeBoundaryCollection& subMeshEdgeBCs, IndexBoundingBox& subMeshArrayLimits);
};

#endif /* SRC_MESH_H_ */


