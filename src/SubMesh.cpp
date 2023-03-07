/*
 * SubMesh.cpp
 *
 *  Created on: Mar 6, 2023
 *      Author: frederk
 */

#include "SubMesh.h"

// Constructor. Taking in no. of nodes in each direction.
SubMesh::SubMesh(Vector3_i subMeshSize)
: NI{subMeshSize.i}, NJ{subMeshSize.j}, NK{subMeshSize.k},
  nNodesTotal{NI*NJ*NK},
  dx{-1}, dy{-1}, dz{-1},
  conservedVariables(NI, NJ, NK),
  conservedVariablesOld(NI, NJ, NK),
  primitiveVariables(NI, NJ, NK),
  transportProperties(NI, NJ, NK),
  flowVariableReferences(conservedVariables, primitiveVariables, transportProperties),
  RK4slopes(NI, NJ, NK),
  nodeType(NI, NJ, NK)
{

}

// Categorize mesh nodes based on boundary conditions. This includes finding image point positions.
void SubMesh::categorizeNodes(const ConfigSettings& params)
{
	const Vector3_i nMeshNodes(NI, NJ, NK);
	const Vector3_d gridSpacing(dx, dy, dz);
	nodeType.setAll(NodeTypeEnum::FluidInterior);

	IndexBoundingBox unclaimedNodes(NI-1, NJ-1, NK-1);						// At first, all nodes are unclaimed.
	for(auto&& boundary : edgeBoundaries)									// Loop through mesh edge boundaries.
		boundary->identifyOwnedNodes(unclaimedNodes, nMeshNodes, nodeType);	// Each boundary flags and claims nodes.
	std::reverse( edgeBoundaries.begin(), edgeBoundaries.end() );			// Reverse order, so last added is first applied.

	// Then each immersed boundary finds its solid and ghost nodes, and body intercept and image points:
	for(auto&& boundary : immersedBoundaries)
		boundary->identifyRelatedNodes(params, gridSpacing, nMeshNodes, positionOffset, nodeType);

	// Finally we note the indices of interior fluid nodes and solid ghost nodes, so we can loop through only those:
	for(int index1D{0}; index1D<nNodesTotal; ++index1D)
		if(nodeType(index1D) == NodeTypeEnum::FluidInterior)
			indexByType.fluidInterior.push_back(index1D);
		else if(nodeType(index1D) == NodeTypeEnum::SolidGhost)
			indexByType.solidGhost.push_back(index1D);
}

// Swap the contents of all the arrays of conserved variables and the intermediate arrays, by move-semantics.
// This operation is super fast and needs no extra copy. Only the ownership of the data is changed.
void SubMesh::swapConservedVariables()
{
	conservedVariables.rho  .dataSwap(conservedVariablesOld.rho  );
	conservedVariables.rho_u.dataSwap(conservedVariablesOld.rho_u);
	conservedVariables.rho_v.dataSwap(conservedVariablesOld.rho_v);
	conservedVariables.rho_w.dataSwap(conservedVariablesOld.rho_w);
	conservedVariables.rho_E.dataSwap(conservedVariablesOld.rho_E);
}

// Loop through boundaries, who are derived classes of MeshEdgeBoundary and ImmersedBoundary,
// and apply their boundary conditions.
void SubMesh::applyAllBoundaryConditions(double t, const ConfigSettings& params)
{
	const Vector3_i nMeshNodes(NI, NJ, NK);
	const Vector3_d gridSpacing(dx, dy, dz);
	MeshDescriptor meshData(nMeshNodes, gridSpacing, positionOffset, nodeType, flowVariableReferences);

	for(auto&& boundary : edgeBoundaries)
		boundary->applyBoundaryCondition(t, nMeshNodes, params, flowVariableReferences);

	for(auto&& boundary : immersedBoundaries)
		boundary->applyBoundaryCondition(params, meshData);
}
