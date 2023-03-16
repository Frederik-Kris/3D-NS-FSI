/*
 * SubMesh.cpp
 *
 *  Created on: Mar 6, 2023
 *      Author: frederk
 */

#include "SubMesh.h"

// Default constructor. Setting all sizes to zero. Initializes the reference members.
SubMesh::SubMesh() :
  NI{0}, NJ{0}, NK{0},
  nNodesTotal{0},
  dx{0}, dy{0}, dz{0},
  conservedVariables(NI, NJ, NK),
  conservedVariablesOld(NI, NJ, NK),
  primitiveVariables(NI, NJ, NK),
  transportProperties(NI, NJ, NK),
  flowVariableReferences(conservedVariables, primitiveVariables, transportProperties),
  RK4slopes(NI, NJ, NK),
  nodeType(NI, NJ, NK)
{}

// Set the number of mesh nodes and allocate space for this submesh.
void SubMesh::setSize(int nNodesX, int nNodesY, int nNodesZ,
					  const IndexBoundingBox& newArrayLimits,
					  const SpaceBoundingBox& spaceDomain)
{
	NI = nNodesX;
	NJ = nNodesY;
	NK = nNodesZ;
	arrayLimits = newArrayLimits;
	nNodesTotal = arrayLimits.nNodesTotal();
	boundingBox = spaceDomain;
	dx = NI>1 ? (boundingBox.xMax - boundingBox.xMin)/(NI - 1) : 1;
	dy = NJ>1 ? (boundingBox.yMax - boundingBox.yMin)/(NJ - 1) : 1;
	dz = NK>1 ? (boundingBox.zMax - boundingBox.zMin)/(NK - 1) : 1;
	conservedVariables    = ConservedVariablesArrayGroup(arrayLimits);
	conservedVariablesOld = ConservedVariablesArrayGroup(arrayLimits);
	primitiveVariables    = PrimitiveVariablesArrayGroup(arrayLimits);
	transportProperties = TransportPropertiesArrayGroup(arrayLimits);
	RK4slopes = RK4slopesArrayGroup(arrayLimits);
	nodeType = Array3D<NodeTypeEnum>(arrayLimits);
}

// Set the BCs for this submesh
void SubMesh::setBoundaries(EdgeBoundaryCollection& _edgeBoundaries, ImmersedBoundaryCollection& _immersedBoundaries)
{
	edgeBoundaries.clear();
	for(auto&& boundary : _edgeBoundaries)
		edgeBoundaries.push_back( boundary );
	immersedBoundaries.clear();
	for(auto&& boundary : _immersedBoundaries)
		immersedBoundaries.push_back( boundary );
}

// Categorize mesh nodes based on boundary conditions. This includes finding image point positions.
void SubMesh::categorizeNodes(const ConfigSettings& params, Array3D<AllFlowVariablesArrayGroup>& neighborSubMeshReferences)
{
	const Vector3_d gridSpacing(dx, dy, dz);
	nodeType.setAll(NodeTypeEnum::FluidInterior);

	IndexBoundingBox unclaimedNodes = arrayLimits;						// At first, all nodes are unclaimed.
	for(auto&& boundary : edgeBoundaries)									// Loop through mesh edge boundaries.
	{
		boundary->identifyOwnedNodes(unclaimedNodes, arrayLimits, nodeType);	// Each boundary flags and claims nodes.
		boundary->identifyBorrowedNodes(neighborSubMeshReferences);			// Find nodes we need in neighbor submeshes
	}
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
