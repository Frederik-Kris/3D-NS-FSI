/*
 * SubMesh.cpp
 *
 *  Created on: Mar 6, 2023
 *      Author: frederk
 */

#include "SubMesh.h"

// Default constructor. Setting all sizes to zero. Initializes the reference members.
SubMesh::SubMesh() :
  nNodes(0,0,0),
  arrayLimits(0,0,0),
  nNodesTotal{0},
  refinementLevel{0},
  gridSpacings(0,0,0),
  conservedVariables(arrayLimits),
  conservedVariablesOld(arrayLimits),
  primitiveVariables(arrayLimits),
  transportProperties(arrayLimits),
  flowVariableReferences(conservedVariables, primitiveVariables, transportProperties),
  RK4slopes(arrayLimits),
  nodeType(arrayLimits)
{}

// Set the number of mesh nodes and allocate space for this submesh.
void SubMesh::setSize(int nNodesX, int nNodesY, int nNodesZ,
					  const IndexBoundingBox& newArrayLimits,
					  const SpaceBoundingBox& spaceDomain)
{
	nNodes.i = nNodesX;
	nNodes.j = nNodesY;
	nNodes.k = nNodesZ;
	arrayLimits = newArrayLimits;
	nNodesTotal = arrayLimits.nNodesTotal();
	boundingBox = spaceDomain;
	gridSpacings.x = nNodes.i>1 ? (boundingBox.xMax - boundingBox.xMin)/(nNodes.i - 1) : 1;
	gridSpacings.y = nNodes.j>1 ? (boundingBox.yMax - boundingBox.yMin)/(nNodes.j - 1) : 1;
	gridSpacings.z = nNodes.k>1 ? (boundingBox.zMax - boundingBox.zMin)/(nNodes.k - 1) : 1;
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
		edgeBoundaries.push_back( boundary->getUniquePtrToCopy() );
	immersedBoundaries.clear();
	for(auto&& boundary : _immersedBoundaries)
		immersedBoundaries.push_back( boundary->getUniquePtrToCopy() );
}

// Categorize mesh nodes based on boundary conditions. This includes finding image point positions.
void SubMesh::categorizeNodes(const ConfigSettings& params,
							  const Array3D<SubMeshDescriptor>& neighborSubMeshes)
{
	nodeType.setAll(NodeTypeEnum::FluidInterior);

	IndexBoundingBox unclaimedNodes = arrayLimits;							// At first, all nodes are unclaimed.
	for(auto&& boundary : edgeBoundaries)									// Loop through mesh edge boundaries.
	{
		boundary->identifyRelatedNodes(unclaimedNodes, neighborSubMeshes);	// Each boundary flags and claims nodes.
	}
	std::reverse( edgeBoundaries.begin(), edgeBoundaries.end() );			// Reverse order, so last added is first applied.

	// Then each immersed boundary finds its solid and ghost nodes, and body intercept and image points:
	for(auto&& boundary : immersedBoundaries)
		boundary->identifyRelatedNodes(params, gridSpacings, nNodes, boundingBox.getMinPoint(), nodeType);

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


// Get info and flow variables for this sub-mesh region:
SubMeshDescriptor SubMesh::getSubMeshDescriptor()
{
	return SubMeshDescriptor( nNodes, gridSpacings, arrayLimits, boundingBox, nodeType, flowVariableReferences, refinementLevel, regionID);
}

// Loop through boundaries, who are derived classes of MeshEdgeBoundary and ImmersedBoundary,
// and apply their boundary conditions.
void SubMesh::applyAllBoundaryConditions(double t, const ConfigSettings& params)
{

	for(auto&& boundary : edgeBoundaries)
		boundary->applyBoundaryCondition(t, params);

	for(auto&& boundary : immersedBoundaries)
		boundary->applyBoundaryCondition(params);
}

// Compute the 2-norm of the difference between two arrys. Intended to monitor the change between two consecutive time levels,
// .g. to check convergence of solution. Could also be used to check error.
double SubMesh::getNormOfChange(const Array3D_d& oldValue, const Array3D_d& newValue) const
{
	double sumOfSquaredDifferences = 0;
	for(int i : indexByType.fluidInterior)
		sumOfSquaredDifferences += pow(oldValue(i)-newValue(i), 2);
	double dx = gridSpacings.x;
	double dy = gridSpacings.y;
	double dz = gridSpacings.z;
	double normOfChange = sqrt( dx*dy*dz * sumOfSquaredDifferences );
	return normOfChange;
}







