/*
 * Node.h
 *
 *  Created on: Oct 19, 2022
 *      Author: frederk
 */

#ifndef SRC_NODE_H_
#define SRC_NODE_H_

#include "includes_and_names.h"
#include "SmallVectors.h"
#include "BoundingBox.h"

// Struct representing a solid ghost node related to an immersed surface
struct GhostNode
{
	Vector3_i indices;
	Vector3_d bodyInterceptPoint;
	Vector3_d imagePoint;

	GhostNode(Vector3_i indices) : indices(indices) {}
};

typedef std::vector<GhostNode>::iterator GhostNodeVectorIterator;

// Classification of mesh nodes
enum class NodeTypeEnum
{
	FluidInterior,	// ← Active domain. Numerical stencil applied here.
	FluidEdge, 		// ← Claimed by a boundary, but flow variables are valid. (in-/outlet)
	FluidGhost, 	// ← Flow variables may not be valid. (periodic or symmetry BC, etc.)
	SolidInactive,	// ← Should never be accessed. Not part of stencil or BCs
	SolidGhost		// ← Related to an immersed boundary. Flow variables may be invalid.
};

// Package with vectors of indices to certain node types.
// Intent: Looping through all nodes of a given type without checking flags.
struct IndexVectorGroup
{
	vector<int> fluidInterior;
	vector<int> solidGhost;
};

// Given the 1D index to a node, get the 3D indices
inline Vector3_i getIndices3D(int index1D, const IndexBoundingBox& arrayLimits)
{
	int nNodesJ = arrayLimits.jMax - arrayLimits.jMin + 1;
	int nNodesK = arrayLimits.kMax - arrayLimits.kMin + 1;
	int i = index1D / (nNodesJ * nNodesK) + arrayLimits.iMin;
	int j = index1D % (nNodesJ * nNodesK) / nNodesK + arrayLimits.jMin;
	int k = index1D % (nNodesJ * nNodesK) % nNodesK + arrayLimits.kMin;
	return Vector3_i(i,j,k);
}

// Given the 3D indices to a node, get the 1D index
inline int getIndex1D(int i, int j, int k, const IndexBoundingBox& arrayLimits)
{
	int nNodesJ = arrayLimits.jMax - arrayLimits.jMin + 1;
	int nNodesK = arrayLimits.kMax - arrayLimits.kMin + 1;
	return (i-arrayLimits.iMin)*nNodesJ*nNodesK + (j-arrayLimits.jMin)*nNodesK + (k-arrayLimits.kMin);
}

// Given the 3D indices to a node, get the 1D index
inline int getIndex1D(const Vector3_i& indices, const IndexBoundingBox& arrayLimits)
{ return getIndex1D(indices.i, indices.j, indices.k, arrayLimits); }

// Get the position (coordinates) of the node with the given indices
inline Vector3_d getNodePosition(int i, int j, int k,
								 const IndexBoundingBox& arrayLimits,
								 const Vector3_d& gridSpacings,
								 const Vector3_d& meshOriginOffset)
{
	double x = (i + arrayLimits.iMin) * gridSpacings.x + meshOriginOffset.x ;
	double y = (j + arrayLimits.jMin) * gridSpacings.y + meshOriginOffset.y ;
	double z = (k + arrayLimits.kMin) * gridSpacings.z + meshOriginOffset.z ;
	return Vector3_d(x, y, z);
}

// Get the position (coordinates) of the node with the given indices
inline Vector3_d getNodePosition(const Vector3_i& indices,
								 const IndexBoundingBox& arrayLimits,
								 const Vector3_d& gridSpacings,
								 const Vector3_d& meshOriginOffset)
{ return getNodePosition(indices.i, indices.j, indices.k, arrayLimits, gridSpacings, meshOriginOffset); }


#endif /* SRC_NODE_H_ */
