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
	FluidInterior,	// <- Active domain. Numerical stencil applied here.
	FluidEdge, 		// <- Claimed by a boundary, but flow variables are valid. (in-/outlet)
	FluidGhost, 	// <- Flow variables may not be valid. (periodic or symmetry BC, etc.)
	SolidInactive,	// <- Should never be accessed. Not part of stencil or BCs
	SolidGhost		// <- Related to an immersed boundary. Flow variables may be invalid.
};

// Package with vectors of indices to certain node types.
// Intent: Looping through all nodes of a given type without checking flags.
struct IndexVectorGroup
{
	vector<int> fluidInterior;
	vector<int> solidGhost;
};


#endif /* SRC_NODE_H_ */
