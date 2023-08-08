/*
 * BoundingBox.h
 *
 *  Created on: May 25, 2023
 *      Author: frederk
 */

#ifndef SRC_BOUNDINGBOX_H_
#define SRC_BOUNDINGBOX_H_


#include "includes_and_names.h"
#include "SmallVectors.h"

// Struct that represents a bounded region of the mesh, defined by lowest and highest indices
struct IndexBoundingBox
{
	int iMin, iMax;	// Indices in x-direction
	int jMin, jMax;	// Indices in y-direction
	int kMin, kMax;	// Indices in z-direction

	// Constructor. Takes in the upper index limits, and sets the lower limits to zero.
	IndexBoundingBox(int iMax, int jMax, int kMax)
	: iMin{0}, iMax{iMax},
	  jMin{0}, jMax{jMax},
	  kMin{0}, kMax{kMax}
	{}

	// Constructor. Takes in all index limits.
	IndexBoundingBox(int iMin, int iMax, int jMin, int jMax, int kMin, int kMax)
	: iMin{iMin}, iMax{iMax},
	  jMin{jMin}, jMax{jMax},
	  kMin{kMin}, kMax{kMax}
	{}

	// Default constructor.
	IndexBoundingBox() : IndexBoundingBox(0,0,0,0,0,0) {}

	// Get a simple fixed size array with the 1D indices to the 8 corner nodes in the box
	// Ordering: z(k) changes most often, x(i) changes most seldom.
	Array8_u cornersAsIndexList(const IndexBoundingBox& arrayLimits) const;

	// Get a vector to the 1D indices to all nodes.
	// The 1D indices are with respect to an array bounded by "arrayLimits".
	vector<int> allNodesAsIndexList(const IndexBoundingBox& arrayLimits) const;

	// Get a vector to the 1D indices to all nodes, except those who are also entailed by another given box.
	// The 1D indices are with respect to an array bounded by "arrayLimits".
	vector<int> asIndexListExcept(const IndexBoundingBox& other, const IndexBoundingBox& arrayLimits) const;

	Vector3_i nNodes() const
	{	return Vector3_i(iMax-iMin+1, jMax-jMin+1, kMax-kMin+1); }

	// Get number of nodes in the box including the nodes on the borders.
	int nNodesTotal() const
	{ return (iMax-iMin+1)*(jMax-jMin+1)*(kMax-kMin+1); }

	// Get the intersection between this box and another. Can be a "flat box" (plane).
	// If the boxes don't intersect, a box with only zeros is returned.
	IndexBoundingBox intersection(const IndexBoundingBox& other) const
	{
		IndexBoundingBox _intersection(max(this->iMin, other.iMin), min(this->iMax, other.iMax),
									   max(this->jMin, other.jMin), min(this->jMax, other.jMax),
									   max(this->kMin, other.kMin), min(this->kMax, other.kMax) );
		if (_intersection.iMax >= _intersection.iMin
		&&	_intersection.jMax >= _intersection.jMin
		&&	_intersection.kMax >= _intersection.kMin )
			return _intersection;
		else
			return IndexBoundingBox(0, 0, 0);
	}

	// Get node with the lowest indices
	Vector3_i getMinIndices() const
	{	return Vector3_i(iMin, jMin, kMin); }

	// Get node with the highest indices
	Vector3_i getMaxIndices() const
	{	return Vector3_i(iMax, jMax, kMax); }

	// Get box that's centered on the given node, and stretches 'radius' nodes in each direction.
	static IndexBoundingBox boxAroundNode(const Vector3_i& centerNode, int radius)
	{
		return IndexBoundingBox(centerNode.i-radius, centerNode.i+radius,
								centerNode.j-radius, centerNode.j+radius,
								centerNode.k-radius, centerNode.k+radius);
	}
};

// Get a box that represents the 8 nodes that surround an arbitrary point.
// If the point is outside the array limits, the box is still inside.
inline IndexBoundingBox getSurroundingNodesBox(const Vector3_d& point,
											   const Vector3_d& gridSpacings,
											   const Vector3_d& meshOriginOffset,
											   const IndexBoundingBox& arrayLimits)
{
	int iNext = static_cast<int>( ceil( point.x - meshOriginOffset.x ) / gridSpacings.x ) + arrayLimits.iMin;
	int jNext = static_cast<int>( ceil( point.y - meshOriginOffset.y ) / gridSpacings.y ) + arrayLimits.jMin;
	int kNext = static_cast<int>( ceil( point.z - meshOriginOffset.z ) / gridSpacings.z ) + arrayLimits.kMin;
	iNext = max<int>( iNext, arrayLimits.iMin+1 ); // Make sure that none of these are less than 1,
	jNext = max<int>( jNext, arrayLimits.jMin+1 ); // because that would make the lower index negative,
	kNext = max<int>( kNext, arrayLimits.kMin+1 ); // e.g., outside the mesh
	iNext = min<int>( iNext, arrayLimits.iMax ); // Also make sure the box is not outside
	jNext = min<int>( jNext, arrayLimits.jMax ); // on the max side, which could happen due
	kNext = min<int>( kNext, arrayLimits.kMax ); // to machine precision on image point.
IndexBoundingBox surroundingBox(iNext-1, iNext, jNext-1, jNext, kNext-1, kNext);
	return surroundingBox;
}

// Struct that represents a bounded region in space
struct SpaceBoundingBox
{
	double xMin, xMax;
	double yMin, yMax;
	double zMin, zMax;

	Vector3_d getMinPoint() const
	{	return Vector3_d(xMin, yMin, zMin); }

	Vector3_d getMaxPoint() const
	{	return Vector3_d(xMax, yMax, zMax); }
};


#endif /* SRC_BOUNDINGBOX_H_ */
