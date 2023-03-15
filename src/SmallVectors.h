/*
 * SmallVectors.h
 *
 *  Created on: Oct 19, 2022
 *      Author: frederk
 */

#ifndef SRC_SMALLVECTORS_H_
#define SRC_SMALLVECTORS_H_


#include "includes_and_names.h"


// Simple 3D vector with double precision components. Supports basic arithmetic operators.
struct Vector3_d
{
	double x, y, z;

	Vector3_d(double x, double y, double z) : x{x}, y{y}, z{z} {}

	Vector3_d() = default;

	// A.k.a norm, by 3D Pythagoras, squareroot of the sum of squared components
	double length() { return sqrt(pow(x,2)+pow(y,2)+pow(z,2)); }

	inline Vector3_d operator+(Vector3_d other) const { return Vector3_d(x+other.x, y+other.y, z+other.z); }
	inline Vector3_d operator-(Vector3_d other) const { return Vector3_d(x-other.x, y-other.y, z-other.z); }
	inline Vector3_d operator*(double factor) const { return Vector3_d(x*factor, y*factor, z*factor); }
	inline Vector3_d operator/(double divisor) const { return Vector3_d(x/divisor, y/divisor, z/divisor); }
};

inline Vector3_d operator*(double factor, Vector3_d vector) { return vector*factor; }

// Simple 3D vector with integer components
struct Vector3_i
{
	int i, j, k;
	Vector3_i(int i, int j, int k) : i{i}, j{j}, k{k} {}
};

// Given the 1D index to a node, get the 3D indices
inline Vector3_i getIndices3D(int index1D, const Vector3_i& nNodes)
{
	int i = index1D / (nNodes.j * nNodes.k);
	int j = index1D % (nNodes.j * nNodes.k) / nNodes.k;
	int k = index1D % (nNodes.j * nNodes.k) % nNodes.k;
	return Vector3_i(i,j,k);
}

// Given the 3D indices to a node, get the 1D index
inline int getIndex1D(int i, int j, int k, const Vector3_i& nNodes)
{
	return i * nNodes.j * nNodes.k + j * nNodes.k + k;
}

// Given the 3D indices to a node, get the 1D index
inline int getIndex1D(const Vector3_i& indices, const Vector3_i& nNodes)
{ return getIndex1D(indices.i, indices.j, indices.k, nNodes); }

// Get the position (coordinates) of the node with the given indices
inline Vector3_d getNodePosition(int i, int j, int k,
								 const Vector3_d& gridSpacings,
								 const Vector3_d& meshOriginOffset)
{
	double x = ( static_cast<int>(i) - meshOriginOffset.x ) * gridSpacings.x ;
	double y = ( static_cast<int>(j) - meshOriginOffset.y ) * gridSpacings.y ;
	double z = ( static_cast<int>(k) - meshOriginOffset.z ) * gridSpacings.z ;
	return Vector3_d(x, y, z);
}

// Get the position (coordinates) of the node with the given indices
inline Vector3_d getNodePosition(const Vector3_i& indices,
								 const Vector3_d& gridSpacings,
								 const Vector3_d& meshOriginOffset)
{ return getNodePosition(indices.i, indices.j, indices.k, gridSpacings, meshOriginOffset); }

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

	// Default constructor. Leaves data uninitialized.
	IndexBoundingBox() = default;

	// Get a simple fixed size array with the 1D indices to the 8 corner nodes in the box
	// Ordering: z(k) changes most often, x(i) changes most seldom.
	Array8_u asIndexList(const Vector3_i& nNodes) const
	{
		Array8_u indices = { getIndex1D(iMin, jMin, kMin, nNodes),
							 getIndex1D(iMin, jMin, kMax, nNodes),
							 getIndex1D(iMin, jMax, kMin, nNodes),
							 getIndex1D(iMin, jMax, kMax, nNodes),
							 getIndex1D(iMax, jMin, kMin, nNodes),
							 getIndex1D(iMax, jMin, kMax, nNodes),
							 getIndex1D(iMax, jMax, kMin, nNodes),
							 getIndex1D(iMax, jMax, kMax, nNodes)
		};
		return indices;
	}

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

	// Get box that's centered on the given node, and stretches 'radius' nodes in each direction.
	static IndexBoundingBox boxAroundNode(const Vector3_i& centerNode, int radius)
	{
		return IndexBoundingBox(centerNode.i-radius, centerNode.i+radius,
								centerNode.j-radius, centerNode.j+radius,
								centerNode.k-radius, centerNode.k+radius);
	}
};

// Get a box that represents the 8 nodes that surround an arbitrary point.
// If the point is outside the mesh edges, the box is still inside.
inline IndexBoundingBox getSurroundingNodesBox(const Vector3_d& point,
											   const Vector3_d& gridSpacings,
											   const Vector3_d& meshOriginOffset,
											   const Vector3_i& nMeshNodes)
{
	int iNext = static_cast<int>( ceil( point.x / gridSpacings.x ) + meshOriginOffset.x );
	int jNext = static_cast<int>( ceil( point.y / gridSpacings.y ) + meshOriginOffset.y );
	int kNext = static_cast<int>( ceil( point.z / gridSpacings.z ) + meshOriginOffset.z );
	iNext = max<int>( iNext, 1 ); // Make sure that none of these are less than 1,
	jNext = max<int>( jNext, 1 ); // because that would make the lower index negative,
	kNext = max<int>( kNext, 1 ); // e.g., outside the mesh
	iNext = min<int>( iNext, nMeshNodes.i ); // Also make sure the box is not outside
	jNext = min<int>( jNext, nMeshNodes.j ); // on the max side, which could happen due
	kNext = min<int>( kNext, nMeshNodes.k ); // to machine precision on image point.
	IndexBoundingBox surroundingBox(iNext, jNext, kNext);
	surroundingBox.iMin = iNext - 1;
	surroundingBox.jMin = jNext - 1;
	surroundingBox.kMin = kNext - 1;
	return surroundingBox;
}

// Struct that represents a bounded region in space
struct SpaceBoundingBox
{
	double xMin, xMax;
	double yMin, yMax;
	double zMin, zMax;
};


#endif /* SRC_SMALLVECTORS_H_ */
