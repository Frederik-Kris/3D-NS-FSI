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
	double length() { return sqrt(pow(x,2)+pow(y,2)+pow(z,2)); }
	inline Vector3_d operator+(Vector3_d other) const { return Vector3_d(x+other.x, y+other.y, z+other.z); }
	inline Vector3_d operator-(Vector3_d other) const { return Vector3_d(x-other.x, y-other.y, z-other.z); }
	inline Vector3_d operator*(double factor) const { return Vector3_d(x*factor, y*factor, z*factor); }
	inline Vector3_d operator/(double divisor) const { return Vector3_d(x/divisor, y/divisor, z/divisor); }
};

inline Vector3_d operator*(double factor, Vector3_d vector) { return vector*factor; }

// Simple 3D vector with size_t (unsigned) components
struct Vector3_u
{
	size_t i, j, k;
	Vector3_u(size_t i, size_t j, size_t k) : i{i}, j{j}, k{k} {}
};

// Given the 1D index to a node, get the 3D indices
inline Vector3_u getIndices3D(size_t index1D, const Vector3_u& nNodes)
{
	size_t i = index1D / (nNodes.j * nNodes.k);
	size_t j = index1D % (nNodes.j * nNodes.k) / nNodes.k;
	size_t k = index1D % (nNodes.j * nNodes.k) % nNodes.k;
	return Vector3_u(i,j,k);
}

// Given the 3D indices to a node, get the 1D index
inline size_t getIndex1D(size_t i, size_t j, size_t k, const Vector3_u& nNodes)
{
	return i * nNodes.j * nNodes.k + j * nNodes.k + k;
}

// Given the 3D indices to a node, get the 1D index
inline size_t getIndex1D(const Vector3_u& indices, const Vector3_u& nNodes)
{ return getIndex1D(indices.i, indices.j, indices.k, nNodes); }

// Get the position (coordinates) of the node with the given indices
inline Vector3_d getNodePosition(size_t i, size_t j, size_t k,
								 const Vector3_d& gridSpacings,
								 const Vector3_d& meshOriginOffset)
{
	double x = ( static_cast<int>(i) - meshOriginOffset.x ) * gridSpacings.x ;
	double y = ( static_cast<int>(j) - meshOriginOffset.y ) * gridSpacings.y ;
	double z = ( static_cast<int>(k) - meshOriginOffset.z ) * gridSpacings.z ;
	return Vector3_d(x, y, z);
}

// Get the position (coordinates) of the node with the given indices
inline Vector3_d getNodePosition(const Vector3_u& indices,
								 const Vector3_d& gridSpacings,
								 const Vector3_d& meshOriginOffset)
{ return getNodePosition(indices.i, indices.j, indices.k, gridSpacings, meshOriginOffset); }

// Struct that represents a bounded region of the mesh, defined by lowest and highest indices
struct IndexBoundingBox
{
	size_t iMin, iMax;	// Indices in x-direction
	size_t jMin, jMax;	// Indices in y-direction
	size_t kMin, kMax;	// Indices in z-direction

	// Constructor. Takes in the upper index limits, and sets the lower limits to zero.
	IndexBoundingBox(size_t iMax, size_t jMax, size_t kMax)
	: iMin{0}, iMax{iMax},
	  jMin{0}, jMax{jMax},
	  kMin{0}, kMax{kMax}
	{}

	// Default constructor. Leaves data uninitialized.
	IndexBoundingBox() = default;

	// Get a simple fixed size array with the 1D indices to the 8 corner nodes in the box
	// Ordering: z(k) changes most often, x(i) changes most seldom.
	Array8_u asIndexList(const Vector3_u& nNodes) const
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
};

// Get a box that represents the 8 nodes that surround an arbitrary point.
// If the point is outside the mesh edges, the box is still inside.
inline IndexBoundingBox getSurroundingNodesBox(const Vector3_d& point,
											   const Vector3_d& gridSpacings,
											   const Vector3_d& meshOriginOffset,
											   const Vector3_u& nMeshNodes)
{
	size_t iNext = static_cast<size_t>( ceil( point.x / gridSpacings.x ) + meshOriginOffset.x );
	size_t jNext = static_cast<size_t>( ceil( point.y / gridSpacings.y ) + meshOriginOffset.y );
	size_t kNext = static_cast<size_t>( ceil( point.z / gridSpacings.z ) + meshOriginOffset.z );
	iNext = max<size_t>( iNext, 1 ); // Make sure that none of these are less than 1,
	jNext = max<size_t>( jNext, 1 ); // because that would make the lower index negative,
	kNext = max<size_t>( kNext, 1 ); // e.g., outside the mesh
	iNext = min<size_t>( iNext, nMeshNodes.i ); // Also make sure the box is not outside
	jNext = min<size_t>( jNext, nMeshNodes.j ); // on the max side, which could happen due
	kNext = min<size_t>( kNext, nMeshNodes.k ); // to machine precision on image point.
	IndexBoundingBox surroundingBox(iNext, jNext, kNext);
	surroundingBox.iMin = iNext - 1;
	surroundingBox.jMin = jNext - 1;
	surroundingBox.kMin = kNext - 1;
	return surroundingBox;
}


#endif /* SRC_SMALLVECTORS_H_ */
