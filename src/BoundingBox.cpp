/*
 * BoundingBox.cpp
 *
 *  Created on: May 25, 2023
 *      Author: frederk
 */


#include "BoundingBox.h"
#include "Node.h"

// Get a simple fixed size array with the 1D indices to the 8 corner nodes in the box
// Ordering: z(k) changes most often, x(i) changes most seldom.
Array8_u IndexBoundingBox::cornersAsIndexList(const IndexBoundingBox& arrayLimits) const
{
	Array8_u indices = { getIndex1D(iMin, jMin, kMin, arrayLimits),
						 getIndex1D(iMin, jMin, kMax, arrayLimits),
						 getIndex1D(iMin, jMax, kMin, arrayLimits),
						 getIndex1D(iMin, jMax, kMax, arrayLimits),
						 getIndex1D(iMax, jMin, kMin, arrayLimits),
						 getIndex1D(iMax, jMin, kMax, arrayLimits),
						 getIndex1D(iMax, jMax, kMin, arrayLimits),
						 getIndex1D(iMax, jMax, kMax, arrayLimits)
	};
	return indices;
}

// Get a vector to the 1D indices to all nodes.
// The 1D indices are with respect to an array bounded by "arrayLimits".
vector<int> IndexBoundingBox::allNodesAsIndexList(const IndexBoundingBox& arrayLimits) const
{
	vector<int> indexList;
	for(int i=iMin; i<=iMax; ++i)
		for(int j=jMin; j<=jMax; ++j)
			for(int k=kMin; k<=kMax; ++k)
				indexList.push_back( getIndex1D(i, j, k, arrayLimits) );
	return indexList;
}

// Get a vector to the 1D indices to all nodes, except those who are also entailed by another given box.
// The 1D indices are with respect to an array bounded by "arrayLimits".
vector<int> IndexBoundingBox::asIndexListExcept(const IndexBoundingBox& other, const IndexBoundingBox& arrayLimits) const
{
	vector<int> indexList;
	for(int i=iMin; i<=iMax; ++i)
		for(int j=jMin; j<=jMax; ++j)
			for(int k=kMin; k<=kMax; ++k)
				if( 	i < other.iMin || i > other.iMax
					||	j < other.jMin || j > other.jMax
					||	k < other.kMin || k > other.kMax )
					indexList.push_back( getIndex1D(i, j, k, arrayLimits) );
	return indexList;
}

