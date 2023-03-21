/*
 * Array3D_d.h
 *
 *  Created on: Apr 8, 2021
 *      Author: frederik
 */

#ifndef SRC_ARRAY3D_H_
#define SRC_ARRAY3D_H_



#include "includes_and_names.h"
#include "SmallVectors.h"
#include "Node.h"


// Simple 3D array template.
// Data is stored internally as 1D std::vector.
template <typename T>
class Array3D
{
public:

	// Constructor. Takes the index limits in all three directions, which decides the size of the array.
	// Default-inserts elements in the array.
	Array3D(const IndexBoundingBox& indexLimits) :
		indexLimits(indexLimits),
		length( indexLimits.iMax - indexLimits.iMin + 1 ),
		width ( indexLimits.jMax - indexLimits.jMin + 1 ),
		height( indexLimits.kMax - indexLimits.kMin + 1 )
	{}

	// Constructor. Takes number of nodes in each direction.
	// Minimum index in all directions is zero
	// Default-inserts the elements in the data-vector.
	Array3D(int length, int width, int height) :
		Array3D( IndexBoundingBox(length-1, width-1, height-1) )
	{}

	// Default constructor. Sets sizes to zero and default-constructs data-vector.
	Array3D() :
		Array3D(0,0,0)
	{}

	// Constructor. Sets index limits (which implicitly gives size) AND initializes all nodes with given value.
	Array3D(const IndexBoundingBox& indexLimits, const T& value) :
		Array3D(indexLimits),
		data(length * width * height, value)
	{}

	// Constructor. Sets number of nodes in each direction AND initializes all nodes with given value.
	// Lower index limit in all directions is zero.
	Array3D(int length, int width, int height, const T& value) :
		Array3D( IndexBoundingBox(length-1, width-1, height-1), value )
	{}

	// Get reference to a node using 3D indices
	inline T& operator()(int i, int j, int k)
	{ return data[ (i-indexLimits.iMin)*width*height + (j-indexLimits.jMin)*height + (k-indexLimits.kMin) ]; }

	// Get value of a node using 3D indices
	inline T operator()(int i, int j, int k) const
	{ return data[ (i-indexLimits.iMin)*width*height + (j-indexLimits.jMin)*height + (k-indexLimits.kMin) ]; }

	// Get reference to a node using 3D indices in Vector3_u
	inline T& operator()(const Vector3_i& indices)
	{ return data[ (indices.i-indexLimits.iMin)*width*height + (indices.j-indexLimits.jMin)*height + (indices.k-indexLimits.kMin) ]; }

	// Get value of a node using 3D indices in Vector3_u
	inline T operator()(const Vector3_i& indices) const
	{ return data[ (indices.i-indexLimits.iMin)*width*height + (indices.j-indexLimits.jMin)*height + (indices.k-indexLimits.kMin) ]; }

	// Get reference to a node using one index
	inline T& operator()(int i)
	{ return data[i]; }

	// Get value of a node using one index
	inline T operator()(int i) const
	{ return data[i]; }

	// Get reference to a node using 3D indices
	// ->at(i,j,k) is more readable than ->operator()(i,j,k), or (*array)(i,j,k)
	inline T& at(int i, int j, int k)
	{ return data[ (i-indexLimits.iMin)*width*height + (j-indexLimits.jMin)*height + (k-indexLimits.kMin) ]; }

	// Get value of a node using 3D indices
	// ->at(i,j,k) is more readable than ->operator()(i,j,k), or (*array)(i,j,k)
	inline T at(int i, int j, int k) const
	{ return data[ (i-indexLimits.iMin)*width*height + (j-indexLimits.jMin)*height + (k-indexLimits.kMin) ]; }

	// Get reference to a node using one index
	// ->at(i) is more readable than ->operator()(i), or (*array)(i)
	inline T& at(int i)
	{ return data[i]; }

	// Get value of a node using one index
	// ->at(i) is more readable than ->operator()(i), or (*array)(i)
	inline T at(int i) const
	{ return data[i]; }

	// Get portion of the array as new array. Index limits are relative to old array.
	Array3D<T> subArray(const IndexBoundingBox& subDomain) const
	{
		Array3D<T> _subArray(subDomain);
		for(int i=subDomain.iMin; i<=subDomain.iMax; ++i)
			for(int j=subDomain.jMin; j<=subDomain.jMax; ++j)
				for(int k=subDomain.kMin; k<=subDomain.kMax; ++k)
					_subArray(i,j,k) = this->at(i, j, k);
	}

	// Set a new size for the array, using index limits in each direction.
	// New underlying array (vector) is constructed. Old data is lost.
	void setSize(const IndexBoundingBox& newIndexLimits)
	{
		indexLimits = newIndexLimits;
		length = indexLimits.iMax - indexLimits.iMin + 1 ;
		width  = indexLimits.jMax - indexLimits.jMin + 1 ;
		height = indexLimits.kMax - indexLimits.kMin + 1 ;
		data = vector<T>(length*width*height);
	}

	// Set a new size for the array, in no. of nodes in each direction.
	// New minimum index in all directions will be zero.
	// New underlying array (vector) is constructed. Old data is lost.
	void setSize(int newLength, int newWidth, int newHeight)
	{ setSize( IndexBoundingBox(newLength-1, newWidth-1, newHeight-1) ); }

	// Set all the nodes in the array to 'value'
	void setAll(const T& value)
	{ data.assign(data.size(), value); }

	// Const iterator to the element with 1D index 0:
	std::vector::const_iterator begin() { return data.begin(); }

	// Const iterator to the past-the-end-element of the underlying 1D vector:
	std::vector::const_iterator end() { return data.end(); }

	// Swap contents of the array with another array, using move-semantics (no copy).
	void dataSwap(Array3D& other)
	{ std::swap(data, other.data); }

	// Check if all of the elements are finite, e.g., not NaN or Inf
	bool allFinite() const
	{
		bool allFinite{true};
		for(T scalar : data)
			allFinite = allFinite && std::isfinite(scalar);
		return allFinite;
	}

	const IndexBoundingBox& getIndexLimits() const
	{ return indexLimits; }

	int getLength() const
	{ return length; }

	int getWidth() const
	{ return width; }

	int getHeight() const
	{ return height; }

	int getSize() const
	{ return data.size(); }

	T getMax() const
	{ return *std::max_element( data.begin(), data.end() ); }

	T getMin() const
	{ return *std::min_element( data.begin(), data.end() ); }

private:
	IndexBoundingBox indexLimits;
	int length, width, height; // Number of nodes in the 3 directions
	vector<T> data;
};

typedef Array3D<double> Array3D_d;




#endif /* SRC_ARRAY3D_H_ */
