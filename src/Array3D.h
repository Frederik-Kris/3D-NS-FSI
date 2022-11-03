/*
 * Array3D_d.h
 *
 *  Created on: Apr 8, 2021
 *      Author: frederik
 */

#ifndef SRC_ARRAY3D_H_
#define SRC_ARRAY3D_H_



#include "includes_and_names.h"



// Simple 3D array type for double precision numbers.
class Array3D_d
{
public:
	// Constructor. Sets the size members and initializes the data-vector with correct size.
	// Does NOT initialize the elements in the data-vector.
	Array3D_d(size_t length, size_t width, size_t height) :
    length(length),
	width(width),
    height(height),
    data(length * width * height)
	{}

	// Get reference to a node using 3D indices
	inline double& operator()(size_t x, size_t y, size_t z)
	{ return data[x * width * height + y * height + z]; }

	// Get value of a node using 3D indices
	inline double operator()(size_t x, size_t y, size_t z) const
	{ return data[x * width * height + y * height + z]; }

	// Get reference to a node using 3D indices in Vector3_u
	inline double& operator()(Vector3_u indices)
	{ return data[indices.i * width * height + indices.j * height + indices.k]; }

	// Get reference to a node using 3D indices in Vector3_u
	inline double operator()(Vector3_u indices) const
	{ return data[indices.i * width * height + indices.j * height + indices.k]; }

	// Get reference to a node using one index
	inline double& operator()(size_t i)
	{ return data[i]; }

	// Get value of a node using one index
	inline double operator()(size_t i) const
	{ return data[i]; }

	// Get reference to a node using one index
	// ->at(i) is more readable than ->operator(i), or (*array)(i)
	inline double& at(size_t i)
	{ return data[i]; }

	// Get value of a node using one index
	// ->at(i) is more readable than ->operator(i), or (*array)(i)
	inline double at(size_t i) const
	{ return data[i]; }

	// Set all the nodes in the array to value 'd'
	void setAll(double d)
	{ data.assign(data.size(), d); }

	// Swap contents of the array with another array, using move-semantics (no copy).
	void dataSwap(Array3D_d& other)
	{ std::swap(data, other.data); }

	// Addition operator overload. Sums elementwise. Intermediate array is created.
	Array3D_d operator+(const Array3D_d& other) const
	{
		Array3D_d sumResult(length, width, height);
		for(size_t i{0}; i<data.size(); ++i)
			sumResult.data[i] = this->data[i] + other.data[i];
		return sumResult;
	}

	// Assign each node as elementwise sum from 2 other arrays. NO intermediate array created.
	void assignSum(const Array3D_d& A1, const Array3D_d& A2)
	{
		for(size_t i{0}; i<data.size(); ++i)
			this->data[i] = A1.data[i] + A2.data[i];
	}

	size_t getLength() const
	{ return length; }

	size_t getWidth() const
	{ return width; }

	size_t getHeight() const
	{ return height; }

	size_t getSize() const
	{ return data.size(); }

	double getMax() const
	{ return *std::max_element( data.begin(), data.end() ); }

	double getMin() const
	{ return *std::min_element( data.begin(), data.end() ); }

private:
	size_t length, width, height;
	vector<double> data;
};

enum class NodeTypeEnum
{
	FluidRegular, FluidEdge, Solid, Ghost
};

// Simple 3D array type for NodeTypeEnum
class Array3D_nodeType
{
public:
	// Constructor. Sets the size members and initializes the data-vector with correct size.
	// Does NOT initialize the elements in the data-vector.
	Array3D_nodeType(size_t length, size_t width, size_t height) :
    length(length),
	width(width),
    height(height),
    data(length * width * height)
	{}

	// Get reference to a node using 3D indices
	inline NodeTypeEnum& operator()(size_t x, size_t y, size_t z)
	{ return data[x * width * height + y * height + z]; }

	// Get value of a node using 3D indices
	inline NodeTypeEnum operator()(size_t x, size_t y, size_t z) const
	{ return data[x * width * height + y * height + z]; }

	// Get reference to a node using 3D indices in vector
	inline NodeTypeEnum& operator()(Vector3_u indices)
	{ return data[indices.i * width * height + indices.j * height + indices.k]; }

	// Get value of a node using 3D indices in vector
	inline NodeTypeEnum operator()(Vector3_u indices) const
	{ return data[indices.i * width * height + indices.j * height + indices.k]; }

	// Get reference to a node using one index
	inline NodeTypeEnum& operator()(size_t i)
	{ return data[i]; }

	// Get value of a node using one index
	inline NodeTypeEnum operator()(size_t i) const
	{ return data[i]; }

	// Set all the nodes in the array to value 'type'
	void setAll(NodeTypeEnum type)
	{ data.assign(data.size(), type); }

	// Swap contents of the array with another array, using move-semantics (no copy).
	void dataSwap(Array3D_nodeType& other)
	{ std::swap(data, other.data); }

private:
	size_t length, width, height;
	vector<NodeTypeEnum> data;
};




#endif /* SRC_ARRAY3D_H_ */
