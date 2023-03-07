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


// Simple 3D array type for double precision numbers.
// Data is stored internally as 1D std::vector.
class Array3D_d
{
public:

	// Constructor. Sets the size members and initializes the data-vector with correct size.
	// Does NOT initialize the elements in the data-vector.
	Array3D_d(int length, int width, int height) :
    length(length),
	width(width),
    height(height),
    data(length * width * height)
	{}

	// Default constructor. Sets sizes to zero and default-constructs data-vector.
	Array3D_d() :
    length(0),
	width(0),
    height(0)
	{}

	// Get reference to a node using 3D indices
	inline double& operator()(int x, int y, int z)
	{ return data[x * width * height + y * height + z]; }

	// Get value of a node using 3D indices
	inline double operator()(int x, int y, int z) const
	{ return data[x * width * height + y * height + z]; }

	// Get reference to a node using 3D indices in Vector3_u
	inline double& operator()(Vector3_i indices)
	{ return data[indices.i * width * height + indices.j * height + indices.k]; }

	// Get value of a node using 3D indices in Vector3_u
	inline double operator()(Vector3_i indices) const
	{ return data[indices.i * width * height + indices.j * height + indices.k]; }

	// Get reference to a node using one index
	inline double& operator()(int i)
	{ return data[i]; }

	// Get value of a node using one index
	inline double operator()(int i) const
	{ return data[i]; }

	// Get reference to a node using 3D indices
	// ->at(i,j,k) is more readable than ->operator()(i,j,k), or (*array)(i,j,k)
	inline double& at(int x, int y, int z)
	{ return data[x * width * height + y * height + z]; }

	// Get value of a node using 3D indices
	// ->at(i,j,k) is more readable than ->operator()(i,j,k), or (*array)(i,j,k)
	inline double at(int x, int y, int z) const
	{ return data[x * width * height + y * height + z]; }

	// Get reference to a node using one index
	// ->at(i) is more readable than ->operator()(i), or (*array)(i)
	inline double& at(int i)
	{ return data[i]; }

	// Get value of a node using one index
	// ->at(i) is more readable than ->operator()(i), or (*array)(i)
	inline double at(int i) const
	{ return data[i]; }

	// Set all the nodes in the array to value 'd'
	void setAll(double d)
	{ data.assign(data.size(), d); }

	// Swap contents of the array with another array, using move-semantics (no copy).
	void dataSwap(Array3D_d& other)
	{ std::swap(data, other.data); }

	// Check if all of the elements are finite, e.g., not NaN or Inf
	bool allFinite() const
	{
		bool allFinite{true};
		for(double d : data)
			allFinite = allFinite && std::isfinite(d);
		return allFinite;
	}

	int getLength() const
	{ return length; }

	int getWidth() const
	{ return width; }

	int getHeight() const
	{ return height; }

	int getSize() const
	{ return data.size(); }

	double getMax() const
	{ return *std::max_element( data.begin(), data.end() ); }

	double getMin() const
	{ return *std::min_element( data.begin(), data.end() ); }

private:
	int length, width, height; // Number of nodes in the 3 directions
	vector<double> data;
};

// Simple 3D array type for NodeTypeEnum
// Data is stored internally as 1D std::vector
class Array3D_nodeType
{
public:
	// Constructor. Sets the size members and initializes the data-vector with correct size.
	// Does NOT initialize the elements in the data-vector.
	Array3D_nodeType(int length, int width, int height) :
    length(length),
	width(width),
    height(height),
    data(length * width * height)
	{}

	// Get reference to a node using 3D indices
	inline NodeTypeEnum& operator()(int x, int y, int z)
	{ return data[x * width * height + y * height + z]; }

	// Get value of a node using 3D indices
	inline NodeTypeEnum operator()(int x, int y, int z) const
	{ return data[x * width * height + y * height + z]; }

	// Get reference to a node using 3D indices in vector
	inline NodeTypeEnum& operator()(Vector3_i indices)
	{ return data[indices.i * width * height + indices.j * height + indices.k]; }

	// Get value of a node using 3D indices in vector
	inline NodeTypeEnum operator()(Vector3_i indices) const
	{ return data[indices.i * width * height + indices.j * height + indices.k]; }

	// Get reference to a node using one index
	inline NodeTypeEnum& operator()(int i)
	{ return data[i]; }

	// Get value of a node using one index
	inline NodeTypeEnum operator()(int i) const
	{ return data[i]; }

	// Set all the nodes in the array to value 'type'
	void setAll(NodeTypeEnum type)
	{ data.assign(data.size(), type); }

	// Swap contents of the array with another array, using move-semantics (no copy).
	void dataSwap(Array3D_nodeType& other)
	{ std::swap(data, other.data); }

private:
	int length, width, height;
	vector<NodeTypeEnum> data;
};




#endif /* SRC_ARRAY3D_H_ */
