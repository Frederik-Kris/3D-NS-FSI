/*
 * SmallVectors.h
 *
 *  Created on: Oct 19, 2022
 *      Author: frederk
 */

#ifndef SRC_SMALLVECTORS_H_
#define SRC_SMALLVECTORS_H_


#include "includes_and_names.h"


// Simple 3D vector with integer components
struct Vector3_i
{
	int i, j, k;
	Vector3_i(int i, int j, int k) : i{i}, j{j}, k{k} {}
	Vector3_i() : Vector3_i(0,0,0) {}

	Vector3_i operator+(const Vector3_i& other) const { return Vector3_i(i+other.i, j+other.j, k+other.k); }
	Vector3_i operator-(const Vector3_i& other) const { return Vector3_i(i-other.i, j-other.j, k-other.k); }
	Vector3_i operator%(const Vector3_i& other)	const { return Vector3_i(i%other.i, j%other.j, k%other.k); }
	Vector3_i operator+(int scalar)		 const { return Vector3_i(i+scalar, j+scalar, k+scalar); }
	Vector3_i operator-(int scalar)		 const { return Vector3_i(i-scalar, j-scalar, k-scalar); }
	Vector3_i operator*(int factor) 	 const { return Vector3_i(i*factor, j*factor, k*factor); }
	Vector3_i operator%(int divisor)	 const { return Vector3_i(i%divisor, j%divisor, k%divisor); }
};

inline Vector3_i operator*(int factor, Vector3_i vector) { return vector*factor; }

// Simple 3D vector with double precision components. Supports basic arithmetic operators.
struct Vector3_d
{
	double x, y, z;

	Vector3_d(double x, double y, double z) : x{x}, y{y}, z{z} {}

	Vector3_d(const Vector3_i& intVector) :
		x{static_cast<double>(intVector.i)},
		y{static_cast<double>(intVector.j)},
		z{static_cast<double>(intVector.k)}
		{}

	Vector3_d() : Vector3_d(0,0,0) {}

	// A.k.a norm, by 3D Pythagoras, squareroot of the sum of squared components
	double length() { return sqrt(pow(x,2)+pow(y,2)+pow(z,2)); }

	Vector3_d operator+(Vector3_d other) const { return Vector3_d(x+other.x, y+other.y, z+other.z); }
	Vector3_d operator-(Vector3_d other) const { return Vector3_d(x-other.x, y-other.y, z-other.z); }
	Vector3_d operator*(double factor) const { return Vector3_d(x*factor, y*factor, z*factor); }
	Vector3_d operator/(double divisor) const { return Vector3_d(x/divisor, y/divisor, z/divisor); }
};

inline Vector3_d operator*(double factor, Vector3_d vector) { return vector*factor; }


#endif /* SRC_SMALLVECTORS_H_ */






