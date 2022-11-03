/*
 * SmallVectors.h
 *
 *  Created on: Oct 19, 2022
 *      Author: frederk
 */

#ifndef SRC_SMALLVECTORS_H_
#define SRC_SMALLVECTORS_H_

#include "includes_and_names.h"

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

struct Vector3_u
{
	size_t i, j, k;
	Vector3_u(size_t i, size_t j, size_t k) : i{i}, j{j}, k{k} {}
};


#endif /* SRC_SMALLVECTORS_H_ */
