/*
 * BoundaryCondition.h
 *
 *  Created on: Oct 12, 2022
 *      Author: frederk
 */

#ifndef SRC_BOUNDARY_H_
#define SRC_BOUNDARY_H_

#include "includes_and_names.h"
#include "CustomExceptions.h"


enum class AxisOrientationEnum
{
	x, y, z
};

enum class EdgeIndexEnum
{
	min, max
};

class MeshEdgeBoundary;
typedef vector<unique_ptr<MeshEdgeBoundary>> EdgeBoundaryCollection;

// Base class for boundaries at the the edge of the Cartesian computational mesh.
// Suggested derived types: Inlet, outlet, walls, periodic, symmetry...
class MeshEdgeBoundary
{
public:
	MeshEdgeBoundary(AxisOrientationEnum normalAxis,
					 EdgeIndexEnum planeIndex);
	virtual ~MeshEdgeBoundary() = default;
	void identifyOwnedNodes(const EdgeBoundaryCollection& existingBoundaries);
	virtual void applyBoundaryCondition();
	const AxisOrientationEnum normalAxis;
	const uint planeIndex;
private:
	vector<uint> nodeIndices;
};

// Class to define inlet boundary condition:
class InletBoundary : public MeshEdgeBoundary
{
public:
	InletBoundary(AxisOrientationEnum normalAxis,
				  EdgeIndexEnum planeIndex,
				  double velocity);
	void applyBoundaryCondition() override;
private:
	double velocity;
};

// Class to define outlet boundary condition:
class OutletBoundary : public MeshEdgeBoundary
{
public:
	OutletBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex);
	void applyBoundaryCondition() override;
};

// Class to define a periodic boundary condition:
class PeriodicBoundary : public MeshEdgeBoundary
{
public:
	PeriodicBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex);
	void applyBoundaryCondition() override;
};

// Class to define a symmetry boundary condition:
class SymmetryBoundary : public MeshEdgeBoundary
{
public:
	SymmetryBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex);
	void applyBoundaryCondition() override;
};

class ImmersedBoundary;
typedef vector<unique_ptr<ImmersedBoundary>> ImmersedBoundaryCollection;

// Base class for immersed boundaries.
// Suggested derived types: Adiabatic wall and iso-thermal wall
// or different shapes.
class ImmersedBoundary
{
public:
	ImmersedBoundary();
	virtual ~ImmersedBoundary() = default;
	virtual void applyBoundaryCondition();
private:
	vector<uint> ghostNodeIndices;
};

// Class to define boundary conditions at an immersed cylinder:
class CylinderBody : public ImmersedBoundary
{
public:
	CylinderBody(sf::Vector2<double> centroidPosition, AxisOrientationEnum axis, double radius);
	void applyBoundaryCondition() override;
private:
	sf::Vector2<double> centroidPosition;
	AxisOrientationEnum axis;
	double radius;
};

// Class to define boundary conditions at an immersed sphere:
class SphereBody : public ImmersedBoundary
{
public:
	SphereBody(sf::Vector3<double> centerPosition, double radius);
	void applyBoundaryCondition() override;
private:
	sf::Vector3<double> centerPosition;
	double radius;
};

#endif /* SRC_BOUNDARY_H_ */




