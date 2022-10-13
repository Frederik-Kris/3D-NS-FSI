/*
 * BoundaryCondition.h
 *
 *  Created on: Oct 12, 2022
 *      Author: frederk
 */

#ifndef SRC_BOUNDARY_H_
#define SRC_BOUNDARY_H_

#include "includes_and_names.h"

enum class AxisOrientationEnum
{
	x, y, z
};

// Base class for boundaries at the the edge of the Cartesian computational mesh.
// Suggested derived types: Inlet, outlet, walls, periodic, symmetry...
class MeshEdgeBoundary
{
public:
	MeshEdgeBoundary(AxisOrientationEnum normalAxis, uint planeIndex) :
					   normalAxis{normalAxis},
					   planeIndex{planeIndex}
	{}
	virtual ~MeshEdgeBoundary();
	virtual void applyBoundaryCondition();
	const AxisOrientationEnum normalAxis;
	const uint planeIndex;
};

// Class to define inlet boundary condition:
class InletBoundary : public MeshEdgeBoundary
{
public:
	InletBoundary(AxisOrientationEnum normalAxis,
				  uint planeIndex,
				  double velocity) :
					  MeshEdgeBoundary(normalAxis, planeIndex),
					  velocity{velocity}
	{}
	void applyBoundaryCondition() override;
private:
	double velocity;
};

// Class to define outlet boundary condition:
class OutletBoundary : public MeshEdgeBoundary
{
public:
	OutletBoundary(AxisOrientationEnum normalAxis,
				   uint planeIndex) :
					   MeshEdgeBoundary(normalAxis, planeIndex)
	{}
	void applyBoundaryCondition() override;
};

// Class to define a periodic boundary condition:
class PeriodicBoundary : public MeshEdgeBoundary
{
public:
	PeriodicBoundary(AxisOrientationEnum normalAxis,
					 uint planeIndex) :
						 MeshEdgeBoundary(normalAxis, planeIndex)
	{}
	void applyBoundaryCondition() override;
};

// Class to define a symmetry boundary condition:
class SymmetryBoundary : public MeshEdgeBoundary
{
public:
	SymmetryBoundary(AxisOrientationEnum normalAxis,
					 uint planeIndex) :
						 MeshEdgeBoundary(normalAxis, planeIndex)
	{}
	void applyBoundaryCondition() override;
};

// Base class for immersed boundaries.
// Suggested derived types: Adiabatic wall and iso-thermal wall
// or different shapes.
class ImmersedBoundary
{
public:
	ImmersedBoundary();
	virtual ~ImmersedBoundary();
	virtual void applyBoundaryCondition();
private:
	vector<uint> ghostNodeIndices;
};

// Class to define boundary conditions at an immersed cylinder:
class CylinderBody : public ImmersedBoundary
{
public:
	CylinderBody(sf::Vector2<double> centroidPosition,
				 AxisOrientationEnum axis,
				 double radius) :
					 centroidPosition{centroidPosition},
					 axis{axis},
					 radius{radius}
	{}
private:
	sf::Vector2<double> centroidPosition;
	AxisOrientationEnum axis;
	double radius;
};

// Class to define boundary conditions at an immersed sphere:
class SphereBody : public ImmersedBoundary
{
public:
	SphereBody(sf::Vector3<double> centerPosition,
			   double radius) :
				   centerPosition{centerPosition},
				   radius{radius}
	{}
private:
	sf::Vector3<double> centerPosition;
	double radius;
};

#endif /* SRC_BOUNDARY_H_ */




