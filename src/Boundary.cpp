/*
 * Boundary.cpp
 *
 *  Created on: Oct 13, 2022
 *      Author: frederk
 */

#include "Boundary.h"

MeshEdgeBoundary::MeshEdgeBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex)
: normalAxis{normalAxis},
  planeIndex{planeIndex}
{}

void MeshEdgeBoundary::identifyOwnedNodes(const EdgeBoundaryCollection& existingBoundaries)
{
	// DETTE MÅ FLYTTES TIL MESH CLASS?
	uint iMin{0}, iMax, jMin, jMax, kMin, kMax;
	switch (normalAxis)
	{
	case AxisOrientationEnum::x:
		iMin = iMax = planeIndex;
		for(MeshEdgeBoundary existing : existingBoundaries)
		{
			if(existing.normalAxis == AxisOrientationEnum::y)
				if(existing.planeIndex == EdgeIndexEnum::min)
					jMin = 1;
		}
		break;
	case AxisOrientationEnum::y:

		break;
	case AxisOrientationEnum::z:

		break;
	default:
		throw UnexpectedEnumValueException<AxisOrientationEnum>(normalAxis);
	}
}

InletBoundary::InletBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex, double velocity)
: MeshEdgeBoundary(normalAxis, planeIndex),
  velocity{velocity}
{}

void InletBoundary::applyBoundaryCondition()
{

}

OutletBoundary::OutletBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex)
: MeshEdgeBoundary(normalAxis, planeIndex)
{}

void OutletBoundary::applyBoundaryCondition()
{

}

PeriodicBoundary::PeriodicBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex)
: MeshEdgeBoundary(normalAxis, planeIndex)
{}

void PeriodicBoundary::applyBoundaryCondition()
{

}

SymmetryBoundary::SymmetryBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex)
: MeshEdgeBoundary(normalAxis, planeIndex)
{}

void SymmetryBoundary::applyBoundaryCondition()
{

}


ImmersedBoundary::ImmersedBoundary()
{

}

CylinderBody::CylinderBody(sf::Vector2<double> centroidPosition, AxisOrientationEnum axis, double radius) :
centroidPosition{centroidPosition},
axis{axis},
radius{radius}
{}

void CylinderBody::applyBoundaryCondition()
{

}

SphereBody::SphereBody(sf::Vector3<double> centerPosition, double radius) :
centerPosition{centerPosition},
radius{radius}
{}

void SphereBody::applyBoundaryCondition()
{

}





