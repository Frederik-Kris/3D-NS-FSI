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

void MeshEdgeBoundary::identifyOwnedNodes(IndexBoundingBox& unclaimedNodes, const IndexBoundingBox meshSize)
{
	IndexBoundingBox indexBounds = unclaimedNodes;
	switch (normalAxis)
	{
	case AxisOrientationEnum::x:
		if(planeIndex == EdgeIndexEnum::min)
		{
			indexBounds.iMax = indexBounds.iMin;
			unclaimedNodes.iMin++;
		}
		else if(planeIndex == EdgeIndexEnum::max)
		{
			indexBounds.iMin = indexBounds.iMax;
			unclaimedNodes.iMax--;
		}
		else throw UnexpectedEnumValueException<EdgeIndexEnum>(planeIndex);
		break;
	case AxisOrientationEnum::y:
		if(planeIndex == EdgeIndexEnum::min)
		{
			indexBounds.jMax = indexBounds.jMin;
			unclaimedNodes.jMin++;
		}
		else if(planeIndex == EdgeIndexEnum::max)
		{
			indexBounds.jMin = indexBounds.jMax;
			unclaimedNodes.jMax--;
		}
		else throw UnexpectedEnumValueException<EdgeIndexEnum>(planeIndex);
		break;
	case AxisOrientationEnum::z:
		if(planeIndex == EdgeIndexEnum::min)
		{
			indexBounds.kMax = indexBounds.kMin;
			unclaimedNodes.kMin++;
		}
		else if(planeIndex == EdgeIndexEnum::max)
		{
			indexBounds.kMin = indexBounds.kMax;
			unclaimedNodes.kMax--;
		}
		else throw UnexpectedEnumValueException<EdgeIndexEnum>(planeIndex);
		break;
	default:
		throw UnexpectedEnumValueException<AxisOrientationEnum>(normalAxis);
	}
	for(uint i{indexBounds.iMin}; i<=indexBounds.iMax; ++i)
		for(uint j{indexBounds.jMin}; j<=indexBounds.jMax; ++j)
			for(uint k{indexBounds.kMin}; k<=indexBounds.kMax; ++k)
				nodeIndices.push_back(i * (meshSize.jMax+1) * (meshSize.kMax+1) + j * (meshSize.kMax+1) + k);
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

CylinderBody::CylinderBody(sf::Vector3<double> centroidPosition, AxisOrientationEnum axis, double radius) :
centroidPosition{centroidPosition},
axis{axis},
radius{radius}
{}

void CylinderBody::identifyGhostNodes(const IndexBoundingBox meshSize, double dx, double dy, double dz)
{
	uint indexRadiusX{static_cast<uint>( ceil(radius/dx) ) + 1};
	uint indexRadiusY{static_cast<uint>( ceil(radius/dy) ) + 1};
	uint indexRadiusZ{static_cast<uint>( ceil(radius/dz) ) + 1};
	uint centroidClosestIndexX{static_cast<uint>( round(centroidPosition.x/dx) )};
	uint centroidClosestIndexY{static_cast<uint>( round(centroidPosition.y/dy) )};
	uint centroidClosestIndexZ{static_cast<uint>( round(centroidPosition.z/dz) )};
	IndexBoundingBox indicesToCheck = meshSize;
	if(axis != AxisOrientationEnum::x)
	{
		indicesToCheck.iMin = centroidClosestIndexX - indexRadiusX;
		indicesToCheck.iMax = centroidClosestIndexX + indexRadiusX;
	}
	if(axis != AxisOrientationEnum::y)
	{
		indicesToCheck.jMin = centroidClosestIndexY - indexRadiusY;
		indicesToCheck.jMax = centroidClosestIndexY + indexRadiusY;
	}
	if(axis != AxisOrientationEnum::z)
	{
		indicesToCheck.kMin = centroidClosestIndexZ - indexRadiusZ;
		indicesToCheck.kMax = centroidClosestIndexZ + indexRadiusZ;
	}
	for(uint i{indicesToCheck.iMin}; i<=indicesToCheck.iMax; ++i)
		for(uint j{indicesToCheck.jMin}; j<=indicesToCheck.jMax; ++j)
			for(uint k{indicesToCheck.kMin}; k<=indicesToCheck.kMax; ++k)
			{

			}
}

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





