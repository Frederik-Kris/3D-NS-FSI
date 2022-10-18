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

void MeshEdgeBoundary::identifyOwnedNodes(IndexBoundingBox& unclaimedNodes, const IndexBoundingBox meshSize,
											Array3D_nodeType& nodeTypes)
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
			{
				nodeIndices.push_back(i * (meshSize.jMax+1) * (meshSize.kMax+1) + j * (meshSize.kMax+1) + k);
				nodeTypes(i,j,k) = NodeTypeEnum::FluidEdge;
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

CylinderBody::CylinderBody(sf::Vector3<double> centroidPosition, AxisOrientationEnum axis, double radius) :
centroidPosition{centroidPosition},
axis{axis},
radius{radius}
{}

IndexBoundingBox CylinderBody::getCylinderBoundingBox(const IndexBoundingBox& meshSize, double dx, double dy, double dz)
{
	uint indexRadiusX { static_cast<uint>(ceil(radius / dx)) + 1 };
	uint indexRadiusY { static_cast<uint>(ceil(radius / dy)) + 1 };
	uint indexRadiusZ { static_cast<uint>(ceil(radius / dz)) + 1 };
	uint centroidClosestIndexX { static_cast<uint>(round(centroidPosition.x / dx)) };
	uint centroidClosestIndexY { static_cast<uint>(round(centroidPosition.y / dy)) };
	uint centroidClosestIndexZ { static_cast<uint>(round(centroidPosition.z / dz)) };
	IndexBoundingBox indicesToCheck = meshSize;
	if (axis != AxisOrientationEnum::x)
	{
		indicesToCheck.iMin = centroidClosestIndexX - indexRadiusX;
		indicesToCheck.iMax = centroidClosestIndexX + indexRadiusX;
	}
	if (axis != AxisOrientationEnum::y)
	{
		indicesToCheck.jMin = centroidClosestIndexY - indexRadiusY;
		indicesToCheck.jMax = centroidClosestIndexY + indexRadiusY;
	}
	if (axis != AxisOrientationEnum::z)
	{
		indicesToCheck.kMin = centroidClosestIndexZ - indexRadiusZ;
		indicesToCheck.kMax = centroidClosestIndexZ + indexRadiusZ;
	}
	return indicesToCheck;
}

vector<uint>& CylinderBody::getSolidNodesInCylinder(const ConfigSettings& params,
		IndexBoundingBox indicesToCheck, double dx, double dy, double dz,
		const IndexBoundingBox& meshSize, Array3D_nodeType &nodeTypes)
{
	vector<uint> solidNodeIndices;
	for (uint i { indicesToCheck.iMin }; i <= indicesToCheck.iMax; ++i)
		for (uint j { indicesToCheck.jMin }; j <= indicesToCheck.jMax; ++j)
			for (uint k { indicesToCheck.kMin }; k <= indicesToCheck.kMax; ++k)
			{
				double x { i * dx }, y { j * dy }, z { k * dz };
				if (axis == AxisOrientationEnum::x)
				{
					double distanceFromCentroid { sqrt(
							pow(y - centroidPosition.y, 2)
									+ pow(z - centroidPosition.z, 2)) };
					if (distanceFromCentroid < radius - params.machinePrecisionBuffer)
					{
						nodeTypes(i, j, k) = NodeTypeEnum::Solid;
						solidNodeIndices.push_back( i * (meshSize.jMax + 1) * (meshSize.kMax + 1)
												  + j * (meshSize.kMax + 1)
												  + k);
					}
				}
			}
	return solidNodeIndices;
}

void CylinderBody::findGhostNodes(const vector<uint>& solidNodeIndices,
		const IndexBoundingBox& meshSize, const Array3D_nodeType& nodeTypes)
{
	for (uint index1D : solidNodeIndices)
	{
		sf::Vector3i indices = Mesh::getIndices3D(index1D);
		bool solidNodeHasFluidNeighbor { false };
		if (indices.x > meshSize.iMin)
			if (nodeTypes(indices.x - 1, indices.y, indices.z) == NodeTypeEnum::FluidRegular)
				solidNodeHasFluidNeighbor = true;

		if (indices.x < meshSize.iMax)
			if (nodeTypes(indices.x + 1, indices.y, indices.z) == NodeTypeEnum::FluidRegular)
				solidNodeHasFluidNeighbor = true;

		if (indices.y > meshSize.jMin)
			if (nodeTypes(indices.x, indices.y - 1, indices.z) == NodeTypeEnum::FluidRegular)
				solidNodeHasFluidNeighbor = true;

		if (indices.y < meshSize.jMax)
			if (nodeTypes(indices.x, indices.y + 1, indices.z) == NodeTypeEnum::FluidRegular)
				solidNodeHasFluidNeighbor = true;

		if (indices.z > meshSize.kMin)
			if (nodeTypes(indices.x, indices.y, indices.z - 1) == NodeTypeEnum::FluidRegular)
				solidNodeHasFluidNeighbor = true;

		if (indices.z < meshSize.kMax)
			if (nodeTypes(indices.x, indices.y, indices.z + 1) == NodeTypeEnum::FluidRegular)
				solidNodeHasFluidNeighbor = true;

		if (solidNodeHasFluidNeighbor)
			ghostNodeIndices.push_back(index1D);
	}
}

void CylinderBody::identifyGhostNodes(const ConfigSettings& params,
									  const IndexBoundingBox meshSize,
									  Array3D_nodeType& nodeTypes,
									  double dx, double dy, double dz)
{
	IndexBoundingBox indicesToCheck = getCylinderBoundingBox(meshSize, dx, dy, dz);
	vector<uint> solidNodeIndices;
	solidNodeIndices = getSolidNodesInCylinder(params, indicesToCheck, dx, dy, dz, meshSize, nodeTypes);
	findGhostNodes(solidNodeIndices, meshSize, nodeTypes);

	// IKKE GLEM DENNE GANGEN Å SJEKKE REKURSIVT SÅNN AT IKKE VI MANGLER GHOSTS SOM LIGGER RUNDT IMAGE POINT!
	// LAG EN FUNC SOM TAR INN VEKTOR MED GHOSTS, OG FINNER IMAGE POINT GREIER, OG RETURNERER NY-OPPSTÅTTE
	// GHOSTS. DA KAN JEG KALLE DEN I WHILE, HELT TIL DET IKKE BLIR FLERE.
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





