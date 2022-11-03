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

void MeshEdgeBoundary::identifyOwnedNodes(IndexBoundingBox& unclaimedNodes, Mesh& mesh)
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
		else throw std::logic_error("Unexpected enum value");
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
		else throw std::logic_error("Unexpected enum value");
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
		else throw std::logic_error("Unexpected enum value");
		break;
	default:
		throw std::logic_error("Unexpected enum value");
	}
	for(uint i{indexBounds.iMin}; i<=indexBounds.iMax; ++i)
		for(uint j{indexBounds.jMin}; j<=indexBounds.jMax; ++j)
			for(uint k{indexBounds.kMin}; k<=indexBounds.kMax; ++k)
			{
				nodeIndices.push_back(mesh.getIndex1D(i,j,k));
				mesh.nodeTypes(i,j,k) = NodeTypeEnum::FluidEdge;
			}
}

void MeshEdgeBoundary::getAdjacentIndices(uint index1D, const Mesh& mesh, uint& boundaryAdjacentIndex, uint& nextToAdjacentIndex)
{
	if(normalAxis == AxisOrientationEnum::x && planeIndex == EdgeIndexEnum::min)
	{
		boundaryAdjacentIndex = index1D + mesh.NJ*mesh.NK;
		nextToAdjacentIndex = index1D + 2*mesh.NJ*mesh.NK;
	}
	else if(normalAxis == AxisOrientationEnum::x && planeIndex == EdgeIndexEnum::max)
	{
		boundaryAdjacentIndex = index1D - mesh.NJ*mesh.NK;
		nextToAdjacentIndex = index1D - 2*mesh.NJ*mesh.NK;
	}
	else if(normalAxis == AxisOrientationEnum::y && planeIndex == EdgeIndexEnum::min)
	{
		boundaryAdjacentIndex = index1D + mesh.NK;
		nextToAdjacentIndex = index1D + 2*mesh.NK;
	}
	else if(normalAxis == AxisOrientationEnum::y && planeIndex == EdgeIndexEnum::max)
	{
		boundaryAdjacentIndex = index1D - mesh.NK;
		nextToAdjacentIndex = index1D - 2*mesh.NK;
	}
	else if(normalAxis == AxisOrientationEnum::z && planeIndex == EdgeIndexEnum::min)
	{
		boundaryAdjacentIndex = index1D + 1;
		nextToAdjacentIndex = index1D + 2*1;
	}
	else if(normalAxis == AxisOrientationEnum::z && planeIndex == EdgeIndexEnum::max)
	{
		boundaryAdjacentIndex = index1D - 1;
		nextToAdjacentIndex = index1D - 2*1;
	}
}

InletBoundary::InletBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex, double velocity)
: MeshEdgeBoundary(normalAxis, planeIndex),
  velocity{velocity}
{}

void InletBoundary::applyBoundaryCondition(double t, const ConfigSettings& params, Mesh& mesh)
{
	double inletVelocity = min(1., t/10.) * velocity; // TODO: move magic const 10 to params
	if(planeIndex == EdgeIndexEnum::max)	// If we're on the highest index, velocity must be
		inletVelocity *= -1;				// negative, for this to be an inlet.

	for(uint index1D : nodeIndices)
	{
		double uScalar, vScalar, wScalar;
		if(normalAxis == AxisOrientationEnum::x)
		{
			uScalar = inletVelocity;
			vScalar = 0;
			wScalar = 0;
		}
		else if(normalAxis == AxisOrientationEnum::y)
		{
			uScalar = 0;
			vScalar = inletVelocity;
			wScalar = 0;
		}
		else if(normalAxis == AxisOrientationEnum::z)
		{
			uScalar = 0;
			vScalar = 0;
			wScalar = inletVelocity;
		}
		else
			throw std::logic_error("Unexpected enum value");
		uint boundaryAdjacentIndex{0}, nextToAdjacentIndex{0};
		getAdjacentIndices(index1D, mesh, boundaryAdjacentIndex, nextToAdjacentIndex);
		double pScalar = 2*mesh.primitiveVariables.p(boundaryAdjacentIndex) // Linear extrapolation
						 - mesh.primitiveVariables.p(nextToAdjacentIndex);
		double TScalar = 0;
		PrimitiveVariablesScalars primitiveVarsLocal(uScalar, vScalar, wScalar, pScalar, TScalar);
		ConservedVariablesScalars   conservedVarsLocal = deriveConservedVariables (primitiveVarsLocal, params);
		TransportPropertiesScalars transportPropsLocal = deriveTransportProperties(primitiveVarsLocal, params);
		mesh.setFlowVariablesAtNode(index1D, conservedVarsLocal, primitiveVarsLocal, transportPropsLocal);
	}
}

OutletBoundary::OutletBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex)
: MeshEdgeBoundary(normalAxis, planeIndex)
{}

void OutletBoundary::applyBoundaryCondition(double t, const ConfigSettings& params, Mesh& mesh)
{
	for(uint index1D : nodeIndices)
	{
		uint boundaryAdjacentIndex{0}, nextToAdjacentIndex{0};
		getAdjacentIndices(index1D, mesh, boundaryAdjacentIndex, nextToAdjacentIndex);
		double uScalar = 2*mesh.primitiveVariables.u(boundaryAdjacentIndex) // Linear extrapolation
						 - mesh.primitiveVariables.u(nextToAdjacentIndex);
		double vScalar = 2*mesh.primitiveVariables.v(boundaryAdjacentIndex) // Linear extrapolation
						 - mesh.primitiveVariables.v(nextToAdjacentIndex);
		double wScalar = 2*mesh.primitiveVariables.w(boundaryAdjacentIndex) // Linear extrapolation
						 - mesh.primitiveVariables.w(nextToAdjacentIndex);
		double pScalar = 0;
		double TScalar = 2*mesh.primitiveVariables.T(boundaryAdjacentIndex) // Linear extrapolation
						 - mesh.primitiveVariables.T(nextToAdjacentIndex);

		PrimitiveVariablesScalars primitiveVarsLocal(uScalar, vScalar, wScalar, pScalar, TScalar);
		ConservedVariablesScalars conservedVarsLocal = deriveConservedVariables(primitiveVarsLocal, params);
		TransportPropertiesScalars transportPropsLocal = deriveTransportProperties(primitiveVarsLocal, params);
		mesh.setFlowVariablesAtNode(index1D, conservedVarsLocal, primitiveVarsLocal, transportPropsLocal);
	}
}

PeriodicBoundary::PeriodicBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex)
: MeshEdgeBoundary(normalAxis, planeIndex)
{}

void PeriodicBoundary::applyBoundaryCondition(double t, const ConfigSettings& params, Mesh& mesh)
{

}

SymmetryBoundary::SymmetryBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex)
: MeshEdgeBoundary(normalAxis, planeIndex)
{}

void SymmetryBoundary::applyBoundaryCondition(double t, const ConfigSettings& params, Mesh& mesh)
{

}


ImmersedBoundary::ImmersedBoundary()
{

}

CylinderBody::CylinderBody(Vector3_d centroidPosition, AxisOrientationEnum axis, double radius) :
centroidPosition(centroidPosition),
axis{axis},
radius{radius}
{}

IndexBoundingBox CylinderBody::getCylinderBoundingBox(Mesh& mesh) const
{
	uint indexRadiusX { static_cast<uint>(ceil(radius / mesh.dx)) + 1 };
	uint indexRadiusY { static_cast<uint>(ceil(radius / mesh.dy)) + 1 };
	uint indexRadiusZ { static_cast<uint>(ceil(radius / mesh.dz)) + 1 };
	uint centroidClosestIndexX { static_cast<uint>(round(centroidPosition.x / mesh.dx)) };
	uint centroidClosestIndexY { static_cast<uint>(round(centroidPosition.y / mesh.dy)) };
	uint centroidClosestIndexZ { static_cast<uint>(round(centroidPosition.z / mesh.dz)) };
	IndexBoundingBox indicesToCheck(mesh.NI-1, mesh.NJ-1, mesh.NK-1);
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

void CylinderBody::getSolidNodesInCylinder(const ConfigSettings& params, vector<uint>& solidNodeIndices, IndexBoundingBox indicesToCheck, Mesh& mesh)
{
	for (uint i { indicesToCheck.iMin }; i <= indicesToCheck.iMax; ++i)
		for (uint j { indicesToCheck.jMin }; j <= indicesToCheck.jMax; ++j)
			for (uint k { indicesToCheck.kMin }; k <= indicesToCheck.kMax; ++k)
			{
				double distanceFromCentroid{0};
				Vector3_d nodePosition = mesh.getNodePosition(i, j, k);
				if (axis == AxisOrientationEnum::x)
					distanceFromCentroid = sqrt( pow(nodePosition.y - centroidPosition.y, 2)
											   + pow(nodePosition.z - centroidPosition.z, 2));
				else if (axis == AxisOrientationEnum::y)
					distanceFromCentroid = sqrt( pow(nodePosition.x - centroidPosition.x, 2)
											   + pow(nodePosition.z - centroidPosition.z, 2));
				else if (axis == AxisOrientationEnum::z)
					distanceFromCentroid = sqrt( pow(nodePosition.x - centroidPosition.x, 2)
											   + pow(nodePosition.y - centroidPosition.y, 2));
				else
					throw std::logic_error("Unexpected enum value");
				if (distanceFromCentroid < radius - params.machinePrecisionBuffer)
				{
					mesh.nodeTypes(i, j, k) = NodeTypeEnum::Solid;
					solidNodeIndices.push_back(mesh.getIndex1D(i, j, k));
				}
			}
}

void CylinderBody::findGhostNodesWithFluidNeighbors(const vector<uint>& solidNodeIndices, Mesh& mesh)
{
	IndexBoundingBox meshSize(mesh.NI-1, mesh.NJ-1, mesh.NK-1);
	for (uint index1D : solidNodeIndices)
	{
		Vector3_u indices = mesh.getIndices3D(index1D);
		bool solidNodeHasFluidNeighbor { false };
		if (indices.i > meshSize.iMin)
			if (mesh.nodeTypes(indices.i - 1, indices.j, indices.k) == NodeTypeEnum::FluidRegular)
				solidNodeHasFluidNeighbor = true;

		if (indices.i < meshSize.iMax)
			if (mesh.nodeTypes(indices.i + 1, indices.j, indices.k) == NodeTypeEnum::FluidRegular)
				solidNodeHasFluidNeighbor = true;

		if (indices.j > meshSize.jMin)
			if (mesh.nodeTypes(indices.i, indices.j - 1, indices.k) == NodeTypeEnum::FluidRegular)
				solidNodeHasFluidNeighbor = true;

		if (indices.j < meshSize.jMax)
			if (mesh.nodeTypes(indices.i, indices.j + 1, indices.k) == NodeTypeEnum::FluidRegular)
				solidNodeHasFluidNeighbor = true;

		if (indices.k > meshSize.kMin)
			if (mesh.nodeTypes(indices.i, indices.j, indices.k - 1) == NodeTypeEnum::FluidRegular)
				solidNodeHasFluidNeighbor = true;

		if (indices.k < meshSize.kMax)
			if (mesh.nodeTypes(indices.i, indices.j, indices.k + 1) == NodeTypeEnum::FluidRegular)
				solidNodeHasFluidNeighbor = true;

		if (solidNodeHasFluidNeighbor)
		{
			ghostNodes.push_back(GhostNode(indices));
			mesh.nodeTypes(index1D) = NodeTypeEnum::Ghost;
		}
	}
}

void CylinderBody::checkIfSurroundingShouldBeGhost(Mesh &mesh, vector<GhostNode>& newGhostNodes, const Vector3_u &surroundingNode)
{
	if (mesh.nodeTypes(surroundingNode) == NodeTypeEnum::Solid) {
		newGhostNodes.emplace_back(surroundingNode);
		mesh.nodeTypes(surroundingNode) = NodeTypeEnum::Ghost;
	}
}

vector<GhostNode> CylinderBody::setImagePointPositions(vector<GhostNode>& ghostNodesToProcess, Mesh &mesh)
{
	vector<GhostNode> newGhostNodes;
	for (GhostNode& ghostNode : ghostNodesToProcess)
	{
		Vector3_d ghostNodePosition = mesh.getNodePosition(ghostNode.indices);
		Vector3_d centroidToGhost = ghostNodePosition - centroidPosition;
		if (axis == AxisOrientationEnum::x)
			centroidToGhost.x = 0;
		else if (axis == AxisOrientationEnum::y)
			centroidToGhost.y = 0;
		else if (axis == AxisOrientationEnum::z)
			centroidToGhost.z = 0;
		else
			throw std::logic_error("Unexpected enum value");

		double lengthFactor = (radius - centroidToGhost.length()) / centroidToGhost.length();
		Vector3_d normalProbe = centroidToGhost * lengthFactor; // from ghost to body intercept point
		ghostNode.bodyInterceptPoint = ghostNodePosition + normalProbe;
		ghostNode.imagePoint = ghostNode.bodyInterceptPoint + normalProbe;
		IndexBoundingBox surroundingNodes = mesh.getSurroundingNodesBox(ghostNode.imagePoint);
		checkIfSurroundingShouldBeGhost(mesh, newGhostNodes,
				Vector3_u(surroundingNodes.iMin, surroundingNodes.jMin, surroundingNodes.kMin));
		checkIfSurroundingShouldBeGhost(mesh, newGhostNodes,
				Vector3_u(surroundingNodes.iMin, surroundingNodes.jMin, surroundingNodes.kMax));
		checkIfSurroundingShouldBeGhost(mesh, newGhostNodes,
				Vector3_u(surroundingNodes.iMin, surroundingNodes.jMax, surroundingNodes.kMin));
		checkIfSurroundingShouldBeGhost(mesh, newGhostNodes,
				Vector3_u(surroundingNodes.iMin, surroundingNodes.jMax, surroundingNodes.kMax));
		checkIfSurroundingShouldBeGhost(mesh, newGhostNodes,
				Vector3_u(surroundingNodes.iMax, surroundingNodes.jMin, surroundingNodes.kMin));
		checkIfSurroundingShouldBeGhost(mesh, newGhostNodes,
				Vector3_u(surroundingNodes.iMax, surroundingNodes.jMin, surroundingNodes.kMax));
		checkIfSurroundingShouldBeGhost(mesh, newGhostNodes,
				Vector3_u(surroundingNodes.iMax, surroundingNodes.jMax, surroundingNodes.kMin));
		checkIfSurroundingShouldBeGhost(mesh, newGhostNodes,
				Vector3_u(surroundingNodes.iMax, surroundingNodes.jMax, surroundingNodes.kMax));
	}
	return newGhostNodes;
}

void CylinderBody::identifyRelatedNodes(const ConfigSettings& params, Mesh& mesh)
{
	IndexBoundingBox indicesToCheck = getCylinderBoundingBox(mesh);
	vector<uint> solidNodeIndices;
	getSolidNodesInCylinder(params, solidNodeIndices, indicesToCheck, mesh);
	findGhostNodesWithFluidNeighbors(solidNodeIndices, mesh);

	vector<GhostNode> newGhostNodes = setImagePointPositions(ghostNodes, mesh);
	while(newGhostNodes.size() > 0)
	{
		newGhostNodes = setImagePointPositions(newGhostNodes, mesh);
	}
}

void CylinderBody::applyBoundaryCondition(Mesh& mesh)
{

}

SphereBody::SphereBody(Vector3_d centerPosition, double radius) :
centerPosition(centerPosition),
radius{radius}
{}

void SphereBody::identifyRelatedNodes(const ConfigSettings& params, Mesh& mesh)
{

}

void SphereBody::applyBoundaryCondition(Mesh& mesh)
{

}





