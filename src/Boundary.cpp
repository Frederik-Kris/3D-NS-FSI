/*
 * Boundary.cpp
 *
 *  Created on: Oct 13, 2022
 *      Author: frederk
 */

#include "Boundary.h"

Array8_u IndexBoundingBox::asIndexList(const Mesh& mesh) const
{
	Array8_u indices = { mesh.getIndex1D(iMin, jMin, kMin),
						 mesh.getIndex1D(iMin, jMin, kMax),
						 mesh.getIndex1D(iMin, jMax, kMin),
						 mesh.getIndex1D(iMin, jMax, kMax),
						 mesh.getIndex1D(iMax, jMin, kMin),
						 mesh.getIndex1D(iMax, jMin, kMax),
						 mesh.getIndex1D(iMax, jMax, kMin),
						 mesh.getIndex1D(iMax, jMax, kMax)
	};
	return indices;
}

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
	for(size_t i{indexBounds.iMin}; i<=indexBounds.iMax; ++i)
		for(size_t j{indexBounds.jMin}; j<=indexBounds.jMax; ++j)
			for(size_t k{indexBounds.kMin}; k<=indexBounds.kMax; ++k)
			{
				nodeIndices.push_back(mesh.getIndex1D(i,j,k));
				mesh.nodeTypes(i,j,k) = NodeTypeEnum::FluidEdge;
			}
}

void MeshEdgeBoundary::getAdjacentIndices(size_t index1D, const Mesh& mesh, size_t& boundaryAdjacentIndex, size_t& nextToAdjacentIndex)
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

size_t MeshEdgeBoundary::getPeriodicIndex(size_t index1D, const Mesh& mesh)
{
	Vector3_u indices = mesh.getIndices3D(index1D);
	if(normalAxis == AxisOrientationEnum::x && planeIndex == EdgeIndexEnum::min)
		indices.i = mesh.NI-2;	// i -> iMax-1

	else if(normalAxis == AxisOrientationEnum::x && planeIndex == EdgeIndexEnum::max)
		indices.i = 1;	// i -> iMin+1

	else if(normalAxis == AxisOrientationEnum::y && planeIndex == EdgeIndexEnum::min)
		indices.j = mesh.NJ-2;	// j -> jMax-1

	else if(normalAxis == AxisOrientationEnum::y && planeIndex == EdgeIndexEnum::max)
		indices.j = 1;	// j -> jMin+1

	else if(normalAxis == AxisOrientationEnum::z && planeIndex == EdgeIndexEnum::min)
		indices.k = mesh.NK-2;	// k -> kMax-1

	else if(normalAxis == AxisOrientationEnum::z && planeIndex == EdgeIndexEnum::max)
		indices.k = 1;	// k -> kMin+1

	return mesh.getIndex1D(indices);
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

	for(size_t index1D : nodeIndices)
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
		size_t boundaryAdjacentIndex{0}, nextToAdjacentIndex{0};
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
	for(size_t index1D : nodeIndices)
	{
		size_t boundaryAdjacentIndex{0}, nextToAdjacentIndex{0};
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
	for(size_t index1D : nodeIndices)
	{
		size_t oppositeSideIndex = getPeriodicIndex(index1D, mesh);
		mesh.conservedVariables.rho  (index1D) = mesh.conservedVariables.rho  (oppositeSideIndex);
		mesh.conservedVariables.rho_u(index1D) = mesh.conservedVariables.rho_u(oppositeSideIndex);
		mesh.conservedVariables.rho_v(index1D) = mesh.conservedVariables.rho_v(oppositeSideIndex);
		mesh.conservedVariables.rho_w(index1D) = mesh.conservedVariables.rho_w(oppositeSideIndex);
		mesh.conservedVariables.rho_E(index1D) = mesh.conservedVariables.rho_E(oppositeSideIndex);
		mesh.primitiveVariables.u(index1D) = mesh.primitiveVariables.u(oppositeSideIndex);
		mesh.primitiveVariables.v(index1D) = mesh.primitiveVariables.v(oppositeSideIndex);
		mesh.primitiveVariables.w(index1D) = mesh.primitiveVariables.w(oppositeSideIndex);
		mesh.primitiveVariables.p(index1D) = mesh.primitiveVariables.p(oppositeSideIndex);
		mesh.primitiveVariables.T(index1D) = mesh.primitiveVariables.T(oppositeSideIndex);
		mesh.transportProperties.mu   (index1D) = mesh.transportProperties.mu   (oppositeSideIndex);
		mesh.transportProperties.kappa(index1D) = mesh.transportProperties.kappa(oppositeSideIndex);
	}
}

SymmetryBoundary::SymmetryBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex)
: MeshEdgeBoundary(normalAxis, planeIndex)
{}

void SymmetryBoundary::applyBoundaryCondition(double t, const ConfigSettings& params, Mesh& mesh)
{
	for(size_t index1D : nodeIndices)
	{
		size_t boundaryAdjacentIndex{0}, nextToAdjacentIndex{0};
		getAdjacentIndices(index1D, mesh, boundaryAdjacentIndex, nextToAdjacentIndex);
		double uScalar, vScalar, wScalar;
		if(normalAxis == AxisOrientationEnum::x)
		{
			uScalar = 0;	// Normal component zero, and other components get zero gradient.
			vScalar = ( 4*mesh.primitiveVariables.v(boundaryAdjacentIndex) - mesh.primitiveVariables.v(nextToAdjacentIndex) ) / 3;
			wScalar = ( 4*mesh.primitiveVariables.w(boundaryAdjacentIndex) - mesh.primitiveVariables.w(nextToAdjacentIndex) ) / 3;
		}
		else if(normalAxis == AxisOrientationEnum::y)
		{
			uScalar = ( 4*mesh.primitiveVariables.u(boundaryAdjacentIndex) - mesh.primitiveVariables.u(nextToAdjacentIndex) ) / 3;
			vScalar = 0;
			wScalar = ( 4*mesh.primitiveVariables.w(boundaryAdjacentIndex) - mesh.primitiveVariables.w(nextToAdjacentIndex) ) / 3;
		}
		else if(normalAxis == AxisOrientationEnum::z)
		{
			uScalar = ( 4*mesh.primitiveVariables.u(boundaryAdjacentIndex) - mesh.primitiveVariables.u(nextToAdjacentIndex) ) / 3;
			vScalar = ( 4*mesh.primitiveVariables.v(boundaryAdjacentIndex) - mesh.primitiveVariables.v(nextToAdjacentIndex) ) / 3;
			wScalar = 0;
		}
		else
			throw std::logic_error("Unexpected enum value");
		double pScalar = ( 4*mesh.primitiveVariables.p(boundaryAdjacentIndex) - mesh.primitiveVariables.p(nextToAdjacentIndex) ) / 3;
		double TScalar = ( 4*mesh.primitiveVariables.T(boundaryAdjacentIndex) - mesh.primitiveVariables.T(nextToAdjacentIndex) ) / 3;
		PrimitiveVariablesScalars primitiveVarsLocal(uScalar, vScalar, wScalar, pScalar, TScalar);
		ConservedVariablesScalars conservedVarsLocal = deriveConservedVariables(primitiveVarsLocal, params);
		TransportPropertiesScalars transportPropsLocal = deriveTransportProperties(primitiveVarsLocal, params);
		mesh.setFlowVariablesAtNode(index1D, conservedVarsLocal, primitiveVarsLocal, transportPropsLocal);
	}
}


ImmersedBoundary::ImmersedBoundary()
{

}

void ImmersedBoundary::applyBoundaryCondition(Mesh& mesh, const ConfigSettings& params)
{
	for(GhostNode ghostNode : ghostNodes)
	{
		InterpolationValues interpolationValues;		// Primitive variables at the 8 surrounding interpolation points.
		InterpolationPositions interpolationPositions;	// Positions of interpolation points
		Array8_b ghostFlag;		// A bool flag for each surrounding node. True if ghost.
		ghostFlag.fill(false);	// <- Sets all 8 values to false
		bool allSurroundingAreFluid{true};
		vector<Vector3_d> unitNormals;	// For each of the surrounding nodes that is ghost, we put its unit normal probe here.
		IndexBoundingBox surroundingNodes = mesh.getSurroundingNodesBox(ghostNode.imagePoint);
		
		setInterpolationValues(
				surroundingNodes, mesh,								// <- Input
				interpolationValues, interpolationPositions,		// <- Output
				ghostFlag, allSurroundingAreFluid, unitNormals);	// <- Output
		
		PrimitiveVariablesScalars imagePointPrimVars = interpolatePrimitiveVariables(
				interpolationValues, interpolationPositions, allSurroundingAreFluid,
				ghostFlag, ghostNode, surroundingNodes, unitNormals, mesh);

		PrimitiveVariablesScalars ghostNodePrimVars = getGhostNodePrimitiveVariables(imagePointPrimVars);
		ConservedVariablesScalars ghostNodeConsVars = deriveConservedVariables(ghostNodePrimVars, params);
		TransportPropertiesScalars ghostNodeTransportProps = deriveTransportProperties(ghostNodePrimVars, params);
		mesh.setFlowVariablesAtNode(ghostNode.indices, ghostNodeConsVars, ghostNodePrimVars, ghostNodeTransportProps);
	}
}

void ImmersedBoundary::findGhostNodesWithFluidNeighbors(const vector<size_t>& solidNodeIndices, Mesh& mesh)
{
	IndexBoundingBox meshSize(mesh.NI-1, mesh.NJ-1, mesh.NK-1);
	for (size_t index1D : solidNodeIndices)
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
			ghostNodeMap.insert( std::pair<size_t,GhostNode*>(index1D, &ghostNodes.back()) );
			mesh.nodeTypes(index1D) = NodeTypeEnum::Ghost;
		}
	}
}

void ImmersedBoundary::checkIfSurroundingShouldBeGhost(Mesh &mesh, vector<GhostNode>& newGhostNodes, const Vector3_u &surroundingNode)
{
	if (mesh.nodeTypes(surroundingNode) == NodeTypeEnum::Solid)
	{
		newGhostNodes.emplace_back(surroundingNode);
		ghostNodes.push_back(GhostNode(surroundingNode));
		ghostNodeMap.insert( std::pair<size_t,GhostNode*>(mesh.getIndex1D(surroundingNode), &ghostNodes.back()) );
		mesh.nodeTypes(surroundingNode) = NodeTypeEnum::Ghost;
	}
}

vector<GhostNode> ImmersedBoundary::setImagePointPositions(vector<GhostNode>& ghostNodesToProcess, Mesh &mesh)
{
	vector<GhostNode> newGhostNodes;
	for (GhostNode& ghostNode : ghostNodesToProcess)
	{
		Vector3_d ghostNodePosition = mesh.getNodePosition(ghostNode.indices);
		Vector3_d normalProbe = getNormalProbe(ghostNodePosition); // from ghost to body intercept point
		ghostNode.bodyInterceptPoint = ghostNodePosition + normalProbe;
		ghostNode.imagePoint = ghostNode.bodyInterceptPoint + normalProbe;
		IndexBoundingBox surroundingNodes = mesh.getSurroundingNodesBox(ghostNode.imagePoint);
		for( size_t surroundingNodeIndex1D : surroundingNodes.asIndexList(mesh) )
			checkIfSurroundingShouldBeGhost(mesh, newGhostNodes, mesh.getIndices3D(surroundingNodeIndex1D));
	}
	return newGhostNodes;
}

void ImmersedBoundary::setInterpolationValuesFluidNode(uint counter, size_t surroundingNodeIndex1D,
													  const Mesh &mesh,
													  InterpolationValues& interpolationValues,			// OUTPUT
													  InterpolationPositions& interpolationPositions)	// OUTPUT
{
	interpolationValues.u[counter] = mesh.primitiveVariables.u(surroundingNodeIndex1D);
	interpolationValues.v[counter] = mesh.primitiveVariables.v(surroundingNodeIndex1D);
	interpolationValues.w[counter] = mesh.primitiveVariables.w(surroundingNodeIndex1D);
	interpolationValues.p[counter] = mesh.primitiveVariables.p(surroundingNodeIndex1D);
	interpolationValues.T[counter] = mesh.primitiveVariables.T(surroundingNodeIndex1D);
	Vector3_d surroundingNodePosition = mesh.getNodePosition(mesh.getIndices3D(surroundingNodeIndex1D));
	interpolationPositions.x[counter] = surroundingNodePosition.x;
	interpolationPositions.y[counter] = surroundingNodePosition.y;
	interpolationPositions.z[counter] = surroundingNodePosition.z;
}

void ImmersedBoundary::setInterpolationValuesGhostNode(
		uint counter,
		size_t surroundingNodeIndex1D,
		vector<Vector3_d>& unitNormals,					// <-
		InterpolationValues& interpolationValues,		// <- Output
		InterpolationPositions& interpolationPositions)	// <-
{
	GhostNode &surroundingGhostNode = *ghostNodeMap.at(surroundingNodeIndex1D);
	Vector3_d &unitNormal = unitNormals.emplace_back();
	unitNormal = surroundingGhostNode.imagePoint - surroundingGhostNode.bodyInterceptPoint;
	unitNormal = unitNormal / unitNormal.length();
	interpolationValues.u[counter] = 0;
	interpolationValues.v[counter] = 0;
	interpolationValues.w[counter] = 0; // NB! This hardcodes zero Dirichlet condition for velocities
	interpolationValues.p[counter] = 0; // and zero gradient Neumann conditions for pressure and temperature.
	interpolationValues.T[counter] = 0;
	interpolationPositions.x[counter] = surroundingGhostNode.bodyInterceptPoint.x;
	interpolationPositions.y[counter] = surroundingGhostNode.bodyInterceptPoint.y;
	interpolationPositions.z[counter] = surroundingGhostNode.bodyInterceptPoint.z;
}

void ImmersedBoundary::setInterpolationValues(
		const IndexBoundingBox& surroundingNodes,
		const Mesh& mesh,
		InterpolationValues& interpolationValues,		// <-
		InterpolationPositions& interpolationPositions,	// <-
		Array8_b& ghostFlag,							// <- Output
		bool& allSurroundingAreFluid,					// <-
		vector<Vector3_d>& unitNormals)					// <-
{
	uint counter { 0 }; // 0, 1, ..., 7.  To index the interpolation point arrays.
	for (size_t surroundingNodeIndex1D : surroundingNodes.asIndexList(mesh))
	{
		if (mesh.nodeTypes(surroundingNodeIndex1D) == NodeTypeEnum::FluidRegular
		||	mesh.nodeTypes(surroundingNodeIndex1D) == NodeTypeEnum::FluidEdge)
		{
			setInterpolationValuesFluidNode(counter, surroundingNodeIndex1D, mesh,	// <- Input
									interpolationValues, interpolationPositions);	// <- Output
		}
		else if (mesh.nodeTypes(surroundingNodeIndex1D) == NodeTypeEnum::Ghost)
		{
			ghostFlag[counter] = true;
			allSurroundingAreFluid = false;
			setInterpolationValuesGhostNode(counter, surroundingNodeIndex1D,	// <- Input
					unitNormals, interpolationValues, interpolationPositions);	// <- Output
		}
		else
			throw std::logic_error("Impossible situation. Found solid node around image point.");

		++counter;
	}
}

double ImmersedBoundary::simplifiedInterpolation(const Vector8_d& interpolationValues, const Vector3_u& lowerIndexNode, const Vector3_d& imagePointPosition, const Mesh& mesh)
{
	Vector3_d lowerIndexNodePosition = mesh.getNodePosition(lowerIndexNode);
	Vector3_d relativePosition = imagePointPosition - lowerIndexNodePosition;
	double variableAt01 = interpolationValues(0) + (interpolationValues(1)-interpolationValues(0)) * relativePosition.z / mesh.dz;
	double variableAt23 = interpolationValues(2) + (interpolationValues(3)-interpolationValues(2)) * relativePosition.z / mesh.dz;
	double variableAt45 = interpolationValues(4) + (interpolationValues(5)-interpolationValues(4)) * relativePosition.z / mesh.dz;
	double variableAt67 = interpolationValues(6) + (interpolationValues(7)-interpolationValues(6)) * relativePosition.z / mesh.dz;
	double variableAt0123 = variableAt01 + (variableAt23-variableAt01) * relativePosition.y / mesh.dy;
	double variableAt4567 = variableAt45 + (variableAt67-variableAt45) * relativePosition.y / mesh.dy;
	return variableAt0123 + (variableAt4567-variableAt0123) * relativePosition.x / mesh.dx;
}

PrimitiveVariablesScalars ImmersedBoundary::simplifiedInterpolationAll(const InterpolationValues& interpolationValues, const Vector3_u& lowerIndexNode, const Vector3_d& imagePointPosition, const Mesh& mesh)
{
	double uInterpolated = simplifiedInterpolation(interpolationValues.u, lowerIndexNode, imagePointPosition, mesh);
	double vInterpolated = simplifiedInterpolation(interpolationValues.v, lowerIndexNode, imagePointPosition, mesh);
	double wInterpolated = simplifiedInterpolation(interpolationValues.w, lowerIndexNode, imagePointPosition, mesh);
	double pInterpolated = simplifiedInterpolation(interpolationValues.p, lowerIndexNode, imagePointPosition, mesh);
	double TInterpolated = simplifiedInterpolation(interpolationValues.T, lowerIndexNode, imagePointPosition, mesh);
	return PrimitiveVariablesScalars(uInterpolated, vInterpolated, wInterpolated, pInterpolated, TInterpolated);
}

PrimitiveVariablesScalars ImmersedBoundary::getGhostNodePrimitiveVariables(const PrimitiveVariablesScalars& imagePointPrimVars)
{
	double uGhost = -imagePointPrimVars.u;	// NB! This hardcodes zero Dirichlet condition for velocities
	double vGhost = -imagePointPrimVars.v;
	double wGhost = -imagePointPrimVars.w;
	double pGhost =  imagePointPrimVars.p;	// and zero gradient Neumann conditions for pressure and temperature.
	double TGhost =  imagePointPrimVars.T;
	return PrimitiveVariablesScalars(uGhost, vGhost, wGhost, pGhost, TGhost);
}

void ImmersedBoundary::populateVandermondeDirichlet(const InterpolationPositions& points, Matrix8x8_d& vandermonde)
{
	vandermonde << Vector8_d::Constant(1), points.x, points.y, points.z, points.x*points.y, points.x*points.z, points.y*points.z, points.x*points.y*points.z;
}

void ImmersedBoundary::populateVandermondeNeumann(const InterpolationPositions& points,
												  const Array8_b& ghostFlags,
												  const vector<Vector3_d>& unitNormals,
												  Matrix8x8_d& vandermonde)
{
	populateVandermondeDirichlet(points, vandermonde);
	size_t counter{0};
	for(size_t rowIndex{0}; rowIndex<8; ++rowIndex)
		if(ghostFlags[rowIndex] == true)
		{
			vandermonde(rowIndex,0) = 0;
			vandermonde(rowIndex,1) = unitNormals[counter].x;
			vandermonde(rowIndex,2) = unitNormals[counter].y;
			vandermonde(rowIndex,3) = unitNormals[counter].z;
			vandermonde(rowIndex,4) = points.x(rowIndex)*unitNormals[counter].y + points.y(rowIndex)*unitNormals[counter].x;
			vandermonde(rowIndex,5) = points.x(rowIndex)*unitNormals[counter].z + points.z(rowIndex)*unitNormals[counter].x;
			vandermonde(rowIndex,6) = points.y(rowIndex)*unitNormals[counter].z + points.z(rowIndex)*unitNormals[counter].y;
			vandermonde(rowIndex,7) = points.x(rowIndex)*points.y(rowIndex)*unitNormals[counter].z
									+ points.x(rowIndex)*points.z(rowIndex)*unitNormals[counter].y
									+ points.y(rowIndex)*points.z(rowIndex)*unitNormals[counter].x;
			++counter;
		}
}

double ImmersedBoundary::trilinearInterpolation(const Vector8_d& interpolationValues,
												const Vector3_d& imagePoint,
												const Matrix8x8_d& vandermonde)
{
	Vector8_d trilinearCoefficients = vandermonde.fullPivLu().solve( interpolationValues.matrix() );
	return trilinearCoefficients(0)
		 + trilinearCoefficients(1)*imagePoint.x
		 + trilinearCoefficients(2)*imagePoint.y
		 + trilinearCoefficients(3)*imagePoint.z
		 + trilinearCoefficients(4)*imagePoint.x*imagePoint.y
		 + trilinearCoefficients(5)*imagePoint.x*imagePoint.z
		 + trilinearCoefficients(6)*imagePoint.y*imagePoint.z
		 + trilinearCoefficients(7)*imagePoint.x*imagePoint.y*imagePoint.z;
}

PrimitiveVariablesScalars ImmersedBoundary::trilinearInterpolationAll(const InterpolationValues& values,
																	  const Vector3_d& imagePoint,
																	  const Matrix8x8_d& vandermondeDirichlet,
																	  const Matrix8x8_d& vandermondeNeumann)
{
	double uInterpolated = trilinearInterpolation(values.u, imagePoint, vandermondeDirichlet);
	double vInterpolated = trilinearInterpolation(values.v, imagePoint, vandermondeDirichlet);
	double wInterpolated = trilinearInterpolation(values.w, imagePoint, vandermondeDirichlet);
	double pInterpolated = trilinearInterpolation(values.p, imagePoint, vandermondeDirichlet);
	double TInterpolated = trilinearInterpolation(values.T, imagePoint, vandermondeDirichlet);
	return PrimitiveVariablesScalars(uInterpolated, vInterpolated, wInterpolated, pInterpolated, TInterpolated);
}

PrimitiveVariablesScalars ImmersedBoundary::interpolatePrimitiveVariables(
		const InterpolationValues& interpolationValues,
		const InterpolationPositions& interpolationPositions,
		bool allSurroundingAreFluid,
		const Array8_b& ghostFlag,
		const GhostNode& ghostNode,
		const IndexBoundingBox& surroundingNodes,
		const vector<Vector3_d>& unitNormals,
		const Mesh& mesh)
{
	PrimitiveVariablesScalars imagePointPrimVars;
	if (allSurroundingAreFluid)
	{ // Then we can use the simplified interpolation method:
		Vector3_u lowerIndexNode(surroundingNodes.iMin, surroundingNodes.jMin, surroundingNodes.kMin);
		imagePointPrimVars = simplifiedInterpolationAll(interpolationValues,
														lowerIndexNode,
														ghostNode.imagePoint,
														mesh);
	}
	else
	{ // Then we must use trilinear interpolation with the Vandermonde matrix:
		Matrix8x8_d vandermondeDirichlet, vandermondeNeumann;
		populateVandermondeDirichlet(interpolationPositions, vandermondeDirichlet);
		populateVandermondeNeumann(interpolationPositions, ghostFlag, unitNormals, vandermondeNeumann);
		imagePointPrimVars = trilinearInterpolationAll(interpolationValues,
													   ghostNode.imagePoint,
													   vandermondeDirichlet,
													   vandermondeNeumann);
	}
	return imagePointPrimVars;
}

CylinderBody::CylinderBody(Vector3_d centroidPosition, AxisOrientationEnum axis, double radius) :
centroidPosition(centroidPosition),
axis{axis},
radius{radius}
{}

void CylinderBody::identifyRelatedNodes(const ConfigSettings& params, Mesh& mesh)
{
	IndexBoundingBox indicesToCheck = getCylinderBoundingBox(mesh);
	vector<size_t> solidNodeIndices;
	getSolidNodesInCylinder(params, solidNodeIndices, indicesToCheck, mesh);
	findGhostNodesWithFluidNeighbors(solidNodeIndices, mesh);

	vector<GhostNode> newGhostNodes = setImagePointPositions(ghostNodes, mesh);
	while(newGhostNodes.size() > 0)
	{
		newGhostNodes = setImagePointPositions(newGhostNodes, mesh);
	}
}

Vector3_d CylinderBody::getNormalProbe(const Vector3_d& ghostNodePosition)
{
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
	Vector3_d normalProbe = centroidToGhost * lengthFactor;
	return normalProbe;
}

IndexBoundingBox CylinderBody::getCylinderBoundingBox(Mesh& mesh) const
{
	size_t indexRadiusX { static_cast<size_t>(ceil(radius / mesh.dx)) + 1 };
	size_t indexRadiusY { static_cast<size_t>(ceil(radius / mesh.dy)) + 1 };
	size_t indexRadiusZ { static_cast<size_t>(ceil(radius / mesh.dz)) + 1 };
	size_t centroidClosestIndexX { static_cast<size_t>(round(centroidPosition.x / mesh.dx)) };
	size_t centroidClosestIndexY { static_cast<size_t>(round(centroidPosition.y / mesh.dy)) };
	size_t centroidClosestIndexZ { static_cast<size_t>(round(centroidPosition.z / mesh.dz)) };
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

void CylinderBody::getSolidNodesInCylinder(const ConfigSettings& params, vector<size_t>& solidNodeIndices, IndexBoundingBox indicesToCheck, Mesh& mesh)
{
	for (size_t i { indicesToCheck.iMin }; i <= indicesToCheck.iMax; ++i)
		for (size_t j { indicesToCheck.jMin }; j <= indicesToCheck.jMax; ++j)
			for (size_t k { indicesToCheck.kMin }; k <= indicesToCheck.kMax; ++k)
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

SphereBody::SphereBody(Vector3_d centerPosition, double radius) :
centerPosition(centerPosition),
radius{radius}
{}

void SphereBody::identifyRelatedNodes(const ConfigSettings& params, Mesh& mesh)
{
	IndexBoundingBox indicesToCheck = getSphereBoundingBox(mesh);
	vector<size_t> solidNodeIndices;
	getSolidNodesInSphere(params, solidNodeIndices, indicesToCheck, mesh);
	findGhostNodesWithFluidNeighbors(solidNodeIndices, mesh);

	vector<GhostNode> newGhostNodes = setImagePointPositions(ghostNodes, mesh);
	while(newGhostNodes.size() > 0)
	{
		newGhostNodes = setImagePointPositions(newGhostNodes, mesh);
	}
}

Vector3_d SphereBody::getNormalProbe(const Vector3_d& ghostNodePosition)
{
	Vector3_d centroidToGhost = ghostNodePosition - centerPosition;
	double lengthFactor = (radius - centroidToGhost.length()) / centroidToGhost.length();
	Vector3_d normalProbe = centroidToGhost * lengthFactor;
	return normalProbe;
}

IndexBoundingBox SphereBody::getSphereBoundingBox(Mesh& mesh) const
{
	size_t indexRadiusX { static_cast<size_t>(ceil(radius / mesh.dx)) + 1 };
	size_t indexRadiusY { static_cast<size_t>(ceil(radius / mesh.dy)) + 1 };
	size_t indexRadiusZ { static_cast<size_t>(ceil(radius / mesh.dz)) + 1 };
	size_t centerClosestIndexX { static_cast<size_t>(round(centerPosition.x / mesh.dx)) };
	size_t centerClosestIndexY { static_cast<size_t>(round(centerPosition.y / mesh.dy)) };
	size_t centerClosestIndexZ { static_cast<size_t>(round(centerPosition.z / mesh.dz)) };
	IndexBoundingBox indicesToCheck;
	indicesToCheck.iMin = centerClosestIndexX - indexRadiusX;
	indicesToCheck.iMax = centerClosestIndexX + indexRadiusX;
	indicesToCheck.jMin = centerClosestIndexY - indexRadiusY;
	indicesToCheck.jMax = centerClosestIndexY + indexRadiusY;
	indicesToCheck.kMin = centerClosestIndexZ - indexRadiusZ;
	indicesToCheck.kMax = centerClosestIndexZ + indexRadiusZ;
	return indicesToCheck;
}

void SphereBody::getSolidNodesInSphere(const ConfigSettings& params,
							 vector<size_t>& solidNodeIndices,
							 IndexBoundingBox indicesToCheck,
							 Mesh& mesh)
{
	for (size_t i { indicesToCheck.iMin }; i <= indicesToCheck.iMax; ++i)
		for (size_t j { indicesToCheck.jMin }; j <= indicesToCheck.jMax; ++j)
			for (size_t k { indicesToCheck.kMin }; k <= indicesToCheck.kMax; ++k)
			{
				Vector3_d nodePosition = mesh.getNodePosition(i, j, k);
				double distanceFromCenter = sqrt( pow(nodePosition.x - centerPosition.x, 2)
												+ pow(nodePosition.y - centerPosition.y, 2)
												+ pow(nodePosition.z - centerPosition.z, 2));
				if (distanceFromCenter < radius - params.machinePrecisionBuffer)
				{
					mesh.nodeTypes(i, j, k) = NodeTypeEnum::Solid;
					solidNodeIndices.push_back(mesh.getIndex1D(i, j, k));
				}
			}
}





