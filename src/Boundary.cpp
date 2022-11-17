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

void MeshEdgeBoundary::identifyOwnedNodes(	IndexBoundingBox& unclaimedNodes,
											const Vector3_u& nMeshNodes,
											Array3D_nodeType& nodeTypeArray )
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
				nodeIndices.push_back(getIndex1D(i,j,k, nMeshNodes));
				nodeTypeArray(i,j,k) = NodeTypeEnum::FluidEdge;
			}
}

void MeshEdgeBoundary::getAdjacentIndices(size_t index1D, const Vector3_u& nMeshNodes,	// <- Input
										  size_t& boundaryAdjacentIndex, 				// <- Output
										  size_t& nextToAdjacentIndex)					// <- Output
{
	if(normalAxis == AxisOrientationEnum::x && planeIndex == EdgeIndexEnum::min)
	{
		boundaryAdjacentIndex = index1D + nMeshNodes.j * nMeshNodes.k;
		nextToAdjacentIndex   = index1D + nMeshNodes.j * nMeshNodes.k * 2;
	}
	else if(normalAxis == AxisOrientationEnum::x && planeIndex == EdgeIndexEnum::max)
	{
		boundaryAdjacentIndex = index1D - nMeshNodes.j * nMeshNodes.k;
		nextToAdjacentIndex   = index1D - nMeshNodes.j * nMeshNodes.k * 2;
	}
	else if(normalAxis == AxisOrientationEnum::y && planeIndex == EdgeIndexEnum::min)
	{
		boundaryAdjacentIndex = index1D + nMeshNodes.k;
		nextToAdjacentIndex   = index1D + nMeshNodes.k * 2;
	}
	else if(normalAxis == AxisOrientationEnum::y && planeIndex == EdgeIndexEnum::max)
	{
		boundaryAdjacentIndex = index1D - nMeshNodes.k;
		nextToAdjacentIndex   = index1D - nMeshNodes.k * 2;
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

size_t MeshEdgeBoundary::getPeriodicIndex(size_t index1D, const Vector3_u& nMeshNodes)
{
	Vector3_u indices = getIndices3D(index1D, nMeshNodes);
	if(normalAxis == AxisOrientationEnum::x && planeIndex == EdgeIndexEnum::min)
		indices.i = nMeshNodes.i - 2;	// i -> iMax-1

	else if(normalAxis == AxisOrientationEnum::x && planeIndex == EdgeIndexEnum::max)
		indices.i = 1;	// i -> iMin+1

	else if(normalAxis == AxisOrientationEnum::y && planeIndex == EdgeIndexEnum::min)
		indices.j = nMeshNodes.j - 2;	// j -> jMax-1

	else if(normalAxis == AxisOrientationEnum::y && planeIndex == EdgeIndexEnum::max)
		indices.j = 1;	// j -> jMin+1

	else if(normalAxis == AxisOrientationEnum::z && planeIndex == EdgeIndexEnum::min)
		indices.k = nMeshNodes.k - 2;	// k -> kMax-1

	else if(normalAxis == AxisOrientationEnum::z && planeIndex == EdgeIndexEnum::max)
		indices.k = 1;	// k -> kMin+1

	return getIndex1D(indices, nMeshNodes);
}

InletBoundary::InletBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex, double velocity)
: MeshEdgeBoundary(normalAxis, planeIndex),
  velocity{velocity}
{}

void InletBoundary::applyBoundaryCondition(double t, const Vector3_u& nMeshNodes,		// <- Input
										   const ConfigSettings& params,				// <- Input
										   AllFlowVariablesArrayGroup& flowVariables)	// <- Output
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
		getAdjacentIndices(index1D, nMeshNodes, boundaryAdjacentIndex, nextToAdjacentIndex);
		double pScalar = 2*flowVariables.primitiveVariables.p(boundaryAdjacentIndex) // Linear extrapolation
						 - flowVariables.primitiveVariables.p(nextToAdjacentIndex);
		double TScalar = 0;
		PrimitiveVariablesScalars primitiveVarsLocal(uScalar, vScalar, wScalar, pScalar, TScalar);
		ConservedVariablesScalars   conservedVarsLocal = deriveConservedVariables (primitiveVarsLocal, params);
		TransportPropertiesScalars transportPropsLocal = deriveTransportProperties(primitiveVarsLocal, params);
		setFlowVariablesAtNode(index1D, conservedVarsLocal, primitiveVarsLocal, transportPropsLocal, flowVariables);
	}
}

OutletBoundary::OutletBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex)
: MeshEdgeBoundary(normalAxis, planeIndex)
{}

void OutletBoundary::applyBoundaryCondition(double t, const Vector3_u& nMeshNodes,		// <- Input
		   	   	   	   	   	   	   	   	   const ConfigSettings& params,				// <- Input
										   AllFlowVariablesArrayGroup& flowVariables)	// <- Output
{
	for(size_t index1D : nodeIndices)
	{
		size_t boundaryAdjacentIndex{0}, nextToAdjacentIndex{0};
		getAdjacentIndices(index1D, nMeshNodes, boundaryAdjacentIndex, nextToAdjacentIndex);
		double rho_Scalar = 2*flowVariables.conservedVariables.rho(boundaryAdjacentIndex) // Linear extrapolation
							- flowVariables.conservedVariables.rho(nextToAdjacentIndex);
		double rho_u_Scalar = 2*flowVariables.conservedVariables.rho_u(boundaryAdjacentIndex) // Linear extrapolation
							  - flowVariables.conservedVariables.rho_u(nextToAdjacentIndex);
		double rho_v_Scalar = 2*flowVariables.conservedVariables.rho_v(boundaryAdjacentIndex) // Linear extrapolation
							  - flowVariables.conservedVariables.rho_v(nextToAdjacentIndex);
		double rho_w_Scalar = 2*flowVariables.conservedVariables.rho_w(boundaryAdjacentIndex) // Linear extrapolation
							  - flowVariables.conservedVariables.rho_w(nextToAdjacentIndex);
		double pScalar = 0;
		double uScalar = rho_u_Scalar / (1+rho_Scalar);
		double vScalar = rho_v_Scalar / (1+rho_Scalar);
		double wScalar = rho_w_Scalar / (1+rho_Scalar);
		double TScalar = ( params.Gamma * pScalar - rho_Scalar ) / ( 1 + rho_Scalar );
		double rho_E_Scalar = pScalar / ( params.Gamma - 1 )
				+ (1 + rho_Scalar)/2 * ( uScalar*uScalar + vScalar*vScalar + wScalar*wScalar );

		PrimitiveVariablesScalars primitiveVarsLocal(uScalar, vScalar, wScalar, pScalar, TScalar);
		ConservedVariablesScalars conservedVarsLocal(rho_Scalar, rho_u_Scalar, rho_v_Scalar, rho_w_Scalar, rho_E_Scalar);
		TransportPropertiesScalars transportPropsLocal = deriveTransportProperties(primitiveVarsLocal, params);
		setFlowVariablesAtNode(index1D, conservedVarsLocal, primitiveVarsLocal, transportPropsLocal, flowVariables);
	}
}

PeriodicBoundary::PeriodicBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex)
: MeshEdgeBoundary(normalAxis, planeIndex)
{}

void PeriodicBoundary::applyBoundaryCondition(double t, const Vector3_u& nMeshNodes,		// <- Input
	   	   	   	   	   	   	   	   	   	   	  const ConfigSettings& params,					// <- Input
											  AllFlowVariablesArrayGroup& flowVariables)	// <- Output
{
	for(size_t index1D : nodeIndices)
	{
		size_t oppositeSideIndex = getPeriodicIndex(index1D, nMeshNodes);
		flowVariables.conservedVariables.rho  (index1D) = flowVariables.conservedVariables.rho  (oppositeSideIndex);
		flowVariables.conservedVariables.rho_u(index1D) = flowVariables.conservedVariables.rho_u(oppositeSideIndex);
		flowVariables.conservedVariables.rho_v(index1D) = flowVariables.conservedVariables.rho_v(oppositeSideIndex);
		flowVariables.conservedVariables.rho_w(index1D) = flowVariables.conservedVariables.rho_w(oppositeSideIndex);
		flowVariables.conservedVariables.rho_E(index1D) = flowVariables.conservedVariables.rho_E(oppositeSideIndex);
		flowVariables.primitiveVariables.u(index1D) = flowVariables.primitiveVariables.u(oppositeSideIndex);
		flowVariables.primitiveVariables.v(index1D) = flowVariables.primitiveVariables.v(oppositeSideIndex);
		flowVariables.primitiveVariables.w(index1D) = flowVariables.primitiveVariables.w(oppositeSideIndex);
		flowVariables.primitiveVariables.p(index1D) = flowVariables.primitiveVariables.p(oppositeSideIndex);
		flowVariables.primitiveVariables.T(index1D) = flowVariables.primitiveVariables.T(oppositeSideIndex);
		flowVariables.transportProperties.mu   (index1D) = flowVariables.transportProperties.mu   (oppositeSideIndex);
		flowVariables.transportProperties.kappa(index1D) = flowVariables.transportProperties.kappa(oppositeSideIndex);
	}
}

SymmetryBoundary::SymmetryBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex)
: MeshEdgeBoundary(normalAxis, planeIndex)
{}

void SymmetryBoundary::applyBoundaryCondition(double t, const Vector3_u& nMeshNodes,		// <- Input
	   	   	   	   	   	  	  	  	  	  	  const ConfigSettings& params,					// <- Input
											  AllFlowVariablesArrayGroup& flowVariables)	// <- Output
{
	for(size_t index1D : nodeIndices)
	{
		size_t boundaryAdjacentIndex{0}, nextToAdjacentIndex{0};
		getAdjacentIndices(index1D, nMeshNodes, boundaryAdjacentIndex, nextToAdjacentIndex);
		double uScalar, vScalar, wScalar;
		if(normalAxis == AxisOrientationEnum::x)
		{
			uScalar = 0;	// Normal component zero, and other components get zero gradient.
			vScalar = ( 4*flowVariables.primitiveVariables.v(boundaryAdjacentIndex)
						- flowVariables.primitiveVariables.v(nextToAdjacentIndex) ) / 3;
			wScalar = ( 4*flowVariables.primitiveVariables.w(boundaryAdjacentIndex)
						- flowVariables.primitiveVariables.w(nextToAdjacentIndex) ) / 3;
		}
		else if(normalAxis == AxisOrientationEnum::y)
		{
			uScalar = ( 4*flowVariables.primitiveVariables.u(boundaryAdjacentIndex)
						- flowVariables.primitiveVariables.u(nextToAdjacentIndex) ) / 3;
			vScalar = 0;
			wScalar = ( 4*flowVariables.primitiveVariables.w(boundaryAdjacentIndex)
						- flowVariables.primitiveVariables.w(nextToAdjacentIndex) ) / 3;
		}
		else if(normalAxis == AxisOrientationEnum::z)
		{
			uScalar = ( 4*flowVariables.primitiveVariables.u(boundaryAdjacentIndex)
						- flowVariables.primitiveVariables.u(nextToAdjacentIndex) ) / 3;
			vScalar = ( 4*flowVariables.primitiveVariables.v(boundaryAdjacentIndex)
						- flowVariables.primitiveVariables.v(nextToAdjacentIndex) ) / 3;
			wScalar = 0;
		}
		else
			throw std::logic_error("Unexpected enum value");
		double pScalar = ( 4*flowVariables.primitiveVariables.p(boundaryAdjacentIndex)
						   - flowVariables.primitiveVariables.p(nextToAdjacentIndex) ) / 3;
		double TScalar = ( 4*flowVariables.primitiveVariables.T(boundaryAdjacentIndex)
						   - flowVariables.primitiveVariables.T(nextToAdjacentIndex) ) / 3;
		PrimitiveVariablesScalars primitiveVarsLocal(uScalar, vScalar, wScalar, pScalar, TScalar);
		ConservedVariablesScalars conservedVarsLocal = deriveConservedVariables(primitiveVarsLocal, params);
		TransportPropertiesScalars transportPropsLocal = deriveTransportProperties(primitiveVarsLocal, params);
		setFlowVariablesAtNode(index1D, conservedVarsLocal, primitiveVarsLocal, transportPropsLocal, flowVariables);
	}
}


ImmersedBoundary::ImmersedBoundary()
{

}

void ImmersedBoundary::applyBoundaryCondition(const Vector3_u& nMeshNodes,
											  const Vector3_d& gridSpacing,
											  const ConfigSettings& params,
											  const Array3D_nodeType& nodeTypeArray,
											  AllFlowVariablesArrayGroup& flowVariables // <- Output
											  )
{
	for(GhostNode ghostNode : ghostNodes)
	{
		InterpolationValues interpolationValues;		// Primitive variables at the 8 surrounding interpolation points.
		InterpolationPositions interpolationPositions;	// Positions of interpolation points
		Array8_b ghostFlag;		// A bool flag for each surrounding node. True if ghost.
		ghostFlag.fill(false);	// <- Sets all 8 values to false
		bool allSurroundingAreFluid{true};
		vector<Vector3_d> unitNormals;	// For each of the surrounding nodes that is ghost, we put its unit normal probe here.
		IndexBoundingBox surroundingNodes = getSurroundingNodesBox(ghostNode.imagePoint, gridSpacing);
		
		setInterpolationValues(
				surroundingNodes, nMeshNodes, gridSpacing,			// <- Input
				nodeTypeArray, flowVariables,						// <- Input
				interpolationValues, interpolationPositions,		// <- Output
				ghostFlag, allSurroundingAreFluid, unitNormals);	// <- Output
		
		PrimitiveVariablesScalars imagePointPrimVars = interpolatePrimitiveVariables(
				interpolationValues, interpolationPositions, allSurroundingAreFluid,
				ghostFlag, ghostNode, surroundingNodes, unitNormals, gridSpacing);

		PrimitiveVariablesScalars ghostNodePrimVars = getGhostNodePrimitiveVariables(imagePointPrimVars);
		ConservedVariablesScalars ghostNodeConsVars = deriveConservedVariables(ghostNodePrimVars, params);
		TransportPropertiesScalars ghostNodeTransportProps = deriveTransportProperties(ghostNodePrimVars, params);
		setFlowVariablesAtNode(ghostNode.indices, ghostNodeConsVars, ghostNodePrimVars, ghostNodeTransportProps,
								flowVariables);
	}
}

void ImmersedBoundary::findGhostNodesWithFluidNeighbors(const vector<size_t>& solidNodeIndices,
														const Vector3_u& nMeshNodes,
														Array3D_nodeType& nodeTypeArray)
{
	IndexBoundingBox meshSize(nMeshNodes.i-1, nMeshNodes.j-1, nMeshNodes.k-1);
	for (size_t index1D : solidNodeIndices)
	{
		Vector3_u solidNode = getIndices3D(index1D, nMeshNodes);
		bool solidNodeHasFluidNeighbor { false };
		if (solidNode.i > meshSize.iMin)
			if (nodeTypeArray(solidNode.i - 1, solidNode.j, solidNode.k) == NodeTypeEnum::FluidActive)
				solidNodeHasFluidNeighbor = true;

		if (solidNode.i < meshSize.iMax)
			if (nodeTypeArray(solidNode.i + 1, solidNode.j, solidNode.k) == NodeTypeEnum::FluidActive)
				solidNodeHasFluidNeighbor = true;

		if (solidNode.j > meshSize.jMin)
			if (nodeTypeArray(solidNode.i, solidNode.j - 1, solidNode.k) == NodeTypeEnum::FluidActive)
				solidNodeHasFluidNeighbor = true;

		if (solidNode.j < meshSize.jMax)
			if (nodeTypeArray(solidNode.i, solidNode.j + 1, solidNode.k) == NodeTypeEnum::FluidActive)
				solidNodeHasFluidNeighbor = true;

		if (solidNode.k > meshSize.kMin)
			if (nodeTypeArray(solidNode.i, solidNode.j, solidNode.k - 1) == NodeTypeEnum::FluidActive)
				solidNodeHasFluidNeighbor = true;

		if (solidNode.k < meshSize.kMax)
			if (nodeTypeArray(solidNode.i, solidNode.j, solidNode.k + 1) == NodeTypeEnum::FluidActive)
				solidNodeHasFluidNeighbor = true;

		if (solidNodeHasFluidNeighbor)
		{
			ghostNodes.emplace_back(solidNode);
			ghostNodeMap[index1D] = &ghostNodes.back();
			nodeTypeArray(solidNode) = NodeTypeEnum::Ghost;
		}
	}
}

void ImmersedBoundary::checkIfSurroundingShouldBeGhost(const Vector3_u &surroundingNode,
													   vector<GhostNode>& newGhostNodes,
													   Array3D_nodeType& nodeTypeArray)
{
	if (nodeTypeArray(surroundingNode) == NodeTypeEnum::Solid)
	{
		newGhostNodes.emplace_back(surroundingNode);
		nodeTypeArray(surroundingNode) = NodeTypeEnum::Ghost;
	}
}

vector<GhostNode> ImmersedBoundary::setImagePointPositions(GhostNodeVectorIterator firstGhostToProcess,
														   const Vector3_d& gridSpacing,
														   const Vector3_u& nMeshNodes,
														   Array3D_nodeType& nodeTypeArray)
{
	vector<GhostNode> newGhostNodes;
	for (GhostNodeVectorIterator ghostIterator{firstGhostToProcess}; ghostIterator!=ghostNodes.end(); ++ghostIterator)
	{
		GhostNode& ghostNode = *ghostIterator;
		Vector3_d ghostNodePosition = getNodePosition(ghostNode.indices, gridSpacing);
		Vector3_d normalProbe = getNormalProbe(ghostNodePosition); // from ghost to body intercept point
		ghostNode.bodyInterceptPoint = ghostNodePosition + normalProbe;
		ghostNode.imagePoint = ghostNode.bodyInterceptPoint + normalProbe;
		IndexBoundingBox surroundingNodes = getSurroundingNodesBox(ghostNode.imagePoint, gridSpacing);
		for( size_t surroundingNodeIndex1D : surroundingNodes.asIndexList(nMeshNodes) )
			checkIfSurroundingShouldBeGhost( getIndices3D(surroundingNodeIndex1D, nMeshNodes),
											 newGhostNodes, nodeTypeArray );
	}
	return newGhostNodes;
}

GhostNodeVectorIterator ImmersedBoundary::appendGhostNodes(const vector<GhostNode>& newGhostNodes,
														   const Vector3_u& nMeshNodes)
{
	for(const GhostNode& newGhostNode : newGhostNodes)
	{
		ghostNodes.push_back(newGhostNode);
		size_t index1D = getIndex1D( newGhostNode.indices, nMeshNodes );
		ghostNodeMap[index1D] = &ghostNodes.back();
	}
	return ghostNodes.end() - newGhostNodes.size();
}

void ImmersedBoundary::setInterpolationValuesFluidNode(uint counter, size_t surroundingNodeIndex1D,
													   const Vector3_u& nMeshNodes, const Vector3_d& gridSpacing,
													   const AllFlowVariablesArrayGroup &flowVariables,
													   InterpolationValues& interpolationValues,		// <- OUTPUT
													   InterpolationPositions& interpolationPositions)	// <- OUTPUT
{
	interpolationValues.u[counter] = flowVariables.primitiveVariables.u(surroundingNodeIndex1D);
	interpolationValues.v[counter] = flowVariables.primitiveVariables.v(surroundingNodeIndex1D);
	interpolationValues.w[counter] = flowVariables.primitiveVariables.w(surroundingNodeIndex1D);
	interpolationValues.p[counter] = flowVariables.primitiveVariables.p(surroundingNodeIndex1D);
	interpolationValues.T[counter] = flowVariables.primitiveVariables.T(surroundingNodeIndex1D);
	Vector3_d surroundingNodePosition = getNodePosition( getIndices3D(surroundingNodeIndex1D, nMeshNodes), gridSpacing );
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
		const Vector3_u& nMeshNodes,
		const Vector3_d& gridSpacing,
		const Array3D_nodeType& nodeTypeArray,
		const AllFlowVariablesArrayGroup& flowVariables,
		InterpolationValues& interpolationValues,		// <-
		InterpolationPositions& interpolationPositions,	// <-
		Array8_b& ghostFlag,							// <- Output
		bool& allSurroundingAreFluid,					// <-
		vector<Vector3_d>& unitNormals)					// <-
{
	uint counter { 0 }; // 0, 1, ..., 7.  To index the interpolation point arrays.
	for (size_t surroundingNodeIndex1D : surroundingNodes.asIndexList(nMeshNodes))
	{
		if (nodeTypeArray(surroundingNodeIndex1D) == NodeTypeEnum::FluidActive
		||	nodeTypeArray(surroundingNodeIndex1D) == NodeTypeEnum::FluidEdge)
		{
			setInterpolationValuesFluidNode(counter, surroundingNodeIndex1D, nMeshNodes,	// <- Input
											gridSpacing, flowVariables,						// <- Input
											interpolationValues, interpolationPositions);	// <- Output
		}
		else if (nodeTypeArray(surroundingNodeIndex1D) == NodeTypeEnum::Ghost)
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

double ImmersedBoundary::simplifiedInterpolation(const Vector8_d& interpolationValues,
												 const Vector3_u& lowerIndexNode,
												 const Vector3_d& imagePointPosition,
												 const Vector3_d& gridSpacing)
{
	Vector3_d lowerIndexNodePosition = getNodePosition(lowerIndexNode, gridSpacing);
	Vector3_d relativePosition = imagePointPosition - lowerIndexNodePosition;

	double variableAt01 = interpolationValues(0)
			+ (interpolationValues(1)-interpolationValues(0)) * relativePosition.z / gridSpacing.z;
	double variableAt23 = interpolationValues(2)
			+ (interpolationValues(3)-interpolationValues(2)) * relativePosition.z / gridSpacing.z;
	double variableAt45 = interpolationValues(4)
			+ (interpolationValues(5)-interpolationValues(4)) * relativePosition.z / gridSpacing.z;
	double variableAt67 = interpolationValues(6)
			+ (interpolationValues(7)-interpolationValues(6)) * relativePosition.z / gridSpacing.z;

	double variableAt0123 = variableAt01 + (variableAt23-variableAt01) * relativePosition.y / gridSpacing.y;
	double variableAt4567 = variableAt45 + (variableAt67-variableAt45) * relativePosition.y / gridSpacing.y;

	return variableAt0123 + (variableAt4567-variableAt0123) * relativePosition.x / gridSpacing.x;
}

PrimitiveVariablesScalars ImmersedBoundary::simplifiedInterpolationAll(
		const InterpolationValues& interpolationValues,
		const Vector3_u& lowerIndexNode,
		const Vector3_d& imagePointPosition,
		const Vector3_d& gridSpacing )
{
	double uInterpolated = simplifiedInterpolation(interpolationValues.u, lowerIndexNode, imagePointPosition, gridSpacing);
	double vInterpolated = simplifiedInterpolation(interpolationValues.v, lowerIndexNode, imagePointPosition, gridSpacing);
	double wInterpolated = simplifiedInterpolation(interpolationValues.w, lowerIndexNode, imagePointPosition, gridSpacing);
	double pInterpolated = simplifiedInterpolation(interpolationValues.p, lowerIndexNode, imagePointPosition, gridSpacing);
	double TInterpolated = simplifiedInterpolation(interpolationValues.T, lowerIndexNode, imagePointPosition, gridSpacing);
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
		const Vector3_d& gridSpacing)
{
	PrimitiveVariablesScalars imagePointPrimVars;
	if (allSurroundingAreFluid)
	{ // Then we can use the simplified interpolation method:
		Vector3_u lowerIndexNode(surroundingNodes.iMin, surroundingNodes.jMin, surroundingNodes.kMin);
		imagePointPrimVars = simplifiedInterpolationAll(interpolationValues,
														lowerIndexNode,
														ghostNode.imagePoint,
														gridSpacing);
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

void CylinderBody::identifyRelatedNodes(const ConfigSettings& params,
		   	   	   	   	   	   	   	    const Vector3_d& gridSpacing,
										const Vector3_u& nMeshNodes,
										Array3D_nodeType& nodeTypeArray	// <- Output
										)
{
	IndexBoundingBox indicesToCheck = getCylinderBoundingBox(gridSpacing, nMeshNodes);
	vector<size_t> solidNodeIndices;
	getSolidNodesInCylinder(params, indicesToCheck, gridSpacing, nMeshNodes, solidNodeIndices, nodeTypeArray);
	findGhostNodesWithFluidNeighbors(solidNodeIndices, nMeshNodes, nodeTypeArray);
	vector<GhostNode> newGhostNodes = setImagePointPositions(ghostNodes.begin(), gridSpacing, nMeshNodes, nodeTypeArray);
	while(newGhostNodes.size() > 0)
	{
		GhostNodeVectorIterator nextToProcess = appendGhostNodes(newGhostNodes, nMeshNodes);
		newGhostNodes = setImagePointPositions(nextToProcess, gridSpacing, nMeshNodes, nodeTypeArray);
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

IndexBoundingBox CylinderBody::getCylinderBoundingBox(const Vector3_d& gridSpacing,
													  const Vector3_u& nMeshNodes) const
{
	size_t indexRadiusX { static_cast<size_t>(ceil(radius / gridSpacing.x)) + 1 };
	size_t indexRadiusY { static_cast<size_t>(ceil(radius / gridSpacing.y)) + 1 };
	size_t indexRadiusZ { static_cast<size_t>(ceil(radius / gridSpacing.z)) + 1 };
	size_t centroidClosestIndexX { static_cast<size_t>(round(centroidPosition.x / gridSpacing.x)) };
	size_t centroidClosestIndexY { static_cast<size_t>(round(centroidPosition.y / gridSpacing.y)) };
	size_t centroidClosestIndexZ { static_cast<size_t>(round(centroidPosition.z / gridSpacing.z)) };
	IndexBoundingBox indicesToCheck(nMeshNodes.i-1, nMeshNodes.j-1, nMeshNodes.k-1);
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

void CylinderBody::getSolidNodesInCylinder(const ConfigSettings& params,
										   const IndexBoundingBox& indicesToCheck,
										   const Vector3_d& gridSpacing,
										   const Vector3_u& nMeshNodes,
										   vector<size_t>& solidNodeIndices, // <- Output
										   Array3D_nodeType& nodeTypeArray   // <- Output
										   )
{
	for (size_t i { indicesToCheck.iMin }; i <= indicesToCheck.iMax; ++i)
		for (size_t j { indicesToCheck.jMin }; j <= indicesToCheck.jMax; ++j)
			for (size_t k { indicesToCheck.kMin }; k <= indicesToCheck.kMax; ++k)
			{
				double distanceFromCentroid{0};
				Vector3_d nodePosition = getNodePosition(i, j, k, gridSpacing);
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
					nodeTypeArray(i, j, k) = NodeTypeEnum::Solid;
					solidNodeIndices.push_back( getIndex1D(i, j, k, nMeshNodes) );
				}
			}
}

SphereBody::SphereBody(Vector3_d centerPosition, double radius) :
centerPosition(centerPosition),
radius{radius}
{}

void SphereBody::identifyRelatedNodes(const ConfigSettings& params,
	   	   	   	    				  const Vector3_d& gridSpacing,
									  const Vector3_u& nMeshNodes,
									  Array3D_nodeType& nodeTypeArray	// <- Output
									  )
{
	IndexBoundingBox indicesToCheck = getSphereBoundingBox(gridSpacing);
	vector<size_t> solidNodeIndices;
	getSolidNodesInSphere(params, indicesToCheck, gridSpacing, nMeshNodes, solidNodeIndices, nodeTypeArray);
	findGhostNodesWithFluidNeighbors(solidNodeIndices, nMeshNodes, nodeTypeArray);

	vector<GhostNode> newGhostNodes = setImagePointPositions(ghostNodes.begin(), gridSpacing, nMeshNodes, nodeTypeArray);
	while(newGhostNodes.size() > 0)
	{
		appendGhostNodes(newGhostNodes, nMeshNodes);
		newGhostNodes = setImagePointPositions(newGhostNodes.begin(), gridSpacing, nMeshNodes, nodeTypeArray);
	}
}

Vector3_d SphereBody::getNormalProbe(const Vector3_d& ghostNodePosition)
{
	Vector3_d centroidToGhost = ghostNodePosition - centerPosition;
	double lengthFactor = (radius - centroidToGhost.length()) / centroidToGhost.length();
	Vector3_d normalProbe = centroidToGhost * lengthFactor;
	return normalProbe;
}

IndexBoundingBox SphereBody::getSphereBoundingBox(const Vector3_d& gridSpacing) const
{
	size_t indexRadiusX { static_cast<size_t>(ceil(radius / gridSpacing.x)) + 1 };
	size_t indexRadiusY { static_cast<size_t>(ceil(radius / gridSpacing.y)) + 1 };
	size_t indexRadiusZ { static_cast<size_t>(ceil(radius / gridSpacing.z)) + 1 };
	size_t centerClosestIndexX { static_cast<size_t>(round(centerPosition.x / gridSpacing.x)) };
	size_t centerClosestIndexY { static_cast<size_t>(round(centerPosition.y / gridSpacing.y)) };
	size_t centerClosestIndexZ { static_cast<size_t>(round(centerPosition.z / gridSpacing.z)) };
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
		   	   	   	   	   	   	   	   const IndexBoundingBox& indicesToCheck,
									   const Vector3_d& gridSpacing,
									   const Vector3_u& nMeshNodes,
									   vector<size_t>& solidNodeIndices, // <- Output
									   Array3D_nodeType& nodeTypeArray   // <- Output
		   	   	   	   	   	   	   	   )
{
	for (size_t i { indicesToCheck.iMin }; i <= indicesToCheck.iMax; ++i)
		for (size_t j { indicesToCheck.jMin }; j <= indicesToCheck.jMax; ++j)
			for (size_t k { indicesToCheck.kMin }; k <= indicesToCheck.kMax; ++k)
			{
				Vector3_d nodePosition = getNodePosition(i, j, k, gridSpacing);
				double distanceFromCenter = sqrt( pow(nodePosition.x - centerPosition.x, 2)
												+ pow(nodePosition.y - centerPosition.y, 2)
												+ pow(nodePosition.z - centerPosition.z, 2));
				if (distanceFromCenter < radius - params.machinePrecisionBuffer)
				{
					nodeTypeArray(i, j, k) = NodeTypeEnum::Solid;
					solidNodeIndices.push_back( getIndex1D(i, j, k, nMeshNodes) );
				}
			}
}





