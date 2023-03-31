/*
 * Boundary.cpp
 *
 *  Created on: Oct 13, 2022
 *      Author: frederk
 */

#include "Boundary.h"

// Constructor. Should be called only by derived classes.
MeshEdgeBoundary::MeshEdgeBoundary(AxisOrientationEnum normalAxis,
								   EdgeIndexEnum planeIndex,
								   NodeTypeEnum ownedNodesType,
								   const SubMeshDescriptor& subMeshData)
: normalAxis{normalAxis},
  planeIndex{planeIndex},
  ownedNodesType{ownedNodesType},
  subMeshData(subMeshData)
{}

// Take an index bounding-box and peel off the specified outer (boundary) plane.
// Returns the off-peeled plane as a "flat" IndexBoundingBox.
IndexBoundingBox MeshEdgeBoundary::peelOffBoundaryPlane(IndexBoundingBox& inputBox,
														AxisOrientationEnum normalAxis,
														EdgeIndexEnum planeIndex)
{
	IndexBoundingBox boundaryPlane = inputBox;
	switch (normalAxis)
	{
	case AxisOrientationEnum::x:
		if(planeIndex == EdgeIndexEnum::min)
		{
			boundaryPlane.iMax = boundaryPlane.iMin;
			inputBox.iMin++;
		}
		else if(planeIndex == EdgeIndexEnum::max)
		{
			boundaryPlane.iMin = boundaryPlane.iMax;
			inputBox.iMax--;
		}
		else throw std::logic_error("Unexpected enum value");
		break;
	case AxisOrientationEnum::y:
		if(planeIndex == EdgeIndexEnum::min)
		{
			boundaryPlane.jMax = boundaryPlane.jMin;
			inputBox.jMin++;
		}
		else if(planeIndex == EdgeIndexEnum::max)
		{
			boundaryPlane.jMin = boundaryPlane.jMax;
			inputBox.jMax--;
		}
		else throw std::logic_error("Unexpected enum value");
		break;
	case AxisOrientationEnum::z:
		if(planeIndex == EdgeIndexEnum::min)
		{
			boundaryPlane.kMax = boundaryPlane.kMin;
			inputBox.kMin++;
		}
		else if(planeIndex == EdgeIndexEnum::max)
		{
			boundaryPlane.kMin = boundaryPlane.kMax;
			inputBox.kMax--;
		}
		else throw std::logic_error("Unexpected enum value");
		break;
	default:
		throw std::logic_error("Unexpected enum value");
	}
	return boundaryPlane;
}

// Claim all unclaimed nodes in the boundary plane, and assign them the appropriate type.
void MeshEdgeBoundary::identifyRelatedNodes(IndexBoundingBox& unclaimedNodes,
											Array3D<NodeTypeEnum>& nodeTypeArray,
											const Vector3_i& regionID,
											const Array3D<SubMeshDescriptor>& neighborSubMeshes)
{
	relatedNodes.regular = peelOffBoundaryPlane(unclaimedNodes, normalAxis, planeIndex);
	for(int i{relatedNodes.regular.iMin}; i<=relatedNodes.regular.iMax; ++i)
		for(int j{relatedNodes.regular.jMin}; j<=relatedNodes.regular.jMax; ++j)
			for(int k{relatedNodes.regular.kMin}; k<=relatedNodes.regular.kMax; ++k)
				nodeTypeArray(i,j,k) = ownedNodesType;
}

// Given a node in the boundary plane, find the two adjacent nodes in the normal direction
void MeshEdgeBoundary::getAdjacentIndices(const Vector3_i& boundaryNode, 		// ← Input
										  int& boundaryAdjacentIndex, 			// ← Output
										  int& nextToAdjacentIndex)				// ← Output
{
	Vector3_i inwardNormalVector = -1*getOutwardNormal();
	boundaryAdjacentIndex = getIndex1D(boundaryNode +   inwardNormalVector, subMeshData.arrayLimits);
	nextToAdjacentIndex   = getIndex1D(boundaryNode + 2*inwardNormalVector, subMeshData.arrayLimits);
}

// Given a node in the boundary plane, find the boundary adjacent node on the opposite side
Vector3_i MeshEdgeBoundary::getPeriodicIndex(const Vector3_i& boundaryNode)
{
	const Vector3_i& nNodes = subMeshData.arrayLimits.nNodes();
	Vector3_i shiftedIndex = boundaryNode + nNodes + getOutwardNormal();
	return shiftedIndex % nNodes;
}

Vector3_i MeshEdgeBoundary::getOutwardNormal() const
{
	Vector3_i normalVector(0,0,0);
	switch(normalAxis)
	{
	case AxisOrientationEnum::x:
		normalVector.i = static_cast<int>(planeIndex);
		break;
	case AxisOrientationEnum::y:
		normalVector.j = static_cast<int>(planeIndex);
		break;
	case AxisOrientationEnum::z:
		normalVector.k = static_cast<int>(planeIndex);
		break;
	}
	return normalVector;
}

// Constructor, taking in the target inlet velocity (after ramp-up).
InletBoundary::InletBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex, const SubMeshDescriptor& subMeshData, double velocity)
: MeshEdgeBoundary(normalAxis, planeIndex, NodeTypeEnum::FluidEdge, subMeshData),
  velocity{velocity}
{}

// Linearly extrapolate density, then filter it. Set normal velocity component, respecting ramp-up.
// Set other velocity components, and temperature fluctuation to zero. Dervive the rest.
// The ramp-up starts at zero when t=0, and then gradually increases to prescribed value.
void InletBoundary::applyBoundaryCondition(double t, const ConfigSettings& params)
{
	double inletVelocity = min(1., t/50.) * velocity; // TODO: move magic const to params as 'rampUpDuration' etc.
	if(planeIndex == EdgeIndexEnum::max)	// If we're on the highest index, velocity must be
		inletVelocity *= -1;				// negative, for this to be an inlet.

	double u{0}, v{0}, w{0};
	if(normalAxis == AxisOrientationEnum::x)
		u = inletVelocity;
	else if(normalAxis == AxisOrientationEnum::y)
		v = inletVelocity;
	else if(normalAxis == AxisOrientationEnum::z)
		w = inletVelocity;
	else
		throw std::logic_error("Unexpected enum value");

	for(int i=relatedNodes.regular.iMin; i<=relatedNodes.regular.iMax; ++i)
		for(int j=relatedNodes.regular.jMin; j<=relatedNodes.regular.jMax; ++j)
			for(int k=relatedNodes.regular.kMin; k<=relatedNodes.regular.kMax; ++k)
			{
				Vector3_i boundaryNode(i,j,k);
				int boundaryAdjacentIndex{0}, nextToAdjacentIndex{0};
				getAdjacentIndices(boundaryNode, boundaryAdjacentIndex, nextToAdjacentIndex);
				double rho = 2*subMeshData.flowVariables.conservedVariables.rho(boundaryAdjacentIndex) // Linear extrapolation
							 - subMeshData.flowVariables.conservedVariables.rho(nextToAdjacentIndex);
				double T = 0;
				double p = (T*(1+rho) + rho) / params.Gamma;
				double uTripFlow=u;
				if(params.Re>50 && j<params.NJ/2 && t<50)
					uTripFlow *=0.9;
				else if(params.Re>50 && j<params.NJ/2 && t<60)
					uTripFlow = u - (60.-t)/10.*0.1*u;
				PrimitiveVariablesScalars primitiveVarsScalars(uTripFlow, v, w, p, T);
				ConservedVariablesScalars   conservedVarsScalars = deriveConservedVariables (primitiveVarsScalars, params);
				TransportPropertiesScalars transportPropsScalars = deriveTransportProperties(primitiveVarsScalars, params);
				setFlowVariablesAtNode(boundaryNode, conservedVarsScalars, primitiveVarsScalars, transportPropsScalars, subMeshData.flowVariables);
			}
}

// Constructor, taking no extra arguments
OutletBoundary::OutletBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex, const SubMeshDescriptor& subMeshData)
: MeshEdgeBoundary(normalAxis, planeIndex, NodeTypeEnum::FluidEdge, subMeshData)
{}

// Set pressure fluctuation to zero. Linearly extrapolate density and momentum. Derive rest.
void OutletBoundary::applyBoundaryCondition(double t, const ConfigSettings& params)
{
	for(int i=relatedNodes.regular.iMin; i<=relatedNodes.regular.iMax; ++i)
		for(int j=relatedNodes.regular.jMin; j<=relatedNodes.regular.jMax; ++j)
			for(int k=relatedNodes.regular.kMin; k<=relatedNodes.regular.kMax; ++k)
			{
				Vector3_i boundaryNode(i,j,k);
				int boundaryAdjacentIndex{0}, nextToAdjacentIndex{0};
				getAdjacentIndices(boundaryNode, boundaryAdjacentIndex, nextToAdjacentIndex);
				double rho = 2*subMeshData.flowVariables.conservedVariables.rho(boundaryAdjacentIndex) // Linear extrapolation
							 - subMeshData.flowVariables.conservedVariables.rho(nextToAdjacentIndex);
				double rho_u = 2*subMeshData.flowVariables.conservedVariables.rho_u(boundaryAdjacentIndex) // Linear extrapolation
							   - subMeshData.flowVariables.conservedVariables.rho_u(nextToAdjacentIndex);
				double rho_v = 2*subMeshData.flowVariables.conservedVariables.rho_v(boundaryAdjacentIndex) // Linear extrapolation
							   - subMeshData.flowVariables.conservedVariables.rho_v(nextToAdjacentIndex);
				double rho_w = 2*subMeshData.flowVariables.conservedVariables.rho_w(boundaryAdjacentIndex) // Linear extrapolation
							   - subMeshData.flowVariables.conservedVariables.rho_w(nextToAdjacentIndex);
				double p = 0;

				double u = rho_u / (1+rho);
				double v = rho_v / (1+rho);
				double w = rho_w / (1+rho);
				double T = ( params.Gamma * p - rho ) / ( 1 + rho );
				double rho_E = p / ( params.Gamma - 1 ) + (1 + rho)/2 * ( u*u + v*v + w*w );

				PrimitiveVariablesScalars primitiveVarsScalars(u, v, w, p, T);
				ConservedVariablesScalars conservedVarsScalars(rho, rho_u, rho_v, rho_w, rho_E);
				TransportPropertiesScalars transportPropsScalars = deriveTransportProperties(primitiveVarsScalars, params);
				setFlowVariablesAtNode(boundaryNode, conservedVarsScalars, primitiveVarsScalars, transportPropsScalars, subMeshData.flowVariables);
			}
}

// Constructor, taking no extra arguments
PeriodicBoundary::PeriodicBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex, const SubMeshDescriptor& subMeshData)
: MeshEdgeBoundary(normalAxis, planeIndex, NodeTypeEnum::FluidGhost, subMeshData)
{}

// Set all flow variables to the value of the boundary-adjacent node at the opposite side of the mesh.
void PeriodicBoundary::applyBoundaryCondition(double t, const ConfigSettings& params)
{
	for(int i=relatedNodes.regular.iMin; i<=relatedNodes.regular.iMax; ++i)
		for(int j=relatedNodes.regular.jMin; j<=relatedNodes.regular.jMax; ++j)
			for(int k=relatedNodes.regular.kMin; k<=relatedNodes.regular.kMax; ++k)
			{
				Vector3_i boundaryNode(i,j,k);
				Vector3_i periodicNode = getPeriodicIndex(boundaryNode);
				subMeshData.flowVariables.conservedVariables.rho  (boundaryNode) = subMeshData.flowVariables.conservedVariables.rho  (periodicNode);
				subMeshData.flowVariables.conservedVariables.rho_u(boundaryNode) = subMeshData.flowVariables.conservedVariables.rho_u(periodicNode);
				subMeshData.flowVariables.conservedVariables.rho_v(boundaryNode) = subMeshData.flowVariables.conservedVariables.rho_v(periodicNode);
				subMeshData.flowVariables.conservedVariables.rho_w(boundaryNode) = subMeshData.flowVariables.conservedVariables.rho_w(periodicNode);
				subMeshData.flowVariables.conservedVariables.rho_E(boundaryNode) = subMeshData.flowVariables.conservedVariables.rho_E(periodicNode);
				subMeshData.flowVariables.primitiveVariables.u(boundaryNode) = subMeshData.flowVariables.primitiveVariables.u(periodicNode);
				subMeshData.flowVariables.primitiveVariables.v(boundaryNode) = subMeshData.flowVariables.primitiveVariables.v(periodicNode);
				subMeshData.flowVariables.primitiveVariables.w(boundaryNode) = subMeshData.flowVariables.primitiveVariables.w(periodicNode);
				subMeshData.flowVariables.primitiveVariables.p(boundaryNode) = subMeshData.flowVariables.primitiveVariables.p(periodicNode);
				subMeshData.flowVariables.primitiveVariables.T(boundaryNode) = subMeshData.flowVariables.primitiveVariables.T(periodicNode);
				subMeshData.flowVariables.transportProperties.mu   (boundaryNode) = subMeshData.flowVariables.transportProperties.mu   (periodicNode);
				subMeshData.flowVariables.transportProperties.kappa(boundaryNode) = subMeshData.flowVariables.transportProperties.kappa(periodicNode);
			}
}

// Constructor, taking no extra arguments
SymmetryBoundary::SymmetryBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex, const SubMeshDescriptor& subMeshData)
: MeshEdgeBoundary(normalAxis, planeIndex, NodeTypeEnum::FluidGhost, subMeshData)
{}

// Symmetry plane is the adjacent interior nodes. The owned boundary nodes are ghost nodes.
// In these ghost nodes all primitive variables are mirrored from the second adjacent interior nodes,
// (mirrored across the symmetry plane), except the normal velocity component, which changes sign.
void SymmetryBoundary::applyBoundaryCondition(double t, const ConfigSettings& params)
{
	const PrimitiveVariablesArrayGroup& primitiveVars = subMeshData.flowVariables.primitiveVariables;
	for(int i=relatedNodes.regular.iMin; i<=relatedNodes.regular.iMax; ++i)
		for(int j=relatedNodes.regular.jMin; j<=relatedNodes.regular.jMax; ++j)
			for(int k=relatedNodes.regular.kMin; k<=relatedNodes.regular.kMax; ++k)
			{
				Vector3_i boundaryNode(i,j,k);
				int boundaryAdjacentIndex{0}, nextToAdjacentIndex{0};
				getAdjacentIndices(boundaryNode, boundaryAdjacentIndex, nextToAdjacentIndex);
				double u, v, w;
				if(normalAxis == AxisOrientationEnum::x)
				{	// Normal component zero, and other components get zero gradient.
					u = -primitiveVars.u(nextToAdjacentIndex);
					v =  primitiveVars.v(nextToAdjacentIndex);
					w =  primitiveVars.w(nextToAdjacentIndex);
				}
				else if(normalAxis == AxisOrientationEnum::y)
				{
					u =  primitiveVars.u(nextToAdjacentIndex);
					v = -primitiveVars.v(nextToAdjacentIndex);
					w =  primitiveVars.w(nextToAdjacentIndex);
				}
				else if(normalAxis == AxisOrientationEnum::z)
				{
					u =  primitiveVars.u(nextToAdjacentIndex);
					v =  primitiveVars.v(nextToAdjacentIndex);
					w = -primitiveVars.w(nextToAdjacentIndex);
				}
				else
					throw std::logic_error("Unexpected enum value");

				double p = primitiveVars.p(nextToAdjacentIndex);
				double T = primitiveVars.T(nextToAdjacentIndex);
				PrimitiveVariablesScalars  primitiveVarsScalars(u, v, w, p, T);
				ConservedVariablesScalars conservedVarsScalars   = deriveConservedVariables (primitiveVarsScalars, params);
				TransportPropertiesScalars transportPropsScalars = deriveTransportProperties(primitiveVarsScalars, params);
				setFlowVariablesAtNode(boundaryNode, conservedVarsScalars, primitiveVarsScalars, transportPropsScalars, subMeshData.flowVariables);
			}
}

// Constructor, taking no extra arguments
ExtrapolationBoundary::ExtrapolationBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex, const SubMeshDescriptor& subMeshData)
: MeshEdgeBoundary(normalAxis, planeIndex, NodeTypeEnum::FluidGhost, subMeshData)
{}

// Everything is extrapolated
void ExtrapolationBoundary::applyBoundaryCondition(double t, const ConfigSettings& params)
{
	const ConservedVariablesArrayGroup& conservedVars = subMeshData.flowVariables.conservedVariables;
	for(int i=relatedNodes.regular.iMin; i<=relatedNodes.regular.iMax; ++i)
		for(int j=relatedNodes.regular.jMin; j<=relatedNodes.regular.jMax; ++j)
			for(int k=relatedNodes.regular.kMin; k<=relatedNodes.regular.kMax; ++k)
			{
				Vector3_i boundaryNode(i,j,k);
				int boundaryAdjacentIndex{0}, nextToAdjacentIndex{0};
				getAdjacentIndices(boundaryNode, boundaryAdjacentIndex, nextToAdjacentIndex);
				double rho = 2*conservedVars.rho(boundaryAdjacentIndex) // Linear extrapolation
							 - conservedVars.rho(nextToAdjacentIndex);
				double rho_u = 2*conservedVars.rho_u(boundaryAdjacentIndex) // Linear extrapolation
							   - conservedVars.rho_u(nextToAdjacentIndex);
				double rho_v = 2*conservedVars.rho_v(boundaryAdjacentIndex) // Linear extrapolation
							   - conservedVars.rho_v(nextToAdjacentIndex);
				double rho_w = 2*conservedVars.rho_w(boundaryAdjacentIndex) // Linear extrapolation
							   - conservedVars.rho_w(nextToAdjacentIndex);
				double rho_E = 2*conservedVars.rho_E(boundaryAdjacentIndex) // Linear extrapolation
							   - conservedVars.rho_E(nextToAdjacentIndex);
				ConservedVariablesScalars conservedVarsScalars(rho, rho_u, rho_v, rho_w, rho_E);
				PrimitiveVariablesScalars primitiveVarsScalars = derivePrimitiveVariables(conservedVarsScalars, params);
				TransportPropertiesScalars transportPropsScalars = deriveTransportProperties(primitiveVarsScalars, params);
				setFlowVariablesAtNode(boundaryNode, conservedVarsScalars, primitiveVarsScalars, transportPropsScalars, subMeshData.flowVariables);
			}
}

int SubmeshInterfaceBoundary::nExtraNodeLayer()
{
	if(neighborSubMesh.refinementLevel > subMeshData.refinementLevel)
		return 0;
	else if(neighborSubMesh.refinementLevel < subMeshData.refinementLevel)
		return 1;
	else if(planeIndex == EdgeIndexEnum::min)
		return 0;
	else
		return 1;
}

// Check if "otherRegion" applies its boundary conditions before "thisRegion". If so, all nodes in otherRegion are safe to copy.
// Otherwise only its interior nodes are safe.
bool SubmeshInterfaceBoundary::subMeshHasAppliedBC(const Vector3_i& otherRegionID, const Vector3_i& thisRegionID) const
{
	return  otherRegionID.i <  thisRegionID.i
		|| (otherRegionID.i == thisRegionID.i && otherRegionID.j <  thisRegionID.j)
		|| (otherRegionID.i == thisRegionID.i && otherRegionID.j == thisRegionID.j && otherRegionID.k < thisRegionID.k);
}

// Get the index that another submesh would use to access the given node.
// It might be outside the other submesh's array limits.
// It might have a decimal part, if the other submesh is coarser than 'thisSubMesh'.
Vector3_d SubmeshInterfaceBoundary::indexInOtherSubMesh(const Vector3_i& node,
														const SubMeshDescriptor& otherSubMesh,
														const SubMeshDescriptor& thisSubMesh) const
{
	Vector3_i neighborDirection = otherSubMesh.regionID - thisSubMesh.regionID;
	Vector3_d indexInOther(0,0,0);
	double levelFactor = pow(2, otherSubMesh.refinementLevel - thisSubMesh.refinementLevel);
	if(otherSubMesh.regionID.i < thisSubMesh.regionID.i)
		indexInOther.x = otherSubMesh.nNodes.i-1 + node.i * levelFactor;
	else if(otherSubMesh.regionID.i > thisSubMesh.regionID.i)
		indexInOther.x = (-thisSubMesh.nNodes.i-1 + node.i) * levelFactor;
	else
		indexInOther.x = node.i * levelFactor;

	if(otherSubMesh.regionID.j < thisSubMesh.regionID.j)
		indexInOther.y = otherSubMesh.nNodes.j-1 + node.j * levelFactor;
	else if(otherSubMesh.regionID.j > thisSubMesh.regionID.j)
		indexInOther.y = (-thisSubMesh.nNodes.j-1 + node.j) * levelFactor;
	else
		indexInOther.y = node.j * levelFactor;

	if(otherSubMesh.regionID.k < thisSubMesh.regionID.k)
		indexInOther.z = otherSubMesh.nNodes.k-1 + node.k * levelFactor;
	else if(otherSubMesh.regionID.k > thisSubMesh.regionID.k)
		indexInOther.z = (-thisSubMesh.nNodes.k-1 + node.k) * levelFactor;
	else
		indexInOther.z = node.k * levelFactor;

	return indexInOther;
}

// Get an index vector that tells which of the neighbor submeshes own a specified node.
// "Owning" here means that it's safe to copy that node value from that submesh.
// "node" indices should be with respect to the submesh that has the boundary that is calling this func.
// Returned vector: regionID of the submesh where it's safe to copy from.
Vector3_i SubmeshInterfaceBoundary::whoOwnsThisNode(const Vector3_i& node,
													const Vector3_i& regionID,
													const Array3D<SubMeshDescriptor>& subMeshes) const
{
	const SubMeshDescriptor& thisSubMesh = subMeshes(regionID);
	if( node.i > thisSubMesh.arrayLimits.iMin
	&&	node.i < thisSubMesh.arrayLimits.iMax
	&&	node.j > thisSubMesh.arrayLimits.jMin
	&&	node.j < thisSubMesh.arrayLimits.jMax
	&&	node.k > thisSubMesh.arrayLimits.kMin
	&&	node.k < thisSubMesh.arrayLimits.kMax )
		return regionID;

	for(const SubMeshDescriptor& loopSubMesh : subMeshes)
	{
		Vector3_d indexInOther = indexInOtherSubMesh(node, loopSubMesh, thisSubMesh);
		if( subMeshHasAppliedBC(loopSubMesh.regionID, thisSubMesh.regionID) )
		{
			if(		indexInOther.x >= loopSubMesh.arrayLimits.iMin
				&&	indexInOther.x <= loopSubMesh.arrayLimits.iMax
				&&	indexInOther.y >= loopSubMesh.arrayLimits.jMin
				&&	indexInOther.y <= loopSubMesh.arrayLimits.jMax
				&&	indexInOther.z >= loopSubMesh.arrayLimits.kMin
				&&	indexInOther.z <= loopSubMesh.arrayLimits.kMax )
				return loopSubMesh.regionID;
		}
		else
		{
			if(		indexInOther.x > loopSubMesh.arrayLimits.iMin
				&&	indexInOther.x < loopSubMesh.arrayLimits.iMax
				&&	indexInOther.y > loopSubMesh.arrayLimits.jMin
				&&	indexInOther.y < loopSubMesh.arrayLimits.jMax
				&&	indexInOther.z > loopSubMesh.arrayLimits.kMin
				&&	indexInOther.z < loopSubMesh.arrayLimits.kMax )
				return loopSubMesh.regionID;
		}
	}
	std::stringstream ss;
	ss 		<< "Could not find any submesh with a valid value for node "
			<< node.i << "," << node.j << "," << node.k
			<< " in mesh "
			<< regionID.i << "," << regionID.j << "," << regionID.k
			<< endl;
	throw std::logic_error( ss.str() );
	return Vector3_i(-1,-1,-1);
}

// Get the node(s) surrounding a point. It can be 1, 2, 4, or 8 nodes, depending on how many grid lines the point coincides with.
vector<int> SubmeshInterfaceBoundary::getNodesAroundPoint(const Vector3_d& point, const IndexBoundingBox& arrayLimits) const
{
	IndexBoundingBox interpolationNodes( static_cast<int>(std::floor(point.x)),
										 static_cast<int>(std::ceil (point.x)),
										 static_cast<int>(std::floor(point.y)),
										 static_cast<int>(std::ceil (point.y)),
										 static_cast<int>(std::floor(point.z)),
										 static_cast<int>(std::ceil (point.z))
										 );
	return interpolationNodes.allNodesAsIndexList(arrayLimits);
}

void SubmeshInterfaceBoundary::identifyRelatedNodes(IndexBoundingBox& unclaimedNodes,
													Array3D<NodeTypeEnum>& nodeTypeArray,
													const Vector3_i& regionID,
													const Array3D<SubMeshDescriptor>& subMeshes)
{
	const SubMeshDescriptor& thisSubMesh = subMeshes(regionID);
	IndexBoundingBox allOwnedNodes = peelOffBoundaryPlane(unclaimedNodes, normalAxis, planeIndex);
	for(int i{allOwnedNodes.iMin}; i<=allOwnedNodes.iMax; ++i)
		for(int j{allOwnedNodes.jMin}; j<=allOwnedNodes.jMax; ++j)
			for(int k{allOwnedNodes.kMin}; k<=allOwnedNodes.kMax; ++k)
				nodeTypeArray(i,j,k) = ownedNodesType;

	IndexBoundingBox regularNodes = allOwnedNodes;
	if(normalAxis != AxisOrientationEnum::x)
	{
		regularNodes.iMin = 1;
		regularNodes.iMax = thisSubMesh.nNodes.i - 2;
	}
	if(normalAxis != AxisOrientationEnum::y)
	{
		regularNodes.jMin = 1;
		regularNodes.jMax = thisSubMesh.nNodes.j - 2;
	}
	if(normalAxis != AxisOrientationEnum::z)
	{
		regularNodes.kMin = 1;
		regularNodes.kMax = thisSubMesh.nNodes.k - 2;
	}
	relatedNodes.regular = regularNodes;

	neighborSubMesh = subMeshes( regionID + getOutwardNormal() );

	for( int index1D : allOwnedNodes.asIndexListExcept(regularNodes, thisSubMesh.arrayLimits) )
	{
		Vector3_i node = getIndices3D(index1D, thisSubMesh.arrayLimits);
		Vector3_i borrowMeshID = whoOwnsThisNode(node, regionID, subMeshes);
		Vector3_d indicesInOther = indexInOtherSubMesh(node, subMeshes(borrowMeshID), thisSubMesh);
		vector<int> nodesInOther = getNodesAroundPoint(indicesInOther, subMeshes(borrowMeshID).arrayLimits);
		relatedNodes.specialTreatment.push_back( SpecialTreatmentNodeInfo( index1D,
																		   nodesInOther,
																		   subMeshes(borrowMeshID).flowVariables
																		   ) );
	}
}

double SubmeshInterfaceBoundary::getMeanNodeValue(vector<int> nodes, const Array3D<double>& flowVariable) const
{
	vector<double> nodeValues = flowVariable(nodes);
	double sumOfNodeValues = 0;
	for(double nodeValue : nodeValues)
		sumOfNodeValues += nodeValue;
	double meanNodeValue = sumOfNodeValues / nodeValues.size();
	return meanNodeValue;
}


void SubmeshInterfaceBoundary::applyBoundaryCondition(double t, const ConfigSettings& params)
{
	ConservedVariablesArrayGroup& conservedVarsNeighbor = neighborSubMesh.flowVariables.conservedVariables;
	for(int i=relatedNodes.regular.iMin; i<=relatedNodes.regular.iMax; i+=2)
		for(int j=relatedNodes.regular.jMin; j<=relatedNodes.regular.jMax; j+=2)
			for(int k=relatedNodes.regular.kMin; k<=relatedNodes.regular.kMax; k+=2)
			{
				Vector3_i node(i,j,k);
				Vector3_d indicesInOther = indexInOtherSubMesh(node, neighborSubMesh, subMeshData);
				vector<int> nodesInOther = getNodesAroundPoint(indicesInOther, neighborSubMesh.arrayLimits);
				double rho   = getMeanNodeValue(nodesInOther, conservedVarsNeighbor.rho  );
				double rho_u = getMeanNodeValue(nodesInOther, conservedVarsNeighbor.rho_u);
				double rho_v = getMeanNodeValue(nodesInOther, conservedVarsNeighbor.rho_v);
				double rho_w = getMeanNodeValue(nodesInOther, conservedVarsNeighbor.rho_w);
				double rho_E = getMeanNodeValue(nodesInOther, conservedVarsNeighbor.rho_E);
				ConservedVariablesScalars consVarsThisNode(rho, rho_u, rho_v, rho_w, rho_E);
				PrimitiveVariablesScalars primVarsThisNode = derivePrimitiveVariables(consVarsThisNode, params);
				TransportPropertiesScalars transPropsThisNode = deriveTransportProperties(primVarsThisNode, params);
				setFlowVariablesAtNode(node, consVarsThisNode, primVarsThisNode, transPropsThisNode, subMeshData.flowVariables);
			}
}

//////////////////////////////////////
// Mesh edge BCs above ↑
//////////////////////////////////////
// Immersed BCs below ↓
//////////////////////////////////////

// Get interpolation values, to compute values in image points, and then the ghost nodes.
// Also filters density and energy in the closest fluid nodes.
void ImmersedBoundary::applyBoundaryCondition(const ConfigSettings& params)
{
	for(GhostNode ghostNode : ghostNodes)
	{
		InterpolationValues interpolationValues;		// Conserved variables at the 8 surrounding interpolation points.
		InterpolationPositions interpolationPositions;	// Positions of interpolation points
		Array8_b ghostFlag;					// A bool flag for each surrounding node. True if ghost.
		ghostFlag.fill(false);				// ← Sets all 8 values to false
		bool allSurroundingAreFluid{true};	// Until otherwise proven
		vector<Vector3_d> unitNormals;		// For each of the surrounding nodes that is ghost, we put its unit normal probe here.
		IndexBoundingBox surroundingNodes = getSurroundingNodesBox(ghostNode.imagePoint, subMeshData.spacings, subMeshData.boundingBox.getMinPoint(), subMeshData.arrayLimits);
		
		setInterpolationValues(
				surroundingNodes,									// ← Input
				interpolationValues, interpolationPositions,		// ← Output
				ghostFlag, allSurroundingAreFluid, unitNormals);	// ← Output
		
		PrimitiveVariablesScalars imagePointPrimVars = interpolateImagePointVariables(
				interpolationValues, interpolationPositions, allSurroundingAreFluid,
				ghostFlag, ghostNode, surroundingNodes, unitNormals);

		PrimitiveVariablesScalars ghostNodePrimVars = getGhostNodeBCVariables(imagePointPrimVars);
		ConservedVariablesScalars ghostNodeConsVars = deriveConservedVariables(ghostNodePrimVars, params);
		TransportPropertiesScalars ghostNodeTransportProps = deriveTransportProperties(ghostNodePrimVars, params);
		setFlowVariablesAtNode(ghostNode.indices, ghostNodeConsVars, ghostNodePrimVars, ghostNodeTransportProps, subMeshData.flowVariables);
	}
}

// Given a vector of indices to solid nodes, find all who will be part of the numerical stencil
// of any fluid node. Those solid nodes are marked as solidGhost, and added to the ghostNodes vector.
void ImmersedBoundary::findGhostNodesWithFluidNeighbors(const vector<int>& solidNodeIndices)
{
	Array3D<NodeTypeEnum>& nodeType = subMeshData.nodeType;
	const IndexBoundingBox& arrayLimits = subMeshData.arrayLimits;
	for (int index1D : solidNodeIndices)
	{
		// In 3D, each node has 26 neighbors. All except the 8 corner nodes (i+1,j+1,k+1), (i+1,j+1,k-1), etc,
		// are part of the 2nd order stencil for the N-S equations. Thus we must check 18 nodes:
		Vector3_i solidNode = getIndices3D(index1D, subMeshData.arrayLimits);
		bool solidNodeHasFluidNeighbor { false };

		if (solidNode.i > arrayLimits.iMin) // Then we can test all nodes at i-1:
		{
			if (nodeType(solidNode.i - 1, solidNode.j, solidNode.k) == NodeTypeEnum::FluidInterior)
				solidNodeHasFluidNeighbor = true; // (-1,0,0)

			if (solidNode.j > arrayLimits.jMin)
				if (nodeType(solidNode.i - 1, solidNode.j - 1, solidNode.k) == NodeTypeEnum::FluidInterior)
					solidNodeHasFluidNeighbor = true; // (-1,-1,0)

			if (solidNode.j < arrayLimits.jMax)
				if (nodeType(solidNode.i - 1, solidNode.j + 1, solidNode.k) == NodeTypeEnum::FluidInterior)
					solidNodeHasFluidNeighbor = true; // (-1,1,0)

			if (solidNode.k > arrayLimits.kMin)
				if (nodeType(solidNode.i - 1, solidNode.j, solidNode.k - 1) == NodeTypeEnum::FluidInterior)
					solidNodeHasFluidNeighbor = true; // (-1,0,-1)

			if (solidNode.k < arrayLimits.kMax)
				if (nodeType(solidNode.i - 1, solidNode.j, solidNode.k + 1) == NodeTypeEnum::FluidInterior)
					solidNodeHasFluidNeighbor = true; // (-1,0,1)
		}

		if (solidNode.i < arrayLimits.iMax) // Then we can test all nodes at i+1:
		{
			if (nodeType(solidNode.i + 1, solidNode.j, solidNode.k) == NodeTypeEnum::FluidInterior)
				solidNodeHasFluidNeighbor = true; // (1,0,0)

			if (solidNode.j > arrayLimits.jMin)
				if (nodeType(solidNode.i + 1, solidNode.j - 1, solidNode.k) == NodeTypeEnum::FluidInterior)
					solidNodeHasFluidNeighbor = true; // (1,-1,0)

			if (solidNode.j < arrayLimits.jMax)
				if (nodeType(solidNode.i + 1, solidNode.j + 1, solidNode.k) == NodeTypeEnum::FluidInterior)
					solidNodeHasFluidNeighbor = true; // (1,1,0)

			if (solidNode.k > arrayLimits.kMin)
				if (nodeType(solidNode.i + 1, solidNode.j, solidNode.k - 1) == NodeTypeEnum::FluidInterior)
					solidNodeHasFluidNeighbor = true; // (1,0,-1)

			if (solidNode.k < arrayLimits.kMax)
				if (nodeType(solidNode.i + 1, solidNode.j, solidNode.k + 1) == NodeTypeEnum::FluidInterior)
					solidNodeHasFluidNeighbor = true; // (1,0,1)
		}

		// And then we check the nodes at the same i-index (x-coordinate) as the solid node:

		if (solidNode.j > arrayLimits.jMin)
		{
			if (nodeType(solidNode.i, solidNode.j - 1, solidNode.k) == NodeTypeEnum::FluidInterior)
				solidNodeHasFluidNeighbor = true; // (0,-1,0

			if (solidNode.k > arrayLimits.kMin)
				if (nodeType(solidNode.i, solidNode.j - 1, solidNode.k - 1) == NodeTypeEnum::FluidInterior)
					solidNodeHasFluidNeighbor = true; // (0,-1,-1)

			if (solidNode.k < arrayLimits.kMax)
				if (nodeType(solidNode.i, solidNode.j - 1, solidNode.k + 1) == NodeTypeEnum::FluidInterior)
					solidNodeHasFluidNeighbor = true; // (0,-1,1)
		}

		if (solidNode.j < arrayLimits.jMax)
		{
			if (nodeType(solidNode.i, solidNode.j + 1, solidNode.k) == NodeTypeEnum::FluidInterior)
				solidNodeHasFluidNeighbor = true; // (0,1,0)

			if (solidNode.k > arrayLimits.kMin)
				if (nodeType(solidNode.i, solidNode.j + 1, solidNode.k - 1) == NodeTypeEnum::FluidInterior)
					solidNodeHasFluidNeighbor = true; // (0,1,-1)

			if (solidNode.k < arrayLimits.kMax)
				if (nodeType(solidNode.i, solidNode.j + 1, solidNode.k + 1) == NodeTypeEnum::FluidInterior)
					solidNodeHasFluidNeighbor = true; // (0,1,1)
		}

		if (solidNode.k > arrayLimits.kMin)
			if (nodeType(solidNode.i, solidNode.j, solidNode.k - 1) == NodeTypeEnum::FluidInterior)
				solidNodeHasFluidNeighbor = true; // (0,0,-1)

		if (solidNode.k < arrayLimits.kMax)
			if (nodeType(solidNode.i, solidNode.j, solidNode.k + 1) == NodeTypeEnum::FluidInterior)
				solidNodeHasFluidNeighbor = true; // (0,0,1)

		if (solidNodeHasFluidNeighbor)
		{
			ghostNodeMap[index1D] = ghostNodes.size();
			ghostNodes.emplace_back(solidNode);
			nodeType(solidNode) = NodeTypeEnum::SolidGhost;
		}
	}
}

// Given one of the 8 nodes surrounding an image point, make sure it is not inactive.
// Mostly the 'findGhostNodesWithFluidNeighbors' function will find all solids that should
// be ghost, but when the immersed boundary intersects a mesh edge boundary (and possibly
// other cases), an image point can be surrounded by a solid that's not part of a stencil.
void ImmersedBoundary::checkIfSurroundingShouldBeGhost(const Vector3_i &surroundingNode,
													   vector<GhostNode>& newGhostNodes)
{
	if (subMeshData.nodeType(surroundingNode) == NodeTypeEnum::SolidInactive)
	{
		newGhostNodes.emplace_back(surroundingNode);
		subMeshData.nodeType(surroundingNode) = NodeTypeEnum::SolidGhost;
	}
}

// Given an iterator to a GhostNode in a vector, go through the ghosts after the iterator, and find
// their image points (IP). Also check if any surrounding node of the IP are solidInactive. If so,
// mark them as ghost ad add them to the return-vector.
vector<GhostNode> ImmersedBoundary::setImagePointPositions(GhostNodeVectorIterator firstGhostToProcess)
{
	vector<GhostNode> newGhostNodes;
	for (GhostNodeVectorIterator ghostIterator{firstGhostToProcess}; ghostIterator!=ghostNodes.end(); ++ghostIterator)
	{
		GhostNode& ghostNode = *ghostIterator;
		Vector3_d ghostNodePosition = getNodePosition( ghostNode.indices, subMeshData.arrayLimits, subMeshData.spacings, subMeshData.boundingBox.getMinPoint() );
		Vector3_d normalProbe = getNormalProbe(ghostNodePosition); // from ghost to body intercept point
		ghostNode.bodyInterceptPoint = ghostNodePosition + normalProbe;
		ghostNode.imagePoint = ghostNode.bodyInterceptPoint + normalProbe;
		IndexBoundingBox surroundingNodes = getSurroundingNodesBox(ghostNode.imagePoint, subMeshData.spacings,
																	subMeshData.boundingBox.getMinPoint(), subMeshData.arrayLimits);
		for( int surroundingNodeIndex1D : surroundingNodes.cornersAsIndexList(IndexBoundingBox()) )
			checkIfSurroundingShouldBeGhost( getIndices3D(surroundingNodeIndex1D, subMeshData.arrayLimits),
											 newGhostNodes );
		// Because of the way that lift/drag etc is computed we need also BI to be only surrounded by ghost or fluid:
		surroundingNodes = getSurroundingNodesBox(ghostNode.bodyInterceptPoint, subMeshData.spacings, subMeshData.boundingBox.getMinPoint(), subMeshData.arrayLimits);
		for( int surroundingNodeIndex1D : surroundingNodes.cornersAsIndexList(subMeshData.arrayLimits) )
			checkIfSurroundingShouldBeGhost( getIndices3D(surroundingNodeIndex1D, subMeshData.arrayLimits),
											 newGhostNodes );
	}
	return newGhostNodes;
}

// Add the nodes from newGhostNodes at the end of ghostNodes. Return iterator to the first appended node.
GhostNodeVectorIterator ImmersedBoundary::appendGhostNodes(const vector<GhostNode>& newGhostNodes)
{
	for(const GhostNode& newGhostNode : newGhostNodes)
	{
		int index1D = getIndex1D( newGhostNode.indices, subMeshData.arrayLimits );
		ghostNodeMap[index1D] = ghostNodes.size();
		ghostNodes.push_back(newGhostNode);
	}
	return ghostNodes.end() - newGhostNodes.size();
}

// Set interpolation values for one of the 8 nodes surrounding an image point (IP). Surrounding node must be fluid.
void ImmersedBoundary::setInterpolationValuesFluidNode(int counter,	// ← Which of the 8 surrounding nodes [0-7]
													   int surroundingNodeIndex1D,
													   InterpolationValues& interpolationValues,		// ← OUTPUT
													   InterpolationPositions& interpolationPositions)	// ← OUTPUT
{
	interpolationValues.u[counter] = subMeshData.flowVariables.primitiveVariables.u(surroundingNodeIndex1D);
	interpolationValues.v[counter] = subMeshData.flowVariables.primitiveVariables.v(surroundingNodeIndex1D);
	interpolationValues.w[counter] = subMeshData.flowVariables.primitiveVariables.w(surroundingNodeIndex1D);
	interpolationValues.p[counter] = subMeshData.flowVariables.primitiveVariables.p(surroundingNodeIndex1D);
	interpolationValues.T[counter] = subMeshData.flowVariables.primitiveVariables.T(surroundingNodeIndex1D);
	Vector3_d surroundingNodePosition = getNodePosition( getIndices3D(surroundingNodeIndex1D, subMeshData.arrayLimits), subMeshData.arrayLimits,
															subMeshData.spacings, subMeshData.boundingBox.getMinPoint() );
	interpolationPositions.x[counter] = surroundingNodePosition.x;
	interpolationPositions.y[counter] = surroundingNodePosition.y;
	interpolationPositions.z[counter] = surroundingNodePosition.z;
}

// Set interpolation values for one of the 8 nodes surrounding an image point (IP). Surrounding node must be ghost.
void ImmersedBoundary::setInterpolationValuesGhostNode(
		int counter,	// ← Which of the 8 surrounding nodes [0-7]
		int surroundingNodeIndex1D,
		vector<Vector3_d>& unitNormals,					// ↰
		InterpolationValues& interpolationValues,		// ← OUTPUT
		InterpolationPositions& interpolationPositions)	// ↲
{
	GhostNode &surroundingGhostNode = ghostNodes.at( ghostNodeMap.at(surroundingNodeIndex1D) );
	Vector3_d &unitNormal = unitNormals.emplace_back();
	unitNormal = surroundingGhostNode.imagePoint - surroundingGhostNode.bodyInterceptPoint;
	unitNormal = unitNormal / unitNormal.length();
	interpolationValues.u[counter] = 0;
	interpolationValues.v[counter] = 0;
	interpolationValues.w[counter] = 0; // NB! This hardcodes zero Dirichlet condition for velocity
	interpolationValues.p[counter] = 0; // and zero gradient Neumann conditions for pressure and temperature!
	interpolationValues.T[counter] = 0;
	interpolationPositions.x[counter] = surroundingGhostNode.bodyInterceptPoint.x;
	interpolationPositions.y[counter] = surroundingGhostNode.bodyInterceptPoint.y;
	interpolationPositions.z[counter] = surroundingGhostNode.bodyInterceptPoint.z;
}

// Get interpolation values and positions of the 8 surrounding nodes of an image point (IP).
// Also check if any of the 8 nodes are ghost. If so, note which ones and their normal vectors.
void ImmersedBoundary::setInterpolationValues(
		const IndexBoundingBox& surroundingNodes,
		InterpolationValues& interpolationValues,		// ↰
		InterpolationPositions& interpolationPositions,	// ↰
		Array8_b& ghostFlag,							// ← Output
		bool& allSurroundingAreFluid,					// ↲
		vector<Vector3_d>& unitNormals)					// ↲
{
	int counter { 0 }; // 0, 1, ..., 7.  To index the interpolation point arrays.
	for ( int surroundingNodeIndex1D : surroundingNodes.cornersAsIndexList(subMeshData.arrayLimits) )
	{
		if (subMeshData.nodeType(surroundingNodeIndex1D) == NodeTypeEnum::FluidInterior
		||	subMeshData.nodeType(surroundingNodeIndex1D) == NodeTypeEnum::FluidEdge
		||	subMeshData.nodeType(surroundingNodeIndex1D) == NodeTypeEnum::FluidGhost)
		{
			setInterpolationValuesFluidNode(counter, surroundingNodeIndex1D,				// ← Input
											interpolationValues, interpolationPositions);	// ← Output
		}
		else if (subMeshData.nodeType(surroundingNodeIndex1D) == NodeTypeEnum::SolidGhost)
		{
			ghostFlag[counter] = true;
			allSurroundingAreFluid = false;
			setInterpolationValuesGhostNode(counter, surroundingNodeIndex1D,	// ← Input
					unitNormals, interpolationValues, interpolationPositions);	// ← Output
		}
		else
			throw std::logic_error("Impossible situation. Found solid node around image point.");

		++counter;
	}
}

// Get image point value of a flow variable, using the fast simplified interpolation.
// Only valid if all 8 surrounding nodes are fluid.
double ImmersedBoundary::simplifiedInterpolation(const Vector8_d& interpolationValues,
												 const Vector3_i& lowerIndexNode, // ← Node with smallest coordinates
												 const Vector3_d& imagePointPosition)
{
	Vector3_d lowerIndexNodePosition = getNodePosition( lowerIndexNode, subMeshData.arrayLimits, subMeshData.spacings, subMeshData.boundingBox.getMinPoint() );
	Vector3_d relativePosition = imagePointPosition - lowerIndexNodePosition;

	// The ordering [0-7] is such that z(k) changes fastest, x(i) changes the least.
	// See member function 'asIndexList' in IndexBoundingBox.
	// First, interpolate along the 4 lines between nodes with low and high z-value:
	double variableAt01 = interpolationValues(0)
			+ (interpolationValues(1)-interpolationValues(0)) * relativePosition.z / gridSpacing.z;
	double variableAt23 = interpolationValues(2)
			+ (interpolationValues(3)-interpolationValues(2)) * relativePosition.z / gridSpacing.z;
	double variableAt45 = interpolationValues(4)
			+ (interpolationValues(5)-interpolationValues(4)) * relativePosition.z / gridSpacing.z;
	double variableAt67 = interpolationValues(6)
			+ (interpolationValues(7)-interpolationValues(6)) * relativePosition.z / gridSpacing.z;
	// With the 4 resulting values, interpolate in the y-direction:
	double variableAt0123 = variableAt01 + (variableAt23-variableAt01) * relativePosition.y / gridSpacing.y;
	double variableAt4567 = variableAt45 + (variableAt67-variableAt45) * relativePosition.y / gridSpacing.y;
	// With the 2 resulting values, interpolate in x-direction to get final value:
	return variableAt0123 + (variableAt4567-variableAt0123) * relativePosition.x / gridSpacing.x;
}

// Get image point values of conserved variables, using the fast simplified interpolation.
// Only valid if all 8 surrounding nodes are fluid.
PrimitiveVariablesScalars ImmersedBoundary::simplifiedInterpolationAll(
		const InterpolationValues& interpolationValues,
		const Vector3_i& lowerIndexNode,
		const Vector3_d& imagePointPosition,
		const Vector3_d& gridSpacing,
		const Vector3_d& meshOriginOffset )
{
	double u = simplifiedInterpolation(interpolationValues.u, lowerIndexNode, imagePointPosition, gridSpacing, meshOriginOffset);
	double v = simplifiedInterpolation(interpolationValues.v, lowerIndexNode, imagePointPosition, gridSpacing, meshOriginOffset);
	double w = simplifiedInterpolation(interpolationValues.w, lowerIndexNode, imagePointPosition, gridSpacing, meshOriginOffset);
	double p = simplifiedInterpolation(interpolationValues.p, lowerIndexNode, imagePointPosition, gridSpacing, meshOriginOffset);
	double T = simplifiedInterpolation(interpolationValues.T, lowerIndexNode, imagePointPosition, gridSpacing, meshOriginOffset);
	return PrimitiveVariablesScalars(u, v, w, p, T);
}

// Given the conserved variables in an image point (IP), get the values for the ghost node.
PrimitiveVariablesScalars ImmersedBoundary::getGhostNodeBCVariables(const PrimitiveVariablesScalars& imagePointBCVars)
{
	double uGhost = -imagePointBCVars.u;	// NB! This hardcodes zero Dirichlet condition for velocity
	double vGhost = -imagePointBCVars.v;	// and zero gradient Neumann conditions for pressure and temperature.
	double wGhost = -imagePointBCVars.w;
	double pGhost =  imagePointBCVars.p;
	double TGhost =  imagePointBCVars.T;
	return PrimitiveVariablesScalars(uGhost, vGhost, wGhost, pGhost, TGhost);
}

// Compute the elements of the Vandermonde matrix for Dirichlet conditions.
void ImmersedBoundary::populateVandermondeDirichlet(const InterpolationPositions& points, Matrix8x8_d& vandermonde)
{
	vandermonde << Vector8_d::Constant(1), points.x, points.y, points.z, points.x*points.y, points.x*points.z, points.y*points.z, points.x*points.y*points.z;
}

// Compute the elements of the Vandermonde matrix for Neumann conditions.
void ImmersedBoundary::populateVandermondeNeumann(const InterpolationPositions& points,
												  const Array8_b& ghostFlags,
												  const vector<Vector3_d>& unitNormals,
												  Matrix8x8_d& vandermonde)
{
	// Start with the same matrix as for Dirichlet, and change the rows with ghost nodes:
	populateVandermondeDirichlet(points, vandermonde);
	int counter{0};
	for(int rowIndex{0}; rowIndex<8; ++rowIndex)
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

// Get a value in the image point, given the values in the 8 surrounding nodes, and the Vandermonde matrix.
double ImmersedBoundary::trilinearInterpolation(const Vector8_d& interpolationValues,
												const Vector3_d& imagePoint,
												const Matrix8x8_d& vandermonde)
{
	// Solve the 8x8 linear system to get the coefficients, and use them in the interpolation function:
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

// Compute all conserved variables in an image point, given the values in the 8 surrounding nodes,
// and the Vandermonde matrices for Dirichlet and Neumann conditions.
PrimitiveVariablesScalars ImmersedBoundary::trilinearInterpolationAll(const InterpolationValues& values,
																	  const Vector3_d& imagePoint,
																	  const Matrix8x8_d& vandermondeDirichlet,
																	  const Matrix8x8_d& vandermondeNeumann)
{
	double u = trilinearInterpolation(values.u, imagePoint, vandermondeDirichlet);
	double v = trilinearInterpolation(values.v, imagePoint, vandermondeDirichlet);
	double w = trilinearInterpolation(values.w, imagePoint, vandermondeDirichlet);
	double p = trilinearInterpolation(values.p, imagePoint, vandermondeNeumann);
	double T = trilinearInterpolation(values.T, imagePoint, vandermondeNeumann);
	return PrimitiveVariablesScalars(u, v, w, p, T);
}

// Compute variables in an image point using simplified interpolation if possible,
// or the Vandermonde matrix if necessary.
PrimitiveVariablesScalars ImmersedBoundary::interpolateImagePointVariables(
		const InterpolationValues& interpolationValues,
		const InterpolationPositions& interpolationPositions,
		bool allSurroundingAreFluid,
		const Array8_b& ghostFlag,
		const GhostNode& ghostNode,
		const IndexBoundingBox& surroundingNodes,
		const vector<Vector3_d>& unitNormals )
{
	PrimitiveVariablesScalars imagePointBCVars;
	if (allSurroundingAreFluid)
	{ // Then we can use the simplified interpolation method:
		Vector3_i lowerIndexNode(surroundingNodes.iMin, surroundingNodes.jMin, surroundingNodes.kMin);
		imagePointBCVars = simplifiedInterpolationAll(interpolationValues,
														lowerIndexNode,
														ghostNode.imagePoint,
														mesh.spacings, mesh.originOffset);
	}
	else
	{ // Then we must use trilinear interpolation with the Vandermonde matrix:
		Matrix8x8_d vandermondeDirichlet, vandermondeNeumann;
		populateVandermondeDirichlet(interpolationPositions, vandermondeDirichlet);
		populateVandermondeNeumann(interpolationPositions, ghostFlag, unitNormals, vandermondeNeumann);
		imagePointBCVars = trilinearInterpolationAll(interpolationValues,
													   ghostNode.imagePoint,
													   vandermondeDirichlet,
													   vandermondeNeumann);
	}
	return imagePointBCVars;
}

// Constructor. Only two of the coordinates in centroidPosition matter. Not the one corresponding
// to the cylinder axis direction.
CylinderBody::CylinderBody(Vector3_d centroidPosition, AxisOrientationEnum axis, double radius, const SubMeshDescriptor& subMeshData) :
ImmersedBoundary(subMeshData),
centroidPosition(centroidPosition),
axis{axis},
radius{radius}
{}

// Check if cylinder centroid passes through box
bool CylinderBody::isInsideBox(const SpaceBoundingBox& box)
{
	return
			( axis == AxisOrientationEnum::x && centroidPosition.y < box.yMax && centroidPosition.y > box.yMin
											 && centroidPosition.z < box.zMax && centroidPosition.z > box.zMin )
		||	( axis == AxisOrientationEnum::y && centroidPosition.x < box.xMax && centroidPosition.x > box.xMin
											 && centroidPosition.z < box.zMax && centroidPosition.z > box.zMin )
		||	( axis == AxisOrientationEnum::z && centroidPosition.x < box.xMax && centroidPosition.x > box.xMin
											 && centroidPosition.y < box.yMax && centroidPosition.y > box.yMin ) ;
}

// Mark solid nodes (ghost and inactive), find image points (IP) and closest fluid nodes (to filter).
void CylinderBody::identifyRelatedNodes(const ConfigSettings& params,
		   	   	   	   	   	   	   	    const Vector3_d& gridSpacing,
										const Vector3_i& nMeshNodes,
										const Vector3_d& meshOriginOffset,
										Array3D<NodeTypeEnum>& nodeTypeArray	// ← Output
										)
{
	int filterNodesLayerWidth = 2; // TODO: Remove if it seems to be OK without this extra filtering. Otherwise: param?
	IndexBoundingBox indicesToCheck = getCylinderBoundingBox(gridSpacing, nMeshNodes, filterNodesLayerWidth);
	vector<int> solidNodeIndices;
	getSolidAndFilterNodesInCylinder(params, indicesToCheck, gridSpacing, nMeshNodes,
			meshOriginOffset, filterNodesLayerWidth, solidNodeIndices, nodeTypeArray);
	findGhostNodesWithFluidNeighbors(solidNodeIndices, nMeshNodes, nodeTypeArray);
	// At this stage, we have probably found all the ghost nodes. However, when we compute the IP positions we must
	// still ensure that the 8 nodes surrounding the IP are not inactive. This is done iteratively until there are
	// no new ghost nodes:
	vector<GhostNode> newGhostNodes = setImagePointPositions(ghostNodes.begin(), gridSpacing, nMeshNodes,
																meshOriginOffset, nodeTypeArray);
	while(newGhostNodes.size() > 0)
	{
		GhostNodeVectorIterator nextToProcess = appendGhostNodes(newGhostNodes, nMeshNodes);
		newGhostNodes = setImagePointPositions(nextToProcess, gridSpacing, nMeshNodes, meshOriginOffset, nodeTypeArray);
	}
}

// Compute integral properties as lift, drag, separation angle.
IntegralProperties CylinderBody::getIntegralProperties(const ConfigSettings& params, const SubMeshDescriptor& mesh)
{
	struct BI_info
	{
		Vector3_d BI;
		Eigen::Vector2d traction;
		Eigen::Vector2d normal;
	};
	std::map<double, BI_info> angleToBIMap;
	double pi = std::atan2(0, -1);
	for(const GhostNode& ghostNode : ghostNodes)
	{
		if(ghostNode.indices.k != 1)
			continue;
		Vector3_d BI = ghostNode.bodyInterceptPoint;
		IndexBoundingBox surroundingNodesOfBI = getSurroundingNodesBox(BI, mesh.spacings,
																		mesh.originOffset, mesh.nNodes);
		if( ceil(BI.y/mesh.spacings.y) == BI.y/mesh.spacings.y )
			std::cout << "WARNING: BI point coincides with mesh node. This will bias the lift/drag calculation." << endl;
		InterpolationValues interpolationValues;
		InterpolationPositions interpolationPositions;
		Array8_b ghostFlag;
		ghostFlag.fill(false);
		bool allSurroundingAreFluid{true};
		vector<Vector3_d> unitNormals;
		setInterpolationValues(surroundingNodesOfBI, mesh, interpolationValues, interpolationPositions,
								ghostFlag, allSurroundingAreFluid, unitNormals);
		Matrix8x8_d vandermondeDirichlet, vandermondeNeumann;
		populateVandermondeDirichlet(interpolationPositions, vandermondeDirichlet);
		populateVandermondeNeumann(interpolationPositions, ghostFlag, unitNormals, vandermondeNeumann);
		Vector8_d aCoeffU = vandermondeDirichlet.fullPivLu().solve( interpolationValues.u.matrix() );
		Vector8_d aCoeffV = vandermondeDirichlet.fullPivLu().solve( interpolationValues.v.matrix() );
		double p = trilinearInterpolation(interpolationValues.p, BI, vandermondeNeumann);
		double T = trilinearInterpolation(interpolationValues.T, BI, vandermondeNeumann);
		double ScPlusOne = 1 + params.sutherlands_C2 / params.T_0;
		double mu = pow( 1+T, 1.5 ) * ScPlusOne / ( params.Re_0*( T + ScPlusOne ) );
		double dudx = aCoeffU(1) + aCoeffU(4)*BI.y + aCoeffU(5)*BI.z + aCoeffU(7)*BI.y*BI.z;
		double dvdx = aCoeffV(1) + aCoeffV(4)*BI.y + aCoeffV(5)*BI.z + aCoeffV(7)*BI.y*BI.z;
		double dudy = aCoeffU(2) + aCoeffU(4)*BI.x + aCoeffU(6)*BI.z + aCoeffU(7)*BI.x*BI.z;
		double dvdy = aCoeffV(2) + aCoeffV(4)*BI.x + aCoeffV(6)*BI.z + aCoeffV(7)*BI.x*BI.z;
		Eigen::Matrix2d deformation{ {dudx, (dudy+dvdx)/2}, {(dudy+dvdx)/2, dvdy} };
		Eigen::Matrix2d viscousStress;
		viscousStress = 2*mu*deformation - 2./3.*mu*(dudx+dvdy)*Eigen::Matrix2d::Identity();
		Eigen::Matrix2d totalStress;
		totalStress = viscousStress - (p + 1./params.Gamma)*Eigen::Matrix2d::Identity();
		Eigen::Vector2d normal({ (ghostNode.imagePoint-BI).x, (ghostNode.imagePoint-BI).y });
		normal.normalize();
		Eigen::Vector2d traction = totalStress * normal;
		Vector3_d cylinderCentroid( params.L_x / 4, params.L_y / 2, 0 );
		Vector3_d relativePos = BI - cylinderCentroid;
		double angle = std::atan2(relativePos.y, relativePos.x);
		if (angle < 0)
			angle += 2*pi;
		angleToBIMap[angle].BI = BI;
		angleToBIMap[angle].traction = traction;
		angleToBIMap[angle].normal = normal;
	}
	IntegralProperties integralProps;
	double prevAngle = ( angleToBIMap.rbegin() )->first - 2*pi;
	BI_info prevBI = ( angleToBIMap.rbegin() )->second;
	Eigen::Vector2d prevUnitTangent({ -prevBI.normal(1), prevBI.normal(0) });
	double prevWallShearStress = prevBI.traction.transpose() * prevUnitTangent;
	Eigen::Vector2d force({0,0});

	for (const std::pair<const double, BI_info>& currentPair : angleToBIMap)
	{
		double thisAngle = currentPair.first;
		const BI_info& thisBI = currentPair.second;
		Eigen::Vector2d unitTangent({ -thisBI.normal(1), thisBI.normal(0) });
		double wallShearStress = thisBI.traction.transpose() * unitTangent;
		if(std::signbit(wallShearStress) != std::signbit(prevWallShearStress) )
		{
			double separationAngle = prevWallShearStress * (thisAngle-prevAngle) / (prevWallShearStress-wallShearStress) + prevAngle;
			integralProps.separationAngles.push_back(separationAngle);
		}
		double thisSurfaceLength = ( thisBI.BI - prevBI.BI ).length();
		Eigen::Vector2d forceContribution = ( prevBI.traction + thisBI.traction ) / 2 * thisSurfaceLength;
		force += forceContribution;
		prevAngle = thisAngle;
		prevWallShearStress = wallShearStress;
		prevBI = thisBI;
	}
	integralProps.drag = force(0);
	integralProps.lift = force(1);
	return integralProps;
}

// Get vector from given ghost node position to the closest point on the surface.
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

// Get index bounding box that is guaranteed to enclose the cylinder and the layer of fluid nodes to filter.
IndexBoundingBox CylinderBody::getCylinderBoundingBox(const Vector3_d& gridSpacing,
													  const Vector3_i& nMeshNodes,
													  int filterNodesLayerWidth) const
{
	int indexRadiusX { static_cast<int>(ceil(radius / gridSpacing.x)) + 1 + filterNodesLayerWidth };
	int indexRadiusY { static_cast<int>(ceil(radius / gridSpacing.y)) + 1 + filterNodesLayerWidth };
	int indexRadiusZ { static_cast<int>(ceil(radius / gridSpacing.z)) + 1 + filterNodesLayerWidth };
	int centroidClosestIndexX { static_cast<int>(round(centroidPosition.x / gridSpacing.x)) };
	int centroidClosestIndexY { static_cast<int>(round(centroidPosition.y / gridSpacing.y)) };
	int centroidClosestIndexZ { static_cast<int>(round(centroidPosition.z / gridSpacing.z)) };
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

// Check the given part of the mesh, and find the solid nodes and the closest fluid nodes (to filter)
void CylinderBody::getSolidAndFilterNodesInCylinder(const ConfigSettings& params,
										   	   	    const IndexBoundingBox& indicesToCheck,
													const Vector3_d& gridSpacing,
													const Vector3_i& nMeshNodes,
													const Vector3_d& meshOriginOffset,
													int nNodesFilterLayer,
													vector<int>& solidNodeIndices, // ← Output
													Array3D<NodeTypeEnum>& nodeTypeArray   // ← Output
										   	   	    )
{
	for (int i { indicesToCheck.iMin }; i <= indicesToCheck.iMax; ++i)
		for (int j { indicesToCheck.jMin }; j <= indicesToCheck.jMax; ++j)
			for (int k { indicesToCheck.kMin }; k <= indicesToCheck.kMax; ++k)
			{
				double distanceFromCentroid{0};
				Vector3_d nodePosition = getNodePosition(i, j, k, gridSpacing, meshOriginOffset);
				double filterLayerWidth;
				if (axis == AxisOrientationEnum::x)
				{
					distanceFromCentroid = sqrt( pow(nodePosition.y - centroidPosition.y, 2)
											   + pow(nodePosition.z - centroidPosition.z, 2));
					filterLayerWidth = nNodesFilterLayer*max(gridSpacing.y, gridSpacing.z);
				}
				else if (axis == AxisOrientationEnum::y)
				{
					distanceFromCentroid = sqrt( pow(nodePosition.x - centroidPosition.x, 2)
											   + pow(nodePosition.z - centroidPosition.z, 2));
					filterLayerWidth = nNodesFilterLayer*max(gridSpacing.x, gridSpacing.z);
				}
				else if (axis == AxisOrientationEnum::z)
				{
					distanceFromCentroid = sqrt( pow(nodePosition.x - centroidPosition.x, 2)
											   + pow(nodePosition.y - centroidPosition.y, 2));
					filterLayerWidth = nNodesFilterLayer*max(gridSpacing.x, gridSpacing.y);
				}
				else
					throw std::logic_error("Unexpected enum value");
				if (distanceFromCentroid < radius - params.machinePrecisionBuffer)
				{
					nodeTypeArray(i, j, k) = NodeTypeEnum::SolidInactive;
					solidNodeIndices.push_back( getIndex1D(i, j, k, nMeshNodes) );
				}
				else if (distanceFromCentroid < radius + filterLayerWidth
						&& nodeTypeArray(i,j,k) == NodeTypeEnum::FluidInterior)
					filterNodes.push_back( getIndex1D(i, j, k, nMeshNodes) );
			}
}

// Constructor
SphereBody::SphereBody(Vector3_d centerPosition, double radius, const SubMeshDescriptor& subMeshData) :
ImmersedBoundary(subMeshData),
centerPosition(centerPosition),
radius{radius}
{}

// Check if sphere center point is inside the given box:
bool SphereBody::isInsideBox(const SpaceBoundingBox& box)
{
	return	centerPosition.x < box.xMax && centerPosition.x > box.xMin
		&&	centerPosition.y < box.yMax && centerPosition.y > box.yMin
		&&	centerPosition.z < box.zMax && centerPosition.z > box.zMin ;
}

// Mark solid nodes (ghost and inactive), find image points (IP) and closest fluid nodes (to filter).
void SphereBody::identifyRelatedNodes(const ConfigSettings& params,
	   	   	   	    				  const Vector3_d& gridSpacing,
									  const Vector3_i& nMeshNodes,
									  const Vector3_d& meshOriginOffset,
									  Array3D<NodeTypeEnum>& nodeTypeArray	// ← Output
									  )
{
	int filterNodesLayerWidth = 2; // TODO: param
	IndexBoundingBox indicesToCheck = getSphereBoundingBox(gridSpacing, filterNodesLayerWidth);
	vector<int> solidNodeIndices;
	getSolidAndFilterNodesInSphere(params, indicesToCheck, gridSpacing, nMeshNodes, meshOriginOffset, filterNodesLayerWidth, solidNodeIndices, nodeTypeArray);
	findGhostNodesWithFluidNeighbors(solidNodeIndices, nMeshNodes, nodeTypeArray);
	// At this stage, we have probably found all the ghost nodes. However, when we compute the IP positions we must
	// still ensure that the 8 nodes surrounding the IP are not inactive. This is done iteratively until there are
	// no new ghost nodes:
	vector<GhostNode> newGhostNodes = setImagePointPositions(ghostNodes.begin(), gridSpacing, nMeshNodes, meshOriginOffset, nodeTypeArray);
	while(newGhostNodes.size() > 0)
	{
		appendGhostNodes(newGhostNodes, nMeshNodes);
		newGhostNodes = setImagePointPositions(newGhostNodes.begin(), gridSpacing, nMeshNodes, meshOriginOffset, nodeTypeArray);
	}
}

// Compute integral properties as lift, drag, separation angle.
IntegralProperties SphereBody::getIntegralProperties(const ConfigSettings& params, const SubMeshDescriptor& mesh)
{
	struct BI_info
	{
		Vector3_d BI;
		Eigen::Vector2d traction;
		Eigen::Vector2d normal;
	};
	std::map<double, BI_info> angleToBIMap;
	double pi = std::atan2(0, -1);
	for(const GhostNode& ghostNode : ghostNodes)
	{
		if(ghostNode.indices.k != 1)
			continue;
		Vector3_d BI = ghostNode.bodyInterceptPoint;
		IndexBoundingBox surroundingNodesOfBI = getSurroundingNodesBox(BI, mesh.spacings,
																		mesh.originOffset, mesh.nNodes);
		if( ceil(BI.y/mesh.spacings.y) == BI.y/mesh.spacings.y )
			std::cout << "WARNING: BI point coincides with mesh node. This will bias the lift/drag calculation." << endl;
		InterpolationValues interpolationValues;
		InterpolationPositions interpolationPositions;
		Array8_b ghostFlag;
		ghostFlag.fill(false);
		bool allSurroundingAreFluid{true};
		vector<Vector3_d> unitNormals;
		setInterpolationValues(surroundingNodesOfBI, mesh, interpolationValues, interpolationPositions,
								ghostFlag, allSurroundingAreFluid, unitNormals);
		Matrix8x8_d vandermondeDirichlet, vandermondeNeumann;
		populateVandermondeDirichlet(interpolationPositions, vandermondeDirichlet);
		populateVandermondeNeumann(interpolationPositions, ghostFlag, unitNormals, vandermondeNeumann);
		Vector8_d aCoeffU = vandermondeDirichlet.fullPivLu().solve( interpolationValues.u.matrix() );
		Vector8_d aCoeffV = vandermondeDirichlet.fullPivLu().solve( interpolationValues.v.matrix() );
		double p = trilinearInterpolation(interpolationValues.p, BI, vandermondeNeumann);
		double T = trilinearInterpolation(interpolationValues.T, BI, vandermondeNeumann);
		double ScPlusOne = 1 + params.sutherlands_C2 / params.T_0;
		double mu = pow( 1+T, 1.5 ) * ScPlusOne / ( params.Re_0*( T + ScPlusOne ) );
		double dudx = aCoeffU(1) + aCoeffU(4)*BI.y + aCoeffU(5)*BI.z + aCoeffU(7)*BI.y*BI.z;
		double dvdx = aCoeffV(1) + aCoeffV(4)*BI.y + aCoeffV(5)*BI.z + aCoeffV(7)*BI.y*BI.z;
		double dudy = aCoeffU(2) + aCoeffU(4)*BI.x + aCoeffU(6)*BI.z + aCoeffU(7)*BI.x*BI.z;
		double dvdy = aCoeffV(2) + aCoeffV(4)*BI.x + aCoeffV(6)*BI.z + aCoeffV(7)*BI.x*BI.z;
		Eigen::Matrix2d deformation{ {dudx, (dudy+dvdx)/2}, {(dudy+dvdx)/2, dvdy} };
		Eigen::Matrix2d viscousStress;
		viscousStress = 2*mu*deformation - 2./3.*mu*(dudx+dvdy)*Eigen::Matrix2d::Identity();
		Eigen::Matrix2d totalStress;
		totalStress = viscousStress - (p + 1./params.Gamma)*Eigen::Matrix2d::Identity();
		Eigen::Vector2d normal({ (ghostNode.imagePoint-BI).x, (ghostNode.imagePoint-BI).y });
		normal.normalize();
		Eigen::Vector2d traction = totalStress * normal;
		Vector3_d cylinderCentroid( params.L_x / 4, params.L_y / 2, 0 );
		Vector3_d relativePos = BI - cylinderCentroid;
		double angle = std::atan2(relativePos.y, relativePos.x);
		if (angle < 0)
			angle += 2*pi;
		angleToBIMap[angle].BI = BI;
		angleToBIMap[angle].traction = traction;
		angleToBIMap[angle].normal = normal;
	}
	IntegralProperties integralProps;
	double prevAngle = ( angleToBIMap.rbegin() )->first - 2*pi;
	BI_info prevBI = ( angleToBIMap.rbegin() )->second;
	Eigen::Vector2d prevUnitTangent({ -prevBI.normal(1), prevBI.normal(0) });
	double prevWallShearStress = prevBI.traction.transpose() * prevUnitTangent;
	Eigen::Vector2d force({0,0});

	for (const std::pair<const double, BI_info>& currentPair : angleToBIMap)
	{
		double thisAngle = currentPair.first;
		const BI_info& thisBI = currentPair.second;
		Eigen::Vector2d unitTangent({ -thisBI.normal(1), thisBI.normal(0) });
		double wallShearStress = thisBI.traction.transpose() * unitTangent;
		if(std::signbit(wallShearStress) != std::signbit(prevWallShearStress) )
		{
			double separationAngle = prevWallShearStress * (thisAngle-prevAngle) / (prevWallShearStress-wallShearStress) + prevAngle;
			integralProps.separationAngles.push_back(separationAngle);
		}
		double thisSurfaceLength = ( thisBI.BI - prevBI.BI ).length();
		Eigen::Vector2d forceContribution = ( prevBI.traction + thisBI.traction ) / 2 * thisSurfaceLength;
		force += forceContribution;
		prevAngle = thisAngle;
		prevWallShearStress = wallShearStress;
		prevBI = thisBI;
	}
	integralProps.drag = force(0);
	integralProps.lift = force(1);
	return integralProps;
}

// Get vector from given ghost node position to the closest point on the surface.
Vector3_d SphereBody::getNormalProbe(const Vector3_d& ghostNodePosition)
{
	Vector3_d centroidToGhost = ghostNodePosition - centerPosition;
	double lengthFactor = (radius - centroidToGhost.length()) / centroidToGhost.length();
	Vector3_d normalProbe = centroidToGhost * lengthFactor;
	return normalProbe;
}

// Get index bounding box that is guaranteed to enclose the sphere and the layer of fluid nodes to filter.
IndexBoundingBox SphereBody::getSphereBoundingBox(const Vector3_d& gridSpacing, int filterNodeLayerWidth) const
{
	int indexRadiusX { static_cast<int>(ceil(radius / gridSpacing.x)) + 1 + filterNodeLayerWidth };
	int indexRadiusY { static_cast<int>(ceil(radius / gridSpacing.y)) + 1 + filterNodeLayerWidth };
	int indexRadiusZ { static_cast<int>(ceil(radius / gridSpacing.z)) + 1 + filterNodeLayerWidth };
	int centerClosestIndexX { static_cast<int>(round(centerPosition.x / gridSpacing.x)) };
	int centerClosestIndexY { static_cast<int>(round(centerPosition.y / gridSpacing.y)) };
	int centerClosestIndexZ { static_cast<int>(round(centerPosition.z / gridSpacing.z)) };
	IndexBoundingBox indicesToCheck;
	indicesToCheck.iMin = centerClosestIndexX - indexRadiusX;
	indicesToCheck.iMax = centerClosestIndexX + indexRadiusX;
	indicesToCheck.jMin = centerClosestIndexY - indexRadiusY;
	indicesToCheck.jMax = centerClosestIndexY + indexRadiusY;
	indicesToCheck.kMin = centerClosestIndexZ - indexRadiusZ;
	indicesToCheck.kMax = centerClosestIndexZ + indexRadiusZ;
	return indicesToCheck;
}

// Check the given part of the mesh, and find the solid nodes and the closest fluid nodes (to filter)
void SphereBody::getSolidAndFilterNodesInSphere(const ConfigSettings& params,
		   	   	   	   	   	   	   	   	   	    const IndexBoundingBox& indicesToCheck,
												const Vector3_d& gridSpacing,
												const Vector3_i& nMeshNodes,
												const Vector3_d& meshOriginOffset,
												int nNodesFilterLayer,
												vector<int>& solidNodeIndices, // ← Output
												Array3D<NodeTypeEnum>& nodeTypeArray   // ← Output
		   	   	   	   	   	   	   	   	   	    )
{
	double filterLayerWidth = nNodesFilterLayer * max({gridSpacing.x, gridSpacing.y, gridSpacing.z});
	for (int i { indicesToCheck.iMin }; i <= indicesToCheck.iMax; ++i)
		for (int j { indicesToCheck.jMin }; j <= indicesToCheck.jMax; ++j)
			for (int k { indicesToCheck.kMin }; k <= indicesToCheck.kMax; ++k)
			{
				Vector3_d nodePosition = getNodePosition(i, j, k, gridSpacing, meshOriginOffset);
				double distanceFromCenter = sqrt( pow(nodePosition.x - centerPosition.x, 2)
												+ pow(nodePosition.y - centerPosition.y, 2)
												+ pow(nodePosition.z - centerPosition.z, 2));
				if (distanceFromCenter < radius - params.machinePrecisionBuffer)
				{
					nodeTypeArray(i, j, k) = NodeTypeEnum::SolidInactive;
					solidNodeIndices.push_back( getIndex1D(i, j, k, nMeshNodes) );
				}
				else if (distanceFromCenter < radius + filterLayerWidth
						&& nodeTypeArray(i,j,k) == NodeTypeEnum::FluidInterior)
					filterNodes.push_back( getIndex1D(i, j, k, nMeshNodes) );
			}
}





