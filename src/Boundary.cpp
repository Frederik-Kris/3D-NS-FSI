/*
 * Boundary.cpp
 *
 *  Created on: Oct 13, 2022
 *      Author: frederk
 */

#include "Boundary.h"

// Constructor. Should be called only by derived classes.
MeshEdgeBoundary::MeshEdgeBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex, NodeTypeEnum ownedNodesType)
: normalAxis{normalAxis},
  planeIndex{planeIndex},
  ownedNodesType{ownedNodesType}
{}

// Claim all unclaimed nodes in the boundary plane, and assign them the appropriate type.
void MeshEdgeBoundary::identifyOwnedNodes(	IndexBoundingBox& unclaimedNodes,
											const Vector3_i& nMeshNodes,
											Array3D_nodeType& nodeTypeArray )
{
	ownedNodes = unclaimedNodes;
	switch (normalAxis)
	{
	case AxisOrientationEnum::x:
		if(planeIndex == EdgeIndexEnum::min)
		{
			ownedNodes.iMax = ownedNodes.iMin;
			unclaimedNodes.iMin++;
		}
		else if(planeIndex == EdgeIndexEnum::max)
		{
			ownedNodes.iMin = ownedNodes.iMax;
			unclaimedNodes.iMax--;
		}
		else throw std::logic_error("Unexpected enum value");
		break;
	case AxisOrientationEnum::y:
		if(planeIndex == EdgeIndexEnum::min)
		{
			ownedNodes.jMax = ownedNodes.jMin;
			unclaimedNodes.jMin++;
		}
		else if(planeIndex == EdgeIndexEnum::max)
		{
			ownedNodes.jMin = ownedNodes.jMax;
			unclaimedNodes.jMax--;
		}
		else throw std::logic_error("Unexpected enum value");
		break;
	case AxisOrientationEnum::z:
		if(planeIndex == EdgeIndexEnum::min)
		{
			ownedNodes.kMax = ownedNodes.kMin;
			unclaimedNodes.kMin++;
		}
		else if(planeIndex == EdgeIndexEnum::max)
		{
			ownedNodes.kMin = ownedNodes.kMax;
			unclaimedNodes.kMax--;
		}
		else throw std::logic_error("Unexpected enum value");
		break;
	default:
		throw std::logic_error("Unexpected enum value");
	}
	for(int i{ownedNodes.iMin}; i<=ownedNodes.iMax; ++i)
		for(int j{ownedNodes.jMin}; j<=ownedNodes.jMax; ++j)
			for(int k{ownedNodes.kMin}; k<=ownedNodes.kMax; ++k)
				nodeTypeArray(i,j,k) = ownedNodesType;
}

// Given a node in the boundary plane, find the two adjacent nodes in the normal direction
void MeshEdgeBoundary::getAdjacentIndices(int index1D, const Vector3_i& nMeshNodes,	// <- Input
										  int& boundaryAdjacentIndex, 				// <- Output
										  int& nextToAdjacentIndex)					// <- Output
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

// Given a node in the boundary plane, find the boundary adjacent node on the opposite side
int MeshEdgeBoundary::getPeriodicIndex(int index1D, const Vector3_i& nMeshNodes)
{
	Vector3_i indices = getIndices3D(index1D, nMeshNodes);
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

// Constructor, taking in the target inlet velocity (after ramp-up).
InletBoundary::InletBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex, double velocity)
: MeshEdgeBoundary(normalAxis, planeIndex, NodeTypeEnum::FluidEdge),
  velocity{velocity}
{}

// Filter the density in the inlet plane. Does not affect or use any nodes outside this plane.
void InletBoundary::filterInletDensity(const Vector3_i& nMeshNodes,					// <- Input
		   	   	   	   	   	   	   	   const ConfigSettings& params,				// <- Input
									   AllFlowVariablesArrayGroup& flowVariables)	// <- Output
{
	vector<double> rhoFilteredVector;
	vector<int> filteredNodes;
	if(normalAxis == AxisOrientationEnum::x)
		for(int j=ownedNodes.jMin; j<=ownedNodes.jMax; ++j)
			for(int k=ownedNodes.kMin; k<=ownedNodes.kMax; ++k)
			{
				double rhoFiltered = 0;
				int nNeighbors = 0;
				if(j<ownedNodes.jMax)
				{
					rhoFiltered += flowVariables.conservedVariables.rho(ownedNodes.iMin, j+1, k);
					++nNeighbors;
				}
				if(j>ownedNodes.jMin)
				{
					rhoFiltered += flowVariables.conservedVariables.rho(ownedNodes.iMin, j-1, k);
					++nNeighbors;
				}
				if(k<ownedNodes.kMax)
				{
					rhoFiltered += flowVariables.conservedVariables.rho(ownedNodes.iMin, j, k+1);
					++nNeighbors;
				}
				if(k>ownedNodes.kMin)
				{
					rhoFiltered += flowVariables.conservedVariables.rho(ownedNodes.iMin, j, k-1);
					++nNeighbors;
				}
				rhoFiltered += nNeighbors * flowVariables.conservedVariables.rho(ownedNodes.iMin, j, k);
				rhoFiltered /= 2 * nNeighbors;
				rhoFilteredVector.push_back(rhoFiltered);
				filteredNodes.push_back(getIndex1D(ownedNodes.iMin, j, k, nMeshNodes));
			}
	else if(normalAxis == AxisOrientationEnum::y)
		for(int i=ownedNodes.iMin; i<=ownedNodes.iMax; ++i)
			for(int k=ownedNodes.kMin; k<=ownedNodes.kMax; ++k)
			{
				double rhoFiltered = 0;
				int nNeighbors = 0;
				if(i<ownedNodes.iMax)
				{
					rhoFiltered += flowVariables.conservedVariables.rho(i+1, ownedNodes.jMin, k);
					++nNeighbors;
				}
				if(i>ownedNodes.iMin)
				{
					rhoFiltered += flowVariables.conservedVariables.rho(i-1, ownedNodes.jMin, k);
					++nNeighbors;
				}
				if(k<ownedNodes.kMax)
				{
					rhoFiltered += flowVariables.conservedVariables.rho(i, ownedNodes.jMin, k+1);
					++nNeighbors;
				}
				if(k>ownedNodes.kMin)
				{
					rhoFiltered += flowVariables.conservedVariables.rho(i, ownedNodes.jMin, k-1);
					++nNeighbors;
				}
				rhoFiltered += nNeighbors * flowVariables.conservedVariables.rho(i, ownedNodes.jMin, k);
				rhoFiltered /= 2 * nNeighbors;
				rhoFilteredVector.push_back(rhoFiltered);
				filteredNodes.push_back(getIndex1D(i, ownedNodes.jMin, k, nMeshNodes));
			}
	else if(normalAxis == AxisOrientationEnum::z)
		for(int i=ownedNodes.iMin; i<=ownedNodes.iMax; ++i)
			for(int j=ownedNodes.jMin; j<=ownedNodes.jMax; ++j)
			{
				double rhoFiltered = 0;
				int nNeighbors = 0;
				if(i<ownedNodes.iMax)
				{
					rhoFiltered += flowVariables.conservedVariables.rho(i+1, j, ownedNodes.kMin);
					++nNeighbors;
				}
				if(i>ownedNodes.iMin)
				{
					rhoFiltered += flowVariables.conservedVariables.rho(i-1, j, ownedNodes.kMin);
					++nNeighbors;
				}
				if(j<ownedNodes.jMax)
				{
					rhoFiltered += flowVariables.conservedVariables.rho(i, j+1, ownedNodes.kMin);
					++nNeighbors;
				}
				if(j>ownedNodes.jMin)
				{
					rhoFiltered += flowVariables.conservedVariables.rho(i, j-1, ownedNodes.kMin);
					++nNeighbors;
				}
				rhoFiltered += nNeighbors * flowVariables.conservedVariables.rho(i, j, ownedNodes.kMin);
				rhoFiltered /= 2 * nNeighbors;
				rhoFilteredVector.push_back(rhoFiltered);
				filteredNodes.push_back(getIndex1D(i, j, ownedNodes.kMin, nMeshNodes));
			}
	else
		throw std::logic_error("Unexpected enum value");
	for(int index{0}; index<filteredNodes.size(); ++index)
		flowVariables.conservedVariables.rho(filteredNodes[index]) = rhoFilteredVector[index];
}

// Linearly extrapolate density, then filter it. Set normal velocity component, respecting ramp-up.
// Set other velocity components, and temperature fluctuation to zero. Dervive the rest.
// The ramp-up starts at zero when t=0, and then gradually increases to prescribed value.
void InletBoundary::applyBoundaryCondition(double t, const Vector3_i& nMeshNodes,		// <- Input
										   const ConfigSettings& params,				// <- Input
										   AllFlowVariablesArrayGroup& flowVariables)	// <- Output
{
	for(int i=ownedNodes.iMin; i<=ownedNodes.iMax; ++i)
		for(int j=ownedNodes.jMin; j<=ownedNodes.jMax; ++j)
			for(int k=ownedNodes.kMin; k<=ownedNodes.kMax; ++k)
			{
				int index1D = getIndex1D(i, j, k, nMeshNodes);
				int boundaryAdjacentIndex{0}, nextToAdjacentIndex{0};
				getAdjacentIndices(index1D, nMeshNodes, boundaryAdjacentIndex, nextToAdjacentIndex);
				flowVariables.conservedVariables.rho(i,j,k) =
						2*flowVariables.conservedVariables.rho(boundaryAdjacentIndex) // Linear extrapolation
						- flowVariables.conservedVariables.rho(nextToAdjacentIndex);
			}
	filterInletDensity(nMeshNodes, params, flowVariables);

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

	for(int i=ownedNodes.iMin; i<=ownedNodes.iMax; ++i)
		for(int j=ownedNodes.jMin; j<=ownedNodes.jMax; ++j)
			for(int k=ownedNodes.kMin; k<=ownedNodes.kMax; ++k)
			{
				int index1D = getIndex1D(i, j, k, nMeshNodes);
				int boundaryAdjacentIndex{0}, nextToAdjacentIndex{0};
				getAdjacentIndices(index1D, nMeshNodes, boundaryAdjacentIndex, nextToAdjacentIndex);
				double rho = flowVariables.conservedVariables.rho(index1D);
				double T = 0;
				double p = (T*(1+rho) + rho) / params.Gamma;
				double uTripFlow=u;
				if(params.Re>50 && j<nMeshNodes.j/2 && t<50)
					uTripFlow *=0.9;
				else if(params.Re>50 && j<nMeshNodes.j/2 && t<60)
					uTripFlow = u - (60.-t)/10.*0.1*u;
				PrimitiveVariablesScalars primitiveVarsScalars(uTripFlow, v, w, p, T);
				ConservedVariablesScalars   conservedVarsScalars = deriveConservedVariables (primitiveVarsScalars, params);
				TransportPropertiesScalars transportPropsScalars = deriveTransportProperties(primitiveVarsScalars, params);
				setFlowVariablesAtNode(index1D, conservedVarsScalars, primitiveVarsScalars, transportPropsScalars, flowVariables);
			}
}

// Constructor, taking no extra arguments
OutletBoundary::OutletBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex)
: MeshEdgeBoundary(normalAxis, planeIndex, NodeTypeEnum::FluidEdge)
{}

// Set pressure fluctuation to zero. Linearly extrapolate density and momentum. Derive rest.
void OutletBoundary::applyBoundaryCondition(double t, const Vector3_i& nMeshNodes,		// <- Input
		   	   	   	   	   	   	   	   	   const ConfigSettings& params,				// <- Input
										   AllFlowVariablesArrayGroup& flowVariables)	// <- Output
{
	for(int i=ownedNodes.iMin; i<=ownedNodes.iMax; ++i)
		for(int j=ownedNodes.jMin; j<=ownedNodes.jMax; ++j)
			for(int k=ownedNodes.kMin; k<=ownedNodes.kMax; ++k)
			{
				int index1D = getIndex1D(i, j, k, nMeshNodes);
				int boundaryAdjacentIndex{0}, nextToAdjacentIndex{0};
				getAdjacentIndices(index1D, nMeshNodes, boundaryAdjacentIndex, nextToAdjacentIndex);
				double rho = 2*flowVariables.conservedVariables.rho(boundaryAdjacentIndex) // Linear extrapolation
							 - flowVariables.conservedVariables.rho(nextToAdjacentIndex);
				double rho_u = 2*flowVariables.conservedVariables.rho_u(boundaryAdjacentIndex) // Linear extrapolation
							   - flowVariables.conservedVariables.rho_u(nextToAdjacentIndex);
				double rho_v = 2*flowVariables.conservedVariables.rho_v(boundaryAdjacentIndex) // Linear extrapolation
							   - flowVariables.conservedVariables.rho_v(nextToAdjacentIndex);
				double rho_w = 2*flowVariables.conservedVariables.rho_w(boundaryAdjacentIndex) // Linear extrapolation
							   - flowVariables.conservedVariables.rho_w(nextToAdjacentIndex);
				double p = 0;

				double u = rho_u / (1+rho);
				double v = rho_v / (1+rho);
				double w = rho_w / (1+rho);
				double T = ( params.Gamma * p - rho ) / ( 1 + rho );
				double rho_E = p / ( params.Gamma - 1 ) + (1 + rho)/2 * ( u*u + v*v + w*w );

				PrimitiveVariablesScalars primitiveVarsScalars(u, v, w, p, T);
				ConservedVariablesScalars conservedVarsScalars(rho, rho_u, rho_v, rho_w, rho_E);
				TransportPropertiesScalars transportPropsScalars = deriveTransportProperties(primitiveVarsScalars, params);
				setFlowVariablesAtNode(index1D, conservedVarsScalars, primitiveVarsScalars, transportPropsScalars, flowVariables);
			}
}

// Constructor, taking no extra arguments
PeriodicBoundary::PeriodicBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex)
: MeshEdgeBoundary(normalAxis, planeIndex, NodeTypeEnum::FluidGhost)
{}

// Set all flow variables to the value of the boundary-adjacent node at the opposite side of the mesh.
void PeriodicBoundary::applyBoundaryCondition(double t, const Vector3_i& nMeshNodes,		// <- Input
	   	   	   	   	   	   	   	   	   	   	  const ConfigSettings& params,					// <- Input
											  AllFlowVariablesArrayGroup& flowVariables)	// <- Output
{
	for(int i=ownedNodes.iMin; i<=ownedNodes.iMax; ++i)
		for(int j=ownedNodes.jMin; j<=ownedNodes.jMax; ++j)
			for(int k=ownedNodes.kMin; k<=ownedNodes.kMax; ++k)
			{
				int index1D = getIndex1D(i, j, k, nMeshNodes);
				int oppositeSideIndex = getPeriodicIndex(index1D, nMeshNodes);
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

// Constructor, taking no extra arguments
SymmetryBoundary::SymmetryBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex)
: MeshEdgeBoundary(normalAxis, planeIndex, NodeTypeEnum::FluidGhost)
{}

// Symmetry plane is the adjacent interior nodes. The owned boundary nodes are ghost nodes.
// In these ghost nodes all primitive variables are mirrored from the second adjacent interior nodes,
// (mirrored across the symmetry plane), except the normal velocity component, which changes sign.
void SymmetryBoundary::applyBoundaryCondition(double t, const Vector3_i& nMeshNodes,		// <- Input
	   	   	   	   	   	  	  	  	  	  	  const ConfigSettings& params,					// <- Input
											  AllFlowVariablesArrayGroup& flowVariables)	// <- Output
{
	for(int i=ownedNodes.iMin; i<=ownedNodes.iMax; ++i)
		for(int j=ownedNodes.jMin; j<=ownedNodes.jMax; ++j)
			for(int k=ownedNodes.kMin; k<=ownedNodes.kMax; ++k)
			{
				int index1D = getIndex1D(i, j, k, nMeshNodes);
				int boundaryAdjacentIndex{0}, nextToAdjacentIndex{0};
				getAdjacentIndices(index1D, nMeshNodes, boundaryAdjacentIndex, nextToAdjacentIndex);
				double u, v, w;
				if(normalAxis == AxisOrientationEnum::x)
				{	// Normal component zero, and other components get zero gradient.
					u = -flowVariables.primitiveVariables.u(nextToAdjacentIndex);
					v =  flowVariables.primitiveVariables.v(nextToAdjacentIndex);
					w =  flowVariables.primitiveVariables.w(nextToAdjacentIndex);
				}
				else if(normalAxis == AxisOrientationEnum::y)
				{
					u =  flowVariables.primitiveVariables.u(nextToAdjacentIndex);
					v = -flowVariables.primitiveVariables.v(nextToAdjacentIndex);
					w =  flowVariables.primitiveVariables.w(nextToAdjacentIndex);
				}
				else if(normalAxis == AxisOrientationEnum::z)
				{
					u =  flowVariables.primitiveVariables.u(nextToAdjacentIndex);
					v =  flowVariables.primitiveVariables.v(nextToAdjacentIndex);
					w = -flowVariables.primitiveVariables.w(nextToAdjacentIndex);
				}
				else
					throw std::logic_error("Unexpected enum value");

				double p = flowVariables.primitiveVariables.p(nextToAdjacentIndex);
				double T = flowVariables.primitiveVariables.T(nextToAdjacentIndex);
				PrimitiveVariablesScalars  primitiveVarsScalars(u, v, w, p, T);
				ConservedVariablesScalars conservedVarsScalars   = deriveConservedVariables (primitiveVarsScalars, params);
				TransportPropertiesScalars transportPropsScalars = deriveTransportProperties(primitiveVarsScalars, params);
				setFlowVariablesAtNode(index1D, conservedVarsScalars, primitiveVarsScalars, transportPropsScalars, flowVariables);
			}
}

// Constructor, taking no extra arguments
ExtrapolationBoundary::ExtrapolationBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex)
: MeshEdgeBoundary(normalAxis, planeIndex, NodeTypeEnum::FluidGhost)
{}

// Everything is extrapolated
void ExtrapolationBoundary::applyBoundaryCondition(double t, const Vector3_i& nMeshNodes,		// <- Input
												   const ConfigSettings& params,				// <- Input
												   AllFlowVariablesArrayGroup& flowVariables)	// <- Output
{
	for(int i=ownedNodes.iMin; i<=ownedNodes.iMax; ++i)
		for(int j=ownedNodes.jMin; j<=ownedNodes.jMax; ++j)
			for(int k=ownedNodes.kMin; k<=ownedNodes.kMax; ++k)
			{
				int index1D = getIndex1D(i, j, k, nMeshNodes);
				int boundaryAdjacentIndex{0}, nextToAdjacentIndex{0};
				getAdjacentIndices(index1D, nMeshNodes, boundaryAdjacentIndex, nextToAdjacentIndex);
				double rho = 2*flowVariables.conservedVariables.rho(boundaryAdjacentIndex) // Linear extrapolation
							 - flowVariables.conservedVariables.rho(nextToAdjacentIndex);
				double rho_u = 2*flowVariables.conservedVariables.rho_u(boundaryAdjacentIndex) // Linear extrapolation
							   - flowVariables.conservedVariables.rho_u(nextToAdjacentIndex);
				double rho_v = 2*flowVariables.conservedVariables.rho_v(boundaryAdjacentIndex) // Linear extrapolation
							   - flowVariables.conservedVariables.rho_v(nextToAdjacentIndex);
				double rho_w = 2*flowVariables.conservedVariables.rho_w(boundaryAdjacentIndex) // Linear extrapolation
							   - flowVariables.conservedVariables.rho_w(nextToAdjacentIndex);
				double rho_E = 2*flowVariables.conservedVariables.rho_E(boundaryAdjacentIndex) // Linear extrapolation
							   - flowVariables.conservedVariables.rho_E(nextToAdjacentIndex);
				ConservedVariablesScalars conservedVarsScalars(rho, rho_u, rho_v, rho_w, rho_E);
				PrimitiveVariablesScalars primitiveVarsScalars = derivePrimitiveVariables(conservedVarsScalars, params);
				TransportPropertiesScalars transportPropsScalars = deriveTransportProperties(primitiveVarsScalars, params);
				setFlowVariablesAtNode(index1D, conservedVarsScalars, primitiveVarsScalars, transportPropsScalars, flowVariables);
			}
}

// Constructor, taking in info about neighbor region.
SubmeshInterfaceBoundary::SubmeshInterfaceBoundary(AxisOrientationEnum normalAxis,
											 	   EdgeIndexEnum planeIndex,
												   AllFlowVariablesArrayGroup& neighborSubMesh)
: MeshEdgeBoundary(normalAxis, planeIndex, NodeTypeEnum::FluidGhost),
  neighborSubMesh(neighborSubMesh)
{}

// Borrow nodes from neighbor region as BC for this region.
void SubmeshInterfaceBoundary::applyBoundaryCondition(double t,
													  const Vector3_i& nMeshNodes,
													  const ConfigSettings& params,
													  AllFlowVariablesArrayGroup& flowVariables)
{

}

// Constructor. Use derived classes instead.
ImmersedBoundary::ImmersedBoundary()
{

}

// Get interpolation values, to compute values in image points, and then the ghost nodes.
// Also filters density and energy in the closest fluid nodes.
void ImmersedBoundary::applyBoundaryCondition(const ConfigSettings& params,
											  MeshDescriptor& mesh // <- In/Output
  	  	    								  )
{
	for(GhostNode ghostNode : ghostNodes)
	{
		InterpolationValues interpolationValues;		// Conserved variables at the 8 surrounding interpolation points.
		InterpolationPositions interpolationPositions;	// Positions of interpolation points
		Array8_b ghostFlag;					// A bool flag for each surrounding node. True if ghost.
		ghostFlag.fill(false);				// <- Sets all 8 values to false
		bool allSurroundingAreFluid{true};	// Until otherwise proven
		vector<Vector3_d> unitNormals;		// For each of the surrounding nodes that is ghost, we put its unit normal probe here.
		IndexBoundingBox surroundingNodes = getSurroundingNodesBox(ghostNode.imagePoint, mesh.spacings, mesh.originOffset, mesh.nNodes);
		
		setInterpolationValues(
				surroundingNodes, mesh,								// <- Input
				interpolationValues, interpolationPositions,		// <- Output
				ghostFlag, allSurroundingAreFluid, unitNormals);	// <- Output
		
		PrimitiveVariablesScalars imagePointPrimVars = interpolateImagePointVariables(
				interpolationValues, interpolationPositions, allSurroundingAreFluid,
				ghostFlag, ghostNode, surroundingNodes, unitNormals, mesh);

		PrimitiveVariablesScalars ghostNodePrimVars = getGhostNodeBCVariables(imagePointPrimVars);
		ConservedVariablesScalars ghostNodeConsVars = deriveConservedVariables(ghostNodePrimVars, params);
		TransportPropertiesScalars ghostNodeTransportProps = deriveTransportProperties(ghostNodePrimVars, params);
		setFlowVariablesAtNode(ghostNode.indices, ghostNodeConsVars, ghostNodePrimVars, ghostNodeTransportProps,
								mesh.flowVariables);
	}
}

// Given a vector of indices to solid nodes, find all who will be part of the numerical stencil
// of any fluid node. Those solid nodes are marked as solidGhost, and added to the ghostNodes vector.
void ImmersedBoundary::findGhostNodesWithFluidNeighbors(const vector<int>& solidNodeIndices,
														const Vector3_i& nMeshNodes,
														Array3D_nodeType& nodeTypeArray)
{
	IndexBoundingBox meshSize(nMeshNodes.i-1, nMeshNodes.j-1, nMeshNodes.k-1);
	for (int index1D : solidNodeIndices)
	{
		// In 3D, each node has 26 neighbors. All except the 8 corner nodes (i+1,j+1,k+1), (i+1,j+1,k-1), etc,
		// are part of the 2nd order stencil for the N-S equations. Thus we must check 18 nodes:
		Vector3_i solidNode = getIndices3D(index1D, nMeshNodes);
		bool solidNodeHasFluidNeighbor { false };

		if (solidNode.i > meshSize.iMin) // Then we can test all nodes at i-1:
		{
			if (nodeTypeArray(solidNode.i - 1, solidNode.j, solidNode.k) == NodeTypeEnum::FluidInterior)
				solidNodeHasFluidNeighbor = true; // (-1,0,0)

			if (solidNode.j > meshSize.jMin)
				if (nodeTypeArray(solidNode.i - 1, solidNode.j - 1, solidNode.k) == NodeTypeEnum::FluidInterior)
					solidNodeHasFluidNeighbor = true; // (-1,-1,0)

			if (solidNode.j < meshSize.jMax)
				if (nodeTypeArray(solidNode.i - 1, solidNode.j + 1, solidNode.k) == NodeTypeEnum::FluidInterior)
					solidNodeHasFluidNeighbor = true; // (-1,1,0)

			if (solidNode.k > meshSize.kMin)
				if (nodeTypeArray(solidNode.i - 1, solidNode.j, solidNode.k - 1) == NodeTypeEnum::FluidInterior)
					solidNodeHasFluidNeighbor = true; // (-1,0,-1)

			if (solidNode.k < meshSize.kMax)
				if (nodeTypeArray(solidNode.i - 1, solidNode.j, solidNode.k + 1) == NodeTypeEnum::FluidInterior)
					solidNodeHasFluidNeighbor = true; // (-1,0,1)
		}

		if (solidNode.i < meshSize.iMax) // Then we can test all nodes at i+1:
		{
			if (nodeTypeArray(solidNode.i + 1, solidNode.j, solidNode.k) == NodeTypeEnum::FluidInterior)
				solidNodeHasFluidNeighbor = true; // (1,0,0)

			if (solidNode.j > meshSize.jMin)
				if (nodeTypeArray(solidNode.i + 1, solidNode.j - 1, solidNode.k) == NodeTypeEnum::FluidInterior)
					solidNodeHasFluidNeighbor = true; // (1,-1,0)

			if (solidNode.j < meshSize.jMax)
				if (nodeTypeArray(solidNode.i + 1, solidNode.j + 1, solidNode.k) == NodeTypeEnum::FluidInterior)
					solidNodeHasFluidNeighbor = true; // (1,1,0)

			if (solidNode.k > meshSize.kMin)
				if (nodeTypeArray(solidNode.i + 1, solidNode.j, solidNode.k - 1) == NodeTypeEnum::FluidInterior)
					solidNodeHasFluidNeighbor = true; // (1,0,-1)

			if (solidNode.k < meshSize.kMax)
				if (nodeTypeArray(solidNode.i + 1, solidNode.j, solidNode.k + 1) == NodeTypeEnum::FluidInterior)
					solidNodeHasFluidNeighbor = true; // (1,0,1)
		}

		// And then we check the nodes at the same i-index (x-coordinate) as the solid node:

		if (solidNode.j > meshSize.jMin)
		{
			if (nodeTypeArray(solidNode.i, solidNode.j - 1, solidNode.k) == NodeTypeEnum::FluidInterior)
				solidNodeHasFluidNeighbor = true; // (0,-1,0

			if (solidNode.k > meshSize.kMin)
				if (nodeTypeArray(solidNode.i, solidNode.j - 1, solidNode.k - 1) == NodeTypeEnum::FluidInterior)
					solidNodeHasFluidNeighbor = true; // (0,-1,-1)

			if (solidNode.k < meshSize.kMax)
				if (nodeTypeArray(solidNode.i, solidNode.j - 1, solidNode.k + 1) == NodeTypeEnum::FluidInterior)
					solidNodeHasFluidNeighbor = true; // (0,-1,1)
		}

		if (solidNode.j < meshSize.jMax)
		{
			if (nodeTypeArray(solidNode.i, solidNode.j + 1, solidNode.k) == NodeTypeEnum::FluidInterior)
				solidNodeHasFluidNeighbor = true; // (0,1,0)

			if (solidNode.k > meshSize.kMin)
				if (nodeTypeArray(solidNode.i, solidNode.j + 1, solidNode.k - 1) == NodeTypeEnum::FluidInterior)
					solidNodeHasFluidNeighbor = true; // (0,1,-1)

			if (solidNode.k < meshSize.kMax)
				if (nodeTypeArray(solidNode.i, solidNode.j + 1, solidNode.k + 1) == NodeTypeEnum::FluidInterior)
					solidNodeHasFluidNeighbor = true; // (0,1,1)
		}

		if (solidNode.k > meshSize.kMin)
			if (nodeTypeArray(solidNode.i, solidNode.j, solidNode.k - 1) == NodeTypeEnum::FluidInterior)
				solidNodeHasFluidNeighbor = true; // (0,0,-1)

		if (solidNode.k < meshSize.kMax)
			if (nodeTypeArray(solidNode.i, solidNode.j, solidNode.k + 1) == NodeTypeEnum::FluidInterior)
				solidNodeHasFluidNeighbor = true; // (0,0,1)

		if (solidNodeHasFluidNeighbor)
		{
			ghostNodeMap[index1D] = ghostNodes.size();
			ghostNodes.emplace_back(solidNode);
			nodeTypeArray(solidNode) = NodeTypeEnum::SolidGhost;
		}
	}
}

// Given one of the 8 nodes surrounding an image point, make sure it is not inactive.
// Mostly the 'findGhostNodesWithFluidNeighbors' function will find all solids that should
// be ghost, but when the immersed boundary intersects a mesh edge boundary (and possibly
// other cases), an image point can be surrounded by a solid that's not part of a stencil.
void ImmersedBoundary::checkIfSurroundingShouldBeGhost(const Vector3_i &surroundingNode,
													   vector<GhostNode>& newGhostNodes,
													   Array3D_nodeType& nodeTypeArray)
{
	if (nodeTypeArray(surroundingNode) == NodeTypeEnum::SolidInactive)
	{
		newGhostNodes.emplace_back(surroundingNode);
		nodeTypeArray(surroundingNode) = NodeTypeEnum::SolidGhost;
	}
}

// Given an iterator to a GhostNode in a vector, go through the ghosts after the iterator, and find
// their image points (IP). Also check if any surrounding node of the IP are solidInactive. If so,
// mark them as ghost ad add them to the return-vector.
vector<GhostNode> ImmersedBoundary::setImagePointPositions(GhostNodeVectorIterator firstGhostToProcess,
														   const Vector3_d& gridSpacing,
														   const Vector3_i& nMeshNodes,
														   const Vector3_d& meshOriginOffset,
														   Array3D_nodeType& nodeTypeArray)
{
	vector<GhostNode> newGhostNodes;
	for (GhostNodeVectorIterator ghostIterator{firstGhostToProcess}; ghostIterator!=ghostNodes.end(); ++ghostIterator)
	{
		GhostNode& ghostNode = *ghostIterator;
		Vector3_d ghostNodePosition = getNodePosition(ghostNode.indices, gridSpacing, meshOriginOffset);
		Vector3_d normalProbe = getNormalProbe(ghostNodePosition); // from ghost to body intercept point
		ghostNode.bodyInterceptPoint = ghostNodePosition + normalProbe;
		ghostNode.imagePoint = ghostNode.bodyInterceptPoint + normalProbe;
		IndexBoundingBox surroundingNodes = getSurroundingNodesBox(ghostNode.imagePoint, gridSpacing, meshOriginOffset, nMeshNodes);
		for( int surroundingNodeIndex1D : surroundingNodes.asIndexList(nMeshNodes) )
			checkIfSurroundingShouldBeGhost( getIndices3D(surroundingNodeIndex1D, nMeshNodes),
											 newGhostNodes, nodeTypeArray );
		// Because of the way that lift/drag etc is computed we need also BI to be only surrounded by ghost or fluid:
		surroundingNodes = getSurroundingNodesBox(ghostNode.bodyInterceptPoint, gridSpacing, meshOriginOffset, nMeshNodes);
		for( int surroundingNodeIndex1D : surroundingNodes.asIndexList(nMeshNodes) )
			checkIfSurroundingShouldBeGhost( getIndices3D(surroundingNodeIndex1D, nMeshNodes),
											 newGhostNodes, nodeTypeArray );
	}
	return newGhostNodes;
}

// Add the nodes from newGhostNodes at the end of ghostNodes. Return iterator to the first appended node.
GhostNodeVectorIterator ImmersedBoundary::appendGhostNodes(const vector<GhostNode>& newGhostNodes,
														   const Vector3_i& nMeshNodes)
{
	for(const GhostNode& newGhostNode : newGhostNodes)
	{
		int index1D = getIndex1D( newGhostNode.indices, nMeshNodes );
		ghostNodeMap[index1D] = ghostNodes.size();
		ghostNodes.push_back(newGhostNode);
	}
	return ghostNodes.end() - newGhostNodes.size();
}

// Filter density and energy in the fluid nodes closest to the immersed surface.
void ImmersedBoundary::filterClosestFluidNodes(const Vector3_i& nMeshNodes,
							 	 	 	 	   const AllFlowVariablesArrayGroup& flowVariables)
{
	vector<double> rhoFilteredVector;
	vector<double> energyFilteredVector;
	for(int index1D : filterNodes)
	{
		Vector3_i indices = getIndices3D(index1D, nMeshNodes);
		double rhoFiltered = ( 6*flowVariables.conservedVariables.rho(indices)
							   + flowVariables.conservedVariables.rho(indices.i+1, indices.j, indices.k)
							   + flowVariables.conservedVariables.rho(indices.i-1, indices.j, indices.k)
							   + flowVariables.conservedVariables.rho(indices.i, indices.j+1, indices.k)
							   + flowVariables.conservedVariables.rho(indices.i, indices.j-1, indices.k)
							   + flowVariables.conservedVariables.rho(indices.i, indices.j, indices.k+1)
							   + flowVariables.conservedVariables.rho(indices.i, indices.j, indices.k-1)
							   ) / 12;
		double energyFiltered = ( 6*flowVariables.conservedVariables.rho_E(indices)
							      + flowVariables.conservedVariables.rho_E(indices.i+1, indices.j, indices.k)
							      + flowVariables.conservedVariables.rho_E(indices.i-1, indices.j, indices.k)
							      + flowVariables.conservedVariables.rho_E(indices.i, indices.j+1, indices.k)
							      + flowVariables.conservedVariables.rho_E(indices.i, indices.j-1, indices.k)
							      + flowVariables.conservedVariables.rho_E(indices.i, indices.j, indices.k+1)
							      + flowVariables.conservedVariables.rho_E(indices.i, indices.j, indices.k-1)
							      ) / 12;
		rhoFilteredVector	.push_back(rhoFiltered);
		energyFilteredVector.push_back(energyFiltered);
	}
	int counter = 0;
	for(int index1D : filterNodes)
	{
		flowVariables.conservedVariables.rho  (index1D) = rhoFilteredVector   [counter];
		flowVariables.conservedVariables.rho_E(index1D) = energyFilteredVector[counter];
		++counter;
	}
}

// Set interpolation values for one of the 8 nodes surrounding an image point (IP). Surrounding node must be fluid.
void ImmersedBoundary::setInterpolationValuesFluidNode(int counter,	// <- Which of the 8 surrounding nodes [0-7]
													   int surroundingNodeIndex1D,
													   const MeshDescriptor& mesh,
													   InterpolationValues& interpolationValues,		// <- OUTPUT
													   InterpolationPositions& interpolationPositions)	// <- OUTPUT
{
	interpolationValues.u[counter] = mesh.flowVariables.primitiveVariables.u(surroundingNodeIndex1D);
	interpolationValues.v[counter] = mesh.flowVariables.primitiveVariables.v(surroundingNodeIndex1D);
	interpolationValues.w[counter] = mesh.flowVariables.primitiveVariables.w(surroundingNodeIndex1D);
	interpolationValues.p[counter] = mesh.flowVariables.primitiveVariables.p(surroundingNodeIndex1D);
	interpolationValues.T[counter] = mesh.flowVariables.primitiveVariables.T(surroundingNodeIndex1D);
	Vector3_d surroundingNodePosition = getNodePosition( getIndices3D(surroundingNodeIndex1D, mesh.nNodes),
															mesh.spacings, mesh.originOffset );
	interpolationPositions.x[counter] = surroundingNodePosition.x;
	interpolationPositions.y[counter] = surroundingNodePosition.y;
	interpolationPositions.z[counter] = surroundingNodePosition.z;
}

// Set interpolation values for one of the 8 nodes surrounding an image point (IP). Surrounding node must be ghost.
void ImmersedBoundary::setInterpolationValuesGhostNode(
		int counter,	// <- Which of the 8 surrounding nodes [0-7]
		int surroundingNodeIndex1D,
		vector<Vector3_d>& unitNormals,					// <-
		InterpolationValues& interpolationValues,		// <- OUTPUT
		InterpolationPositions& interpolationPositions)	// <-
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
		const MeshDescriptor& mesh,
		InterpolationValues& interpolationValues,		// <-
		InterpolationPositions& interpolationPositions,	// <-
		Array8_b& ghostFlag,							// <- Output
		bool& allSurroundingAreFluid,					// <-
		vector<Vector3_d>& unitNormals)					// <-
{
	int counter { 0 }; // 0, 1, ..., 7.  To index the interpolation point arrays.
	for (int surroundingNodeIndex1D : surroundingNodes.asIndexList(mesh.nNodes))
	{
		if (mesh.nodeType(surroundingNodeIndex1D) == NodeTypeEnum::FluidInterior
		||	mesh.nodeType(surroundingNodeIndex1D) == NodeTypeEnum::FluidEdge
		||	mesh.nodeType(surroundingNodeIndex1D) == NodeTypeEnum::FluidGhost)
		{
			setInterpolationValuesFluidNode(counter, surroundingNodeIndex1D, mesh,	// <- Input
											interpolationValues, interpolationPositions);	// <- Output
		}
		else if (mesh.nodeType(surroundingNodeIndex1D) == NodeTypeEnum::SolidGhost)
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

// Get image point value of a flow variable, using the fast simplified interpolation.
// Only valid if all 8 surrounding nodes are fluid.
double ImmersedBoundary::simplifiedInterpolation(const Vector8_d& interpolationValues,
												 const Vector3_i& lowerIndexNode, // <- Node with smallest coordinates
												 const Vector3_d& imagePointPosition,
												 const Vector3_d& gridSpacing,
												 const Vector3_d& meshOriginOffset)
{
	Vector3_d lowerIndexNodePosition = getNodePosition(lowerIndexNode, gridSpacing, meshOriginOffset);
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
		const vector<Vector3_d>& unitNormals,
		const MeshDescriptor& mesh )
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
CylinderBody::CylinderBody(Vector3_d centroidPosition, AxisOrientationEnum axis, double radius) :
centroidPosition(centroidPosition),
axis{axis},
radius{radius}
{}

// Mark solid nodes (ghost and inactive), find image points (IP) and closest fluid nodes (to filter).
void CylinderBody::identifyRelatedNodes(const ConfigSettings& params,
		   	   	   	   	   	   	   	    const Vector3_d& gridSpacing,
										const Vector3_i& nMeshNodes,
										const Vector3_d& meshOriginOffset,
										Array3D_nodeType& nodeTypeArray	// <- Output
										)
{
	int filterNodesLayerWidth = 2; // TODO: param
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
IntegralProperties CylinderBody::getIntegralProperties(const ConfigSettings& params, const MeshDescriptor& mesh)
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
													vector<int>& solidNodeIndices, // <- Output
													Array3D_nodeType& nodeTypeArray   // <- Output
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
SphereBody::SphereBody(Vector3_d centerPosition, double radius) :
centerPosition(centerPosition),
radius{radius}
{}

// Mark solid nodes (ghost and inactive), find image points (IP) and closest fluid nodes (to filter).
void SphereBody::identifyRelatedNodes(const ConfigSettings& params,
	   	   	   	    				  const Vector3_d& gridSpacing,
									  const Vector3_i& nMeshNodes,
									  const Vector3_d& meshOriginOffset,
									  Array3D_nodeType& nodeTypeArray	// <- Output
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
IntegralProperties SphereBody::getIntegralProperties(const ConfigSettings& params, const MeshDescriptor& mesh)
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
												vector<int>& solidNodeIndices, // <- Output
												Array3D_nodeType& nodeTypeArray   // <- Output
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





