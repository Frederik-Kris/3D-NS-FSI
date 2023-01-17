/*
 * Mesh.cpp
 *
 *  Created on: Oct 7, 2022
 *      Author: frederk
 */

#include "Mesh.h"

Mesh::Mesh(const ConfigSettings& params) :
NI{params.NI}, NJ{params.NJ}, NK{params.NK},
nNodesTotal{NI*NJ*NK},
dx{ -1 }, dy{ -1 }, dz{ -1 }, // <- Depend on boundary condition types, must set later
conservedVariables(NI, NJ, NK),
conservedVariablesOld(NI, NJ, NK),
primitiveVariables(NI, NJ, NK),
transportProperties(NI, NJ, NK),
RK4slopes(NI, NJ, NK),
nodeType(NI, NJ, NK)
{
}


// Construct the objects that define the boundary conditions (BCs) and set the grid spacings based on that.
void Mesh::setupBoundaries(const ConfigSettings &params)
{
	// First define the BCs at the planes bounding the Cartesian mesh.
	// The order here matters, in the sense that the first added BC will claim all the nodes in its plane.
	// The following BCs can only claim the unclaimed nodes in their planes.
	// E.g., if the first added BC is at z=z_min (enum values z and min), then nodes i=[0,imax], j=[0,jmax], k=kmin,
	// are claimed and cannot be part of another BC.
	// The order which the BCs are applied need to be the reverse of this, to be correct.
	edgeBoundaries.push_back(std::make_unique<PeriodicBoundary>(AxisOrientationEnum::z, EdgeIndexEnum::min));
	edgeBoundaries.push_back(std::make_unique<PeriodicBoundary>(AxisOrientationEnum::z, EdgeIndexEnum::max));
	edgeBoundaries.push_back(std::make_unique<SymmetryBoundary>(AxisOrientationEnum::y, EdgeIndexEnum::min));
	edgeBoundaries.push_back(std::make_unique<SymmetryBoundary>(AxisOrientationEnum::y, EdgeIndexEnum::max));
	edgeBoundaries.push_back(std::make_unique<InletBoundary>   (AxisOrientationEnum::x, EdgeIndexEnum::min, params.M_0));
	edgeBoundaries.push_back(std::make_unique<OutletBoundary>  (AxisOrientationEnum::x, EdgeIndexEnum::max));

	// Set grid spacing and origin point offset based on the BC types:
	// In/Out:		spacing=N-1, offset=0
	// Periodic:	spacing=N-2, offset=1 or 0
	// Symmetry:	spacing=N-3, offset=1
	dx = params.L_x / (NI-1);
	dy = params.L_y / (NJ-3);
	dz = params.L_z / (NK-2);
	positionOffset.x = 0;
	positionOffset.y = 1;
	positionOffset.z = 1;

	// Add immersed bodies:
	Vector3_d cylinderCentroidPosition(params.L_x / 4, params.L_y / 2, 0);
	immersedBoundaries.push_back(std::make_unique<CylinderBody>(cylinderCentroidPosition,
																AxisOrientationEnum::z,
																params.L_y/10));
//	Vector3_d sphereCenterPoint(params.L_x/4, params.L_y/2, params.L_z/2);
//	immersedBoundaries.push_back(std::make_unique<SphereBody>(sphereCenterPoint, params.L_y/8.1));
}

// Categorize mesh nodes based on boundary conditions. This includes finding image point positions.
void Mesh::categorizeNodes(const ConfigSettings& params)
{
	const Vector3_u nMeshNodes(NI, NJ, NK);
	const Vector3_d gridSpacing(dx, dy, dz);
	nodeType.setAll(NodeTypeEnum::FluidInterior);

	IndexBoundingBox unclaimedNodes(NI-1, NJ-1, NK-1);						// At first, all nodes are unclaimed.
	for(auto&& boundary : edgeBoundaries)									// Loop through mesh edge boundaries.
		boundary->identifyOwnedNodes(unclaimedNodes, nMeshNodes, nodeType);	// Each boundary flags and claims nodes.
	std::reverse( edgeBoundaries.begin(), edgeBoundaries.end() );			// Reverse order, so last added is first applied.

	// Then each immersed boundary finds its solid and ghost nodes, and body intercept and image points:
	for(auto&& boundary : immersedBoundaries)
		boundary->identifyRelatedNodes(params, gridSpacing, nMeshNodes, positionOffset, nodeType);

	// Finally we note the indices of interior fluid nodes and solid ghost nodes, so we can loop through only those:
	for(size_t index1D{0}; index1D<nNodesTotal; ++index1D)
		if(nodeType(index1D) == NodeTypeEnum::FluidInterior)
			indexByType.fluidInterior.push_back(index1D);
		else if(nodeType(index1D) == NodeTypeEnum::SolidGhost)
			indexByType.solidGhost.push_back(index1D);
}

// Check if the 2nd order filter should be applied at the current time
bool Mesh::checkFilterCondition(const ConfigSettings& params, ulong timeLevel, double t)
{
	// Loop through the list of specified times in the config file, and see if we have surpassed them:
	size_t counter = 0;
	for (double intervalChangeTime : params.filterIntervalChangeTimes)
		if (t > intervalChangeTime)
			++counter;
	bool filterNow = false;
	// If we passed fewer thresholds than there are entries in filterIntervals, then select the correct
	// filter interval and check if timeLevel+1 is divisible by that interval:
	if ( counter < params.filterIntervals.size() )
	{
		uint currentFilterInterval = params.filterIntervals.at(counter);
		if(currentFilterInterval > 0)
			if( (timeLevel+1) % currentFilterInterval == 0 )
				filterNow = true;
	}
	return filterNow;
}

// Filters one variable field, i.e., one solution array, 'filterVariable' and stores the filtered
// result in 'variableTemporaryStorage'. Then, the arrays are swapped by move-semantics.
// Only filters if conditions set in config file are met.
void Mesh::applyFilter_ifAppropriate(Array3D_d& filterVariable, Array3D_d& variableTemporaryStorage,
									const ConfigSettings& params, ulong timeLevel, double t)
{
	if( checkFilterCondition(params, timeLevel, t) )
	{
		// Copy boundary nodes:
		for (auto&& edgeBoundary : edgeBoundaries)
		{
			const IndexBoundingBox& ownedNodes = edgeBoundary->ownedNodes;
			for(size_t i=ownedNodes.iMin; i<=ownedNodes.iMax; ++i)
				for(size_t j=ownedNodes.jMin; j<=ownedNodes.jMax; ++j)
					for(size_t k=ownedNodes.kMin; k<=ownedNodes.kMax; ++k)
						variableTemporaryStorage(i,j,k) = filterVariable(i,j,k);
		}
		for (size_t index1D : indexByType.solidGhost)
			variableTemporaryStorage(index1D) = filterVariable(index1D);
		// Apply filter to interior nodes:
		Vector3_u nNodes(NI, NJ, NK);
		for (size_t index1D : indexByType.fluidInterior)
		{
			Vector3_u indices = getIndices3D(index1D, nNodes);
			variableTemporaryStorage(index1D) = (6*filterVariable(indices)
												 + filterVariable(indices.i+1, indices.j, indices.k)
												 + filterVariable(indices.i-1, indices.j, indices.k)
												 + filterVariable(indices.i, indices.j+1, indices.k)
												 + filterVariable(indices.i, indices.j-1, indices.k)
												 + filterVariable(indices.i, indices.j, indices.k+1)
												 + filterVariable(indices.i, indices.j, indices.k-1)
												 ) / 12;
		}
		filterVariable.dataSwap(variableTemporaryStorage);	// Swap the arrays using move-sematics (super-fast)
	}
}

// Swap the contents of all the arrays of conserved variables and the intermediate arrays, by move-semantics.
// This operation is super fast and needs no extra copy. Only the ownership of the data is changed.
void Mesh::swapConservedVariables()
{
	conservedVariables.rho  .dataSwap(conservedVariablesOld.rho  );
	conservedVariables.rho_u.dataSwap(conservedVariablesOld.rho_u);
	conservedVariables.rho_v.dataSwap(conservedVariablesOld.rho_v);
	conservedVariables.rho_w.dataSwap(conservedVariablesOld.rho_w);
	conservedVariables.rho_E.dataSwap(conservedVariablesOld.rho_E);
}

// Loop through boundaries, who are derived classes of MeshEdgeBoundary and ImmersedBoundary,
// and apply their boundary conditions.
void Mesh::applyAllBoundaryConditions(double t, const ConfigSettings& params)
{
	const Vector3_u nMeshNodes(NI, NJ, NK);
	const Vector3_d gridSpacing(dx, dy, dz);
	AllFlowVariablesArrayGroup flowVariableReferences(conservedVariables, primitiveVariables, transportProperties);
	MeshDescriptor meshData(nMeshNodes, gridSpacing, positionOffset, nodeType, flowVariableReferences);

	for(auto&& boundary : edgeBoundaries)
		boundary->applyBoundaryCondition(t, nMeshNodes, params, flowVariableReferences);

	for(auto&& boundary : immersedBoundaries)
		boundary->applyBoundaryCondition(params, meshData);
}

// Compute the 2-norm of the difference between two arrys. Intended to monitor the change between two consecutive time levels,
// .g. to check convergence of solution. Could also be used to check error.
double Mesh::getNormOfChange(const Array3D_d& oldValue, const Array3D_d& newValue)
{
	double sumOfSquaredDifferences = 0;
	for(size_t i : indexByType.fluidInterior)
		sumOfSquaredDifferences += pow(oldValue(i)-newValue(i), 2);
	double normOfChange = sqrt( dx*dy*dz * sumOfSquaredDifferences );
	return normOfChange;
}





