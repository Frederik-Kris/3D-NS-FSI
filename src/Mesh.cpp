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
conservedVariables(NI, NJ, NK),
primitiveVariables(NI, NJ, NK),
transportProperties(NI, NJ, NK),
intermediateConservedVariables(NI, NJ, NK),
RK4slopes(NI, NJ, NK),
nodeTypes(NI, NJ, NK)
{
	setGridSpacings(params.L_x, params.L_y, params.L_z);
}


// Calculate the space between nodes in the grid, based on domain size and no. of nodes.
void Mesh::setGridSpacings(double domainLengthX,
						   double domainLengthY,
						   double domainLengthZ )
{
	dx = domainLengthX / (NI - 1);
	dy = domainLengthY / (NJ - 1);
	dz = domainLengthZ / (NK - 1);
	cout << "Grid spacings set: dx = " << dx << " , dy = " << dy << " , dz = " << dz << endl;
}

// Construct the objects that define the boundary conditions (BCs)
void Mesh::setupBoundaries(const ConfigSettings &params)
{
	edgeBoundaries.push_back(std::make_unique<InletBoundary>(AxisOrientationEnum::x, EdgeIndexEnum::min, params.M_0));
	edgeBoundaries.push_back(std::make_unique<OutletBoundary>(AxisOrientationEnum::x, EdgeIndexEnum::max));
	edgeBoundaries.push_back(std::make_unique<SymmetryBoundary>(AxisOrientationEnum::y, EdgeIndexEnum::min));
	edgeBoundaries.push_back(std::make_unique<SymmetryBoundary>(AxisOrientationEnum::y, EdgeIndexEnum::max));
	edgeBoundaries.push_back(std::make_unique<SymmetryBoundary>(AxisOrientationEnum::z, EdgeIndexEnum::min));
	edgeBoundaries.push_back(std::make_unique<SymmetryBoundary>(AxisOrientationEnum::z, EdgeIndexEnum::max));

	Vector3_d cylinderCentroidPosition(params.L_x / 5, params.L_y / 2, 0);
	immersedBoundaries.push_back(std::make_unique<CylinderBody>(cylinderCentroidPosition,
																AxisOrientationEnum::z,
																params.L_y/10));
}

// Sanity check for the combination of boundary conditions.
void Mesh::assertBoundaryConditionCompliance()
{
	// Implement when creating general IBM.
}

void Mesh::categorizeNodes(const ConfigSettings& params)
{
	const IndexBoundingBox meshSize(NI-1, NJ-1, NK-1);
	nodeTypes.setAll(NodeTypeEnum::FluidRegular);

	IndexBoundingBox unclaimedNodes = meshSize;
	for(auto&& boundary : edgeBoundaries)
		boundary->identifyOwnedNodes(unclaimedNodes, *this);

	for(auto&& boundary : immersedBoundaries)
		boundary->identifyRelatedNodes(params, *this);

	for(uint index1D{0}; index1D<nNodesTotal; ++index1D)
		if(nodeTypes(index1D) == NodeTypeEnum::FluidRegular)
			activeNodeIndices.push_back(index1D);
}

// Filters one variable field, i.e., one solution array, 'filterVariable' and stores the filtered
// result in 'variableTemporaryStorage'. Then, the arrays are swapped by move-semantics.
// Only filters if the modulo of time level plus one, by the filter interval is zero.
// This causes the first filtering to happen as late as possible.
void Mesh::applyFilter_ifAppropriate(Array3D_d& filterVariable, Array3D_d& variableTemporaryStorage,
									uint filterInterval, uint timeLevel)
{
	if(filterInterval > 0)
		if( (timeLevel+1) % filterInterval == 0 )
		{
			uint iMax{NI-1}, jMax{NJ-1}, kMax{NK-1};
			// Copy boundary nodes:
			for(uint i{0}; i<=iMax; ++i)
				for(uint j{0}; j<=jMax; ++j)
				{
					variableTemporaryStorage(i,j,0   ) = filterVariable(i,j,0   );
					variableTemporaryStorage(i,j,kMax) = filterVariable(i,j,kMax);
				}
			for(uint i{0}; i<=iMax; ++i)
				for(uint k{1}; k<=kMax-1; ++k)
				{
					variableTemporaryStorage(i,0,   k) = filterVariable(i,0,   k);
					variableTemporaryStorage(i,jMax,k) = filterVariable(i,jMax,k);
				}
			for(uint j{1}; j<=jMax-1; ++j)
				for(uint k{1}; k<=kMax-1; ++k)
				{
					variableTemporaryStorage(0,   j,k) = filterVariable(0,   j,k);
					variableTemporaryStorage(iMax,j,k) = filterVariable(iMax,j,k);
				}
			// Apply filter to interior nodes:
			for(uint i{1}; i<=iMax-1; ++i)
				for(uint j{1}; j<=jMax-1; ++j)
					for(uint k{1}; k<=kMax-1; ++k)
					{
						variableTemporaryStorage(i,j,k) = 1./2.  *   filterVariable(i,j,k)
														+ 1./12. * ( filterVariable(i+1,j,k) + filterVariable(i-1,j,k)
																   + filterVariable(i,j+1,k) + filterVariable(i,j-1,k)
																   + filterVariable(i,j,k+1) + filterVariable(i,j,k-1) );
					}
			filterVariable.dataSwap(variableTemporaryStorage);	// Swap the arrays using move-sematics (super-fast)
		}
}

// Compute the norm of change in the conserved variables, and store in the history vectors.
ConservedVariablesScalars Mesh::computeNorms_conservedVariables()
{
	double norm_rho  { getNormOfChange(conservedVariables.rho  , intermediateConservedVariables.rho  ) };
	double norm_rho_u{ getNormOfChange(conservedVariables.rho_u, intermediateConservedVariables.rho_u) };
	double norm_rho_v{ getNormOfChange(conservedVariables.rho_v, intermediateConservedVariables.rho_v) };
	double norm_rho_w{ getNormOfChange(conservedVariables.rho_w, intermediateConservedVariables.rho_w) };
	double norm_E    { getNormOfChange(conservedVariables.rho_E, intermediateConservedVariables.rho_E) };
	return ConservedVariablesScalars(norm_rho, norm_rho_u, norm_rho_v, norm_rho_w, norm_E);
}

// Swap the contents of all the arrays of conserved variables and the intermediate arrays, by move-semantics.
// This operation is super fast and needs no extra copy. Only the ownership of the data is changed.
void Mesh::swapConservedVariables()
{
	conservedVariables.rho  .dataSwap(intermediateConservedVariables.rho  );
	conservedVariables.rho_u.dataSwap(intermediateConservedVariables.rho_u);
	conservedVariables.rho_v.dataSwap(intermediateConservedVariables.rho_v);
	conservedVariables.rho_w.dataSwap(intermediateConservedVariables.rho_w);
	conservedVariables.rho_E.dataSwap(intermediateConservedVariables.rho_E);
}

Vector3_u Mesh::getIndices3D(uint index1D) const
{
	uint i = index1D / NJ*NK;
	uint j = index1D % (NJ*NK) / NK;
	uint k = index1D % (NJ*NK) % NK;
	return Vector3_u(i,j,k);
}

uint Mesh::getIndex1D(uint i, uint j, uint k) const
{
	return i*NJ*NK + j*NK + k;
}

Vector3_d Mesh::getNodePosition(uint i, uint j, uint k) const
{
	double x { i * dx }, y { j * dy }, z { k * dz };
	return Vector3_d(x, y, z);
}

IndexBoundingBox Mesh::getSurroundingNodesBox(Vector3_d point) const
{
	uint iNext = static_cast<uint>( ceil( point.x / dx ) );
	uint jNext = static_cast<uint>( ceil( point.y / dy ) );
	uint kNext = static_cast<uint>( ceil( point.z / dz ) );
	IndexBoundingBox surroundingBox(iNext, jNext, kNext);
	surroundingBox.iMin = iNext - 1;
	surroundingBox.jMin = jNext - 1;
	surroundingBox.kMin = kNext - 1;
	return surroundingBox;
}

void Mesh::applyAllBoundaryConditions(double t, const ConfigSettings& params)
{
	for(auto&& boundary : edgeBoundaries)
		boundary->applyBoundaryCondition(*this);

	for(auto&& boundary : immersedBoundaries)
		boundary->applyBoundaryCondition(*this);
}

// Compute the 2-norm of the difference between two arrys. Intended to monitor the change between two consecutive time levels.
// E.g. to check convergence of solution.
double Mesh::getNormOfChange(const Array3D_d& oldValue, const Array3D_d& newValue)
{
	double sumOfSquaredChanges = 0;
	uint numberOfMeshNodes = NI * NJ * NK;
	for(uint i=0; i<numberOfMeshNodes; ++i)
		sumOfSquaredChanges += pow(oldValue(i)-newValue(i), 2);
	double normOfChange = sqrt( dx*dy*dz * sumOfSquaredChanges );
	return normOfChange;
}





