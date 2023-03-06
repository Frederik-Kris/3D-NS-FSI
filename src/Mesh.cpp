/*
 * Mesh.cpp
 *
 *  Created on: Oct 7, 2022
 *      Author: frederk
 */

#include "Mesh.h"

// Constructor. Taking parameters in ConfigSettings.
Mesh::Mesh(const ConfigSettings& params) :
NI{params.NI}, NJ{params.NJ}, NK{params.NK}
{
	vector<double> xThresholds = params.refineBoundX;
	vector<double> yThresholds = params.refineBoundY;
	vector<double> zThresholds = params.refineBoundZ;
	xThresholds.insert( xThresholds.begin(), 0 );
	yThresholds.insert( yThresholds.begin(), 0 );
	zThresholds.insert( zThresholds.begin(), 0 );
	xThresholds.push_back(params.L_x);
	yThresholds.push_back(params.L_y);
	zThresholds.push_back(params.L_z);
	double dxBase = params.L_x / (NI-1);
	double dyBase = params.L_y / (NJ-1);
	double dzBase = params.L_z / (NK-1);
	for(size_t i=0; i<xThresholds.size()-1; ++i)
	{
		size_t iMin = round(xThresholds[i  ]/dxBase);
		size_t iMax = round(xThresholds[i+1]/dxBase);
		if(iMax-iMin < 2)
			throw std::runtime_error("A submesh is too small in the x direction.");
		for(size_t j=0; j<yThresholds.size()-1; ++j)
		{
			size_t jMin = round(yThresholds[j  ]/dyBase);
			size_t jMax = round(yThresholds[j+1]/dyBase);
			if(jMax-jMin < 2)
				throw std::runtime_error("A submesh is too small in the y direction.");
			for(size_t k=0; k<zThresholds.size()-1; ++k)
			{
				size_t kMin = round(zThresholds[k  ]/dzBase);
				size_t kMax = round(zThresholds[k+1]/dzBase);
				if(kMax-kMin < 2)
					throw std::runtime_error("A submesh is too small in the z direction.");
ØØ; // LAGE SUBMESHES HER?
			}
		}
	}
	if(params.refineBoundX.size() > 0)
	{

	}
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
	edgeBoundaries.push_back(std::make_unique<ExtrapolationBoundary>(AxisOrientationEnum::y, EdgeIndexEnum::min));
	edgeBoundaries.push_back(std::make_unique<ExtrapolationBoundary>(AxisOrientationEnum::y, EdgeIndexEnum::max));
	edgeBoundaries.push_back(std::make_unique<InletBoundary>   (AxisOrientationEnum::x, EdgeIndexEnum::min, params.M_0));
	edgeBoundaries.push_back(std::make_unique<OutletBoundary>  (AxisOrientationEnum::x, EdgeIndexEnum::max));

	// Set grid spacing and origin point offset based on the BC types:
	// In/Out:		spacing=N-1, offset=0
	// Periodic:	spacing=N-2, offset=-gridSpacing or 0
	// Symmetry:	spacing=N-3, offset=-gridSpacing
	dx = params.L_x / (NI-1);
	dy = params.L_y / (NJ-3);
	dz = params.L_z / (NK-2);
	positionOffset.x = 0;
	positionOffset.y = -1 * dy;
	positionOffset.z = -1 * dz;

	// Add immersed bodies:
	Vector3_d cylinderCentroidPosition(params.L_x / 4, params.L_y / 2, 0);
	immersedBoundaries.push_back(std::make_unique<CylinderBody>(cylinderCentroidPosition,
																AxisOrientationEnum::z, 1./2));
//	Vector3_d sphereCenterPoint(params.L_x/4, params.L_y/2, params.L_z/2);
//	immersedBoundaries.push_back(std::make_unique<SphereBody>(sphereCenterPoint, params.L_y/8.1));
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
		applyAllBoundaryConditions(t, params);
	}
}

// Swap the contents of all the arrays of conserved variables and the intermediate arrays, by move-semantics.
// This operation is super fast and needs no extra copy. Only the ownership of the data is changed.
void Mesh::swapConservedVariables()
{
	for(SubMesh& subMesh : subMeshes)
		subMesh.swapConservedVariables();
}

// Loop through boundaries, who are derived classes of MeshEdgeBoundary and ImmersedBoundary,
// and apply their boundary conditions.
void Mesh::applyAllBoundaryConditions(double t, const ConfigSettings& params)
{
	const Vector3_u nMeshNodes(NI, NJ, NK);
	const Vector3_d gridSpacing(dx, dy, dz);
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





