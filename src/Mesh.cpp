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
	double dxBase = NI>1 ? params.L_x/(NI-1) : 1;
	double dyBase = NJ>1 ? params.L_y/(NJ-1) : 1;
	double dzBase = NK>1 ? params.L_z/(NK-1) : 1;
	nRegions.i = xThresholds.size()-1;
	nRegions.j = yThresholds.size()-1;
	nRegions.k = zThresholds.size()-1;
	subMeshes.setSize(nRegions.i, nRegions.j, nRegions.k);
	refinementLevels = Array3D<int>(nRegions.i, nRegions.j, nRegions.k, 0);
	for(const RefinementSpecification& refSpec : params.specifiedRefinementLevels)
	{
		Vector3_i region = refSpec.region;
		for(int indexOffset=refSpec.level-1; indexOffset<=1; ++indexOffset)
		{
			int iMin=max(region.i-indexOffset, 0);
			int iMax=min(region.i+indexOffset, nRegions.i-1);
			int jMin=max(region.j-indexOffset, 0);
			int jMax=min(region.j+indexOffset, nRegions.j-1);
			int kMin=max(region.k-indexOffset, 0);
			int kMax=min(region.k+indexOffset, nRegions.k-1);
			for(int i=iMin; i<=iMax; ++i)
				for(int j=jMin; j<=jMax; ++j)
					for(int k=kMin; k<=kMax; ++k)
						refinementLevels(i,j,k) = max( refinementLevels(i,j,k), refSpec.level-1 );
		}
		refinementLevels(region.i, region.j, region.k) = max( refinementLevels(region.i, region.j, region.k), refSpec.level );
	}
	double BIGG = std::numeric_limits<double>::max();
	smallestGridSpacings = Vector3_d( BIGG, BIGG, BIGG );
	for(int i=0; i<nRegions.i; ++i)
	{
		int iPrev = round(xThresholds[i  ]/dxBase);
		int iNext = round(xThresholds[i+1]/dxBase);
		if(iNext-iPrev < 2)
			throw std::runtime_error("A submesh is too small in the x direction.");
		for(int j=0; j<nRegions.j; ++j)
		{
			int jPrev = round(yThresholds[j  ]/dyBase);
			int jNext = round(yThresholds[j+1]/dyBase);
			if(jNext-jPrev < 2)
				throw std::runtime_error("A submesh is too small in the y direction.");
			for(int k=0; k<nRegions.k; ++k)
			{
				int kPrev = round(zThresholds[k  ]/dzBase);
				int kNext = round(zThresholds[k+1]/dzBase);
				if(kNext-kPrev < 2)
					throw std::runtime_error("A submesh is too small in the z direction.");
				int subMeshSizeX = (iNext-iPrev) * pow(2, refinementLevels(i,j,k)) + 1;
				int subMeshSizeY = (jNext-jPrev) * pow(2, refinementLevels(i,j,k)) + 1;
				int subMeshSizeZ = (kNext-kPrev) * pow(2, refinementLevels(i,j,k)) + 1;
				IndexBoundingBox subMeshArrayLimits(subMeshSizeX-1, subMeshSizeY-1, subMeshSizeZ-1);
				SpaceBoundingBox subMeshVolume;
				subMeshVolume.xMin = xThresholds[i];
				subMeshVolume.xMax = xThresholds[i+1];
				subMeshVolume.yMin = yThresholds[j];
				subMeshVolume.yMax = yThresholds[j+1];
				subMeshVolume.zMin = zThresholds[k];
				subMeshVolume.zMax = zThresholds[k+1];
				EdgeBoundaryCollection subMeshEdgeBCs;
				for(auto&& boundary : edgeBoundaries)
				{
					if(boundary->normalAxis == AxisOrientationEnum::x && boundary->planeIndex == EdgeIndexEnum::min && iPrev == 0)
					{
						subMeshEdgeBCs.push_back( boundary->getUniquePtrToCopy() );
						subMeshArrayLimits.iMin -= subMeshEdgeBCs.back()->nExtraNodeLayer();
					}
					if(boundary->normalAxis == AxisOrientationEnum::x && boundary->planeIndex == EdgeIndexEnum::max && iNext == NI-1)
					{
						subMeshEdgeBCs.push_back( boundary->getUniquePtrToCopy() );
						subMeshArrayLimits.iMax += subMeshEdgeBCs.back()->nExtraNodeLayer();
					}
					if(boundary->normalAxis == AxisOrientationEnum::y && boundary->planeIndex == EdgeIndexEnum::min && jPrev == 0)
					{
						subMeshEdgeBCs.push_back( boundary->getUniquePtrToCopy() );
						subMeshArrayLimits.jMin -= subMeshEdgeBCs.back()->nExtraNodeLayer();
					}
					if(boundary->normalAxis == AxisOrientationEnum::y && boundary->planeIndex == EdgeIndexEnum::max && jNext == NJ-1)
					{
						subMeshEdgeBCs.push_back( boundary->getUniquePtrToCopy() );
						subMeshArrayLimits.jMax += subMeshEdgeBCs.back()->nExtraNodeLayer();
					}
					if(boundary->normalAxis == AxisOrientationEnum::z && boundary->planeIndex == EdgeIndexEnum::min && kPrev == 0)
					{
						subMeshEdgeBCs.push_back( boundary->getUniquePtrToCopy() );
						subMeshArrayLimits.kMin -= subMeshEdgeBCs.back()->nExtraNodeLayer();
					}
					if(boundary->normalAxis == AxisOrientationEnum::z && boundary->planeIndex == EdgeIndexEnum::max && kNext == NK-1)
					{
						subMeshEdgeBCs.push_back( boundary->getUniquePtrToCopy() );
						subMeshArrayLimits.kMax += subMeshEdgeBCs.back()->nExtraNodeLayer();
					}
				}
				if(iPrev > 0)
				{
					if( refinementLevels(i,j,k) > refinementLevels(i-1,j,k) )
						subMeshEdgeBCs.push_back( std::make_unique<InterfaceToCoarserSubMesh>(
								AxisOrientationEnum::x,
								EdgeIndexEnum::min,
								subMeshes(i-1,j,k).flowVariableReferences ) );
					else if( refinementLevels(i,j,k) < refinementLevels(i-1,j,k) )
						subMeshEdgeBCs.push_back( std::make_unique<InterfaceToFinerSubMesh>(
								AxisOrientationEnum::x,
								EdgeIndexEnum::min,
								subMeshes(i-1,j,k).flowVariableReferences ) );
					else // same level
						subMeshEdgeBCs.push_back( std::make_unique<InterfaceToEqualLevelSubMesh>(
								AxisOrientationEnum::x,
								EdgeIndexEnum::min,
								subMeshes(i-1,j,k).flowVariableReferences ) );
					subMeshArrayLimits.iMin -= subMeshEdgeBCs.back()->nExtraNodeLayer();
				}
				if(iNext < NI-1)
				{
					if( refinementLevels(i,j,k) > refinementLevels(i+1,j,k) )
						subMeshEdgeBCs.push_back( std::make_unique<InterfaceToCoarserSubMesh>(
								AxisOrientationEnum::x,
								EdgeIndexEnum::max,
								subMeshes(i+1,j,k).flowVariableReferences ) );
					else if( refinementLevels(i,j,k) < refinementLevels(i+1,j,k) )
						subMeshEdgeBCs.push_back( std::make_unique<InterfaceToFinerSubMesh>(
								AxisOrientationEnum::x,
								EdgeIndexEnum::max,
								subMeshes(i+1,j,k).flowVariableReferences ) );
					else // same level
						subMeshEdgeBCs.push_back( std::make_unique<InterfaceToEqualLevelSubMesh>(
								AxisOrientationEnum::x,
								EdgeIndexEnum::max,
								subMeshes(i+1,j,k).flowVariableReferences ) );
					subMeshArrayLimits.iMax += subMeshEdgeBCs.back()->nExtraNodeLayer();
				}
				if(jPrev > 0)
				{
					if( refinementLevels(i,j,k) > refinementLevels(i,j-1,k) )
						subMeshEdgeBCs.push_back( std::make_unique<InterfaceToCoarserSubMesh>(
								AxisOrientationEnum::y,
								EdgeIndexEnum::min,
								subMeshes(i,j-1,k).flowVariableReferences ) );
					else if( refinementLevels(i,j,k) < refinementLevels(i,j-1,k) )
						subMeshEdgeBCs.push_back( std::make_unique<InterfaceToFinerSubMesh>(
								AxisOrientationEnum::y,
								EdgeIndexEnum::min,
								subMeshes(i,j-1,k).flowVariableReferences ) );
					else // same level
						subMeshEdgeBCs.push_back( std::make_unique<InterfaceToEqualLevelSubMesh>(
								AxisOrientationEnum::y,
								EdgeIndexEnum::min,
								subMeshes(i,j-1,k).flowVariableReferences ) );
					subMeshArrayLimits.jMin -= subMeshEdgeBCs.back()->nExtraNodeLayer();
				}
				if(jNext < NJ-1)
				{
					if( refinementLevels(i,j,k) > refinementLevels(i,j+1,k) )
						subMeshEdgeBCs.push_back( std::make_unique<InterfaceToCoarserSubMesh>(
								AxisOrientationEnum::y,
								EdgeIndexEnum::max,
								subMeshes(i,j+1,k).flowVariableReferences ) );
					else if( refinementLevels(i,j,k) < refinementLevels(i,j+1,k) )
						subMeshEdgeBCs.push_back( std::make_unique<InterfaceToFinerSubMesh>(
								AxisOrientationEnum::y,
								EdgeIndexEnum::max,
								subMeshes(i,j+1,k).flowVariableReferences ) );
					else // same level
						subMeshEdgeBCs.push_back( std::make_unique<InterfaceToEqualLevelSubMesh>(
								AxisOrientationEnum::y,
								EdgeIndexEnum::max,
								subMeshes(i,j+1,k).flowVariableReferences ) );
					subMeshArrayLimits.jMax += subMeshEdgeBCs.back()->nExtraNodeLayer();
				}
				if(kPrev > 0)
				{
					if( refinementLevels(i,j,k) > refinementLevels(i,j,k-1) )
						subMeshEdgeBCs.push_back( std::make_unique<InterfaceToCoarserSubMesh>(
								AxisOrientationEnum::z,
								EdgeIndexEnum::min,
								subMeshes(i,j,k-1).flowVariableReferences ) );
					else if( refinementLevels(i,j,k) < refinementLevels(i,j,k-1) )
						subMeshEdgeBCs.push_back( std::make_unique<InterfaceToFinerSubMesh>(
								AxisOrientationEnum::z,
								EdgeIndexEnum::min,
								subMeshes(i,j,k-1).flowVariableReferences ) );
					else // same level
						subMeshEdgeBCs.push_back( std::make_unique<InterfaceToEqualLevelSubMesh>(
								AxisOrientationEnum::z,
								EdgeIndexEnum::min,
								subMeshes(i,j,k-1).flowVariableReferences ) );
					subMeshArrayLimits.kMin -= subMeshEdgeBCs.back()->nExtraNodeLayer();
				}
				if(kNext < NK-1)
				{
					if( refinementLevels(i,j,k) > refinementLevels(i,j,k+1) )
						subMeshEdgeBCs.push_back( std::make_unique<InterfaceToCoarserSubMesh>(
								AxisOrientationEnum::z,
								EdgeIndexEnum::max,
								subMeshes(i,j,k+1).flowVariableReferences ) );
					else if( refinementLevels(i,j,k) < refinementLevels(i,j,k+1) )
						subMeshEdgeBCs.push_back( std::make_unique<InterfaceToFinerSubMesh>(
								AxisOrientationEnum::z,
								EdgeIndexEnum::max,
								subMeshes(i,j,k+1).flowVariableReferences ) );
					else // same level
						subMeshEdgeBCs.push_back( std::make_unique<InterfaceToEqualLevelSubMesh>(
								AxisOrientationEnum::z,
								EdgeIndexEnum::max,
								subMeshes(i,j,k+1).flowVariableReferences ) );
					subMeshArrayLimits.kMax += subMeshEdgeBCs.back()->nExtraNodeLayer();
				}
				ImmersedBoundaryCollection subMeshImmersedBCs;
				for(auto&& surface : immersedBoundaries)
				{
					if( surface->isInsideBox(subMeshVolume) )
						subMeshImmersedBCs.push_back( surface->getUniquePtrToCopy() );
				}
				subMeshes(i,j,k).regionID = Vector3_i(i,j,k);
				subMeshes(i,j,k).setSize(subMeshSizeX, subMeshSizeY, subMeshSizeZ, subMeshArrayLimits, subMeshVolume);
				subMeshes(i,j,k).setBoundaries(subMeshEdgeBCs, subMeshImmersedBCs);
				smallestGridSpacings.x = min(smallestGridSpacings.x, subMeshes(i,j,k).gridSpacings.x);
				smallestGridSpacings.y = min(smallestGridSpacings.y, subMeshes(i,j,k).gridSpacings.y);
				smallestGridSpacings.z = min(smallestGridSpacings.z, subMeshes(i,j,k).gridSpacings.z);
			}
		}
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

	// Add immersed bodies:
	Vector3_d cylinderCentroidPosition(params.L_x / 4, params.L_y / 2, 0);
	immersedBoundaries.push_back(std::make_unique<CylinderBody>(cylinderCentroidPosition,
																AxisOrientationEnum::z, 1./2));
//	Vector3_d sphereCenterPoint(params.L_x/4, params.L_y/2, params.L_z/2);
//	immersedBoundaries.push_back(std::make_unique<SphereBody>(sphereCenterPoint, params.L_y/8.1));
}

// Get SubMeshDescriptors about the sub-meshes given by the given indices.
Array3D<SubMeshDescriptor> Mesh::getSubMeshDescriptors(const IndexBoundingBox& subMeshIndices)
{
	Array3D<SubMeshDescriptor> neighborRegions(subMeshIndices);
	for(int i=subMeshIndices.iMin; i<=subMeshIndices.iMax; ++i)
		for(int j=subMeshIndices.jMin; j<=subMeshIndices.jMax; ++j)
			for(int k=subMeshIndices.kMin; k<=subMeshIndices.kMax; ++k)
			{
				neighborRegions(i,j,k) = SubMeshDescriptor(subMeshes(i,j,k).nNodes,
														   subMeshes(i,j,k).arrayLimits,
														   subMeshes(i,j,k).gridSpacings,
														   subMeshes(i,j,k).nodeType,
														   subMeshes(i,j,k).flowVariableReferences,
														   refinementLevels(i,j,k)
														   );
			}
	return neighborRegions;
}

void Mesh::categorizeNodes(const ConfigSettings& params)
{
	for(int i=0; i<nRegions.i; ++i)
		for(int j=0; j<nRegions.j; ++j)
			for(int k=0; k<nRegions.k; ++k)
			{
				IndexBoundingBox neighborSubMeshIndices = IndexBoundingBox::boxAroundNode(Vector3_i(i,j,k), 1);
				neighborSubMeshIndices = neighborSubMeshIndices.intersection( IndexBoundingBox(nRegions.i, nRegions.j, nRegions.k) );
				Array3D<SubMeshDescriptor> neighborSubmeshes = getSubMeshDescriptors(neighborSubMeshIndices);
				subMeshes(i,j,k).categorizeNodes(params, neighborSubMeshVariableReferences, neighborSubMeshRefinementLevels);
			}
}

// Check if the 2nd order filter should be applied at the current time
bool Mesh::checkFilterCondition(const ConfigSettings& params, long timeLevel, double t)
{
	// Loop through the list of specified times in the config file, and see if we have surpassed them:
	int counter = 0;
	for (double intervalChangeTime : params.filterIntervalChangeTimes)
		if (t > intervalChangeTime)
			++counter;
	bool filterNow = false;
	// If we passed fewer thresholds than there are entries in filterIntervals, then select the correct
	// filter interval and check if timeLevel+1 is divisible by that interval:
	if ( counter < params.filterIntervals.size() )
	{
		int currentFilterInterval = params.filterIntervals.at(counter);
		if(currentFilterInterval > 0)
			if( (timeLevel+1) % currentFilterInterval == 0 )
				filterNow = true;
	}
	return filterNow;
}

// Filters density, ρ, and stores the filtered result in 'conservedVariablesOld.rho'.
// Then, the arrays are swapped by move-semantics.
// Only filters if conditions set in config file are met.
void Mesh::applyFilter_ifAppropriate(const ConfigSettings& params, long timeLevel, double t)
{
	if( checkFilterCondition(params, timeLevel, t) )
	{
		for(SubMesh& subMesh : subMeshes)
		{
			Array3D<double>& rhoTemporary = subMesh.conservedVariablesOld.rho;
			Array3D<double>& rho = subMesh.conservedVariables.rho;
			// Apply filter to interior nodes:
			for (int index1D : subMesh.indexByType.fluidInterior)
			{
				Vector3_i indices = getIndices3D(index1D, subMesh.arrayLimits);
				rhoTemporary(index1D) = (6*rho(indices)
							 	 	 	 + rho(indices.i+1, indices.j, indices.k)
										 + rho(indices.i-1, indices.j, indices.k)
										 + rho(indices.i, indices.j+1, indices.k)
										 + rho(indices.i, indices.j-1, indices.k)
										 + rho(indices.i, indices.j, indices.k+1)
										 + rho(indices.i, indices.j, indices.k-1)
										 ) / 12;
			}
			rho.dataSwap(rhoTemporary);	// Swap the arrays using move-sematics (super-fast)
			applyAllBoundaryConditions(t, params);
		}
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
	for(SubMesh& subMesh : subMeshes)
	{

		for(auto&& boundary : edgeBoundaries)
			boundary->applyBoundaryCondition(t, params);

		SubMeshDescriptor subMeshData = subMesh.getSubMeshDescriptor();
		for(auto&& boundary : immersedBoundaries)
			boundary->applyBoundaryCondition(params, meshData);
	}
}

// Compute the 2-norm of the difference between two arrys. Intended to monitor the change between two consecutive time levels,
// .g. to check convergence of solution. Could also be used to check error.
double Mesh::getNormOfChange(const Array3D_d& oldValue, const Array3D_d& newValue)
{
	double sumOfSquaredDifferences = 0;
	for(int i : indexByType.fluidInterior)
		sumOfSquaredDifferences += pow(oldValue(i)-newValue(i), 2);
	double normOfChange = sqrt( dx*dy*dz * sumOfSquaredDifferences );
	return normOfChange;
}





