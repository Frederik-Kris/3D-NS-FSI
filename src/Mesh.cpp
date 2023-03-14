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
	int nRegionsX = xThresholds.size()-1;
	int nRegionsY = yThresholds.size()-1;
	int nRegionsZ = zThresholds.size()-1;
	subMeshes.setSize(nRegionsX, nRegionsY, nRegionsZ);
	refinementLevels = Array3D<int>(nRegionsX, nRegionsY, nRegionsZ, 0);
	for(const RefinementSpecification& refSpec : params.specifiedRefinementLevels)
	{
		Vector3_i region = refSpec.region;
		for(int indexOffset=refSpec.level-1; indexOffset<=1; ++indexOffset)
		{
			int iMin=max(region.i-indexOffset, 0);
			int iMax=min(region.i+indexOffset, nRegionsX-1);
			int jMin=max(region.j-indexOffset, 0);
			int jMax=min(region.j+indexOffset, nRegionsY-1);
			int kMin=max(region.k-indexOffset, 0);
			int kMax=min(region.k+indexOffset, nRegionsZ-1);
			for(int i=iMin; i<=iMax; ++i)
				for(int j=jMin; j<=jMax; ++j)
					for(int k=kMin; k<=kMax; ++k)
						refinementLevels(i,j,k) = max( refinementLevels(i,j,k), refSpec.level-1 );
		}
		refinementLevels(region.i, region.j, region.k) = max( refinementLevels(region.i, region.j, region.k), refSpec.level );
	}
	for(int i=0; i<nRegionsX; ++i)
	{
		int iMin = round(xThresholds[i  ]/dxBase);
		int iMax = round(xThresholds[i+1]/dxBase);
		if(iMax-iMin < 2)
			throw std::runtime_error("A submesh is too small in the x direction.");
		for(int j=0; j<nRegionsY; ++j)
		{
			int jMin = round(yThresholds[j  ]/dyBase);
			int jMax = round(yThresholds[j+1]/dyBase);
			if(jMax-jMin < 2)
				throw std::runtime_error("A submesh is too small in the y direction.");
			for(int k=0; k<nRegionsZ; ++k)
			{
				int kMin = round(zThresholds[k  ]/dzBase);
				int kMax = round(zThresholds[k+1]/dzBase);
				if(kMax-kMin < 2)
					throw std::runtime_error("A submesh is too small in the z direction.");
				int subMeshSizeX = (iMax-iMin) * pow(2, refinementLevels(i,j,k)) + 1;
				int subMeshSizeY = (jMax-jMin) * pow(2, refinementLevels(i,j,k)) + 1;
				int subMeshSizeZ = (kMax-kMin) * pow(2, refinementLevels(i,j,k)) + 1;
				IndexBoundingBox subMeshIndexBox(subMeshSizeX-1, subMeshSizeY-1, subMeshSizeZ-1);
				SpaceBoundingBox subMeshVolume;
				subMeshVolume.xMin = xThresholds[i];
				subMeshVolume.xMax = xThresholds[i+1];
				subMeshVolume.yMin = yThresholds[j];
				subMeshVolume.yMax = yThresholds[j+1];
				subMeshVolume.zMin = zThresholds[k];
				subMeshVolume.zMax = zThresholds[k+1];
				subMeshes(i,j,k).dx = subMeshSizeX>1 ? (subMeshVolume.xMax - subMeshVolume.xMin)/(subMeshSizeX - 1) : 1;
				subMeshes(i,j,k).dy = subMeshSizeY>1 ? (subMeshVolume.yMax - subMeshVolume.yMin)/(subMeshSizeY - 1) : 1;
				subMeshes(i,j,k).dz = subMeshSizeZ>1 ? (subMeshVolume.zMax - subMeshVolume.zMin)/(subMeshSizeZ - 1) : 1;
				int subMeshNI=subMeshSizeX,subMeshNJ=subMeshSizeY, subMeshNK=subMeshSizeZ;
				EdgeBoundaryCollection subMeshEdgeBCs;
				for(auto&& boundary : edgeBoundaries)
				{
					if(boundary->normalAxis == AxisOrientationEnum::x && boundary->planeIndex == EdgeIndexEnum::min && iMin == 0)
					{
						subMeshEdgeBCs.push_back( boundary->getUniquePtrToCopy() );
						subMeshNI += subMeshEdgeBCs.back()->nExtraNodeLayer();
					}
					if(boundary->normalAxis == AxisOrientationEnum::x && boundary->planeIndex == EdgeIndexEnum::max && iMax == NI-1)
					{
						subMeshEdgeBCs.push_back( boundary->getUniquePtrToCopy() );
						subMeshNI += subMeshEdgeBCs.back()->nExtraNodeLayer();
					}
					if(boundary->normalAxis == AxisOrientationEnum::y && boundary->planeIndex == EdgeIndexEnum::min && jMin == 0)
					{
						subMeshEdgeBCs.push_back( boundary->getUniquePtrToCopy() );
						subMeshNJ += subMeshEdgeBCs.back()->nExtraNodeLayer();
					}
					if(boundary->normalAxis == AxisOrientationEnum::y && boundary->planeIndex == EdgeIndexEnum::max && jMax == NJ-1)
					{
						subMeshEdgeBCs.push_back( boundary->getUniquePtrToCopy() );
						subMeshNJ += subMeshEdgeBCs.back()->nExtraNodeLayer();
					}
					if(boundary->normalAxis == AxisOrientationEnum::z && boundary->planeIndex == EdgeIndexEnum::min && kMin == 0)
					{
						subMeshEdgeBCs.push_back( boundary->getUniquePtrToCopy() );
						subMeshNK += subMeshEdgeBCs.back()->nExtraNodeLayer();
					}
					if(boundary->normalAxis == AxisOrientationEnum::z && boundary->planeIndex == EdgeIndexEnum::max && kMax == NK-1)
					{
						subMeshEdgeBCs.push_back( boundary->getUniquePtrToCopy() );
						subMeshNK += subMeshEdgeBCs.back()->nExtraNodeLayer();
					}
				}
				if(iMin > 0)
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
					subMeshNI += subMeshEdgeBCs.back()->nExtraNodeLayer();
				}
				if(iMax < NI-1)
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
					subMeshNI += subMeshEdgeBCs.back()->nExtraNodeLayer();
				}
				if(jMin > 0)
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
					subMeshNJ += subMeshEdgeBCs.back()->nExtraNodeLayer();
				}
				if(jMax < NJ-1)
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
					subMeshNJ += subMeshEdgeBCs.back()->nExtraNodeLayer();
				}
				if(kMin > 0)
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
					subMeshNK += subMeshEdgeBCs.back()->nExtraNodeLayer();
				}
				if(kMax < NK-1)
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
					subMeshNK += subMeshEdgeBCs.back()->nExtraNodeLayer();
				}
				ImmersedBoundaryCollection subMeshImmersedBCs;
				for(auto&& surface : immersedBoundaries)
				{
					if( surface->isInsideBox(subMeshVolume) )
						subMeshImmersedBCs.push_back( surface->getUniquePtrToCopy() );
				}
				subMeshes(i,j,k).setSize(subMeshSizeX, subMeshSizeY, subMeshSizeZ);
				subMeshes(i,j,k).setBoundaries(subMeshEdgeBCs, subMeshImmersedBCs);
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

// Filters one variable field, i.e., one solution array, 'filterVariable' and stores the filtered
// result in 'variableTemporaryStorage'. Then, the arrays are swapped by move-semantics.
// Only filters if conditions set in config file are met.
void Mesh::applyFilter_ifAppropriate(Array3D_d& filterVariable, Array3D_d& variableTemporaryStorage,
									const ConfigSettings& params, long timeLevel, double t)
{
	if( checkFilterCondition(params, timeLevel, t) )
	{
		// Apply filter to interior nodes:
		Vector3_i nNodes(NI, NJ, NK);
		for (int index1D : indexByType.fluidInterior)
		{
			Vector3_i indices = getIndices3D(index1D, nNodes);
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
	const Vector3_i nMeshNodes(NI, NJ, NK);
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
	for(int i : indexByType.fluidInterior)
		sumOfSquaredDifferences += pow(oldValue(i)-newValue(i), 2);
	double normOfChange = sqrt( dx*dy*dz * sumOfSquaredDifferences );
	return normOfChange;
}





