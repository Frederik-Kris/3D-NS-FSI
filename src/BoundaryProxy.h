/*
 * EdgeBoundaryProxy.h
 *
 *  Created on: May 24, 2023
 *      Author: frederk
 */

#ifndef SRC_BOUNDARYPROXY_H_
#define SRC_BOUNDARYPROXY_H_


#include "Boundary.h"

class EdgeBoundaryProxy
{
public:

	virtual unique_ptr<MeshEdgeBoundary> createBoundary(const SubMeshDescriptor& subMeshData) const = 0;
	virtual ~EdgeBoundaryProxy() = default;

	const AxisOrientationEnum normalAxis;	// The axis that is normal to the boundary plane
	const EdgeIndexEnum planeIndex;			// Denotes if the plane is at lowest or highest index side

protected:

	EdgeBoundaryProxy(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex)
	: normalAxis(normalAxis),
	  planeIndex(planeIndex)
	{}
};

class InletBoundaryProxy : public EdgeBoundaryProxy
{
public:

	InletBoundaryProxy(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex, double velocity)
	: EdgeBoundaryProxy(normalAxis, planeIndex),
	  velocity{velocity}
	{}

	unique_ptr<MeshEdgeBoundary> createBoundary(const SubMeshDescriptor& subMeshData) const override
	{ return std::make_unique<InletBoundary>(normalAxis, planeIndex, subMeshData, velocity); }

private:

	double velocity;
};

class OutletBoundaryProxy : public EdgeBoundaryProxy
{
public:

	OutletBoundaryProxy(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex)
	: EdgeBoundaryProxy(normalAxis, planeIndex)
	{}

	unique_ptr<MeshEdgeBoundary> createBoundary(const SubMeshDescriptor& subMeshData) const override
	{ return std::make_unique<OutletBoundary>(normalAxis, planeIndex, subMeshData); }
};

class PeriodicBoundaryProxy : public EdgeBoundaryProxy
{
public:

	PeriodicBoundaryProxy(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex)
	: EdgeBoundaryProxy(normalAxis, planeIndex)
	{}

	unique_ptr<MeshEdgeBoundary> createBoundary(const SubMeshDescriptor& subMeshData) const override
	{ return std::make_unique<PeriodicBoundary>(normalAxis, planeIndex, subMeshData); }
};

class SymmetryBoundaryProxy : public EdgeBoundaryProxy
{
public:

	SymmetryBoundaryProxy(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex)
	: EdgeBoundaryProxy(normalAxis, planeIndex)
	{}

	unique_ptr<MeshEdgeBoundary> createBoundary(const SubMeshDescriptor& subMeshData) const override
	{ return std::make_unique<SymmetryBoundary>(normalAxis, planeIndex, subMeshData); }
};

class ExtrapolationBoundaryProxy : public EdgeBoundaryProxy
{
public:

	ExtrapolationBoundaryProxy(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex)
	: EdgeBoundaryProxy(normalAxis, planeIndex)
	{}

	unique_ptr<MeshEdgeBoundary> createBoundary(const SubMeshDescriptor& subMeshData) const override
	{ return std::make_unique<ExtrapolationBoundary>(normalAxis, planeIndex, subMeshData); }
};

typedef vector<unique_ptr<EdgeBoundaryProxy>> EdgeBoundaryProxyCollection;

// ↑ Mesh edge boundaries above ↑
// ↓ Immersed boundaries below ↓

class ImmersedBoundaryProxy
{
public:

	// Check if the immersed body is inside the given bounding box:
	virtual bool isInsideBox(const SpaceBoundingBox& box) = 0;

	virtual unique_ptr<ImmersedBoundary> createBoundary(const SubMeshDescriptor& subMeshData) const = 0;

	virtual ~ImmersedBoundaryProxy() = default;
};

class CylinderBodyProxy : public ImmersedBoundaryProxy
{
public:

	// Check if cylinder centroid passes through box
	bool isInsideBox(const SpaceBoundingBox& box) override
	{
		return
				( axis == AxisOrientationEnum::x && centroidPosition.y < box.yMax && centroidPosition.y > box.yMin
												 && centroidPosition.z < box.zMax && centroidPosition.z > box.zMin )
			||	( axis == AxisOrientationEnum::y && centroidPosition.x < box.xMax && centroidPosition.x > box.xMin
												 && centroidPosition.z < box.zMax && centroidPosition.z > box.zMin )
			||	( axis == AxisOrientationEnum::z && centroidPosition.x < box.xMax && centroidPosition.x > box.xMin
												 && centroidPosition.y < box.yMax && centroidPosition.y > box.yMin ) ;
	}

	CylinderBodyProxy(Vector3_d centroidPosition, AxisOrientationEnum axis, double radius)
	: centroidPosition(centroidPosition),
	  axis{axis},
	  radius{radius}
	{}

	unique_ptr<ImmersedBoundary> createBoundary(const SubMeshDescriptor& subMeshData) const
	{ return std::make_unique<CylinderBody>(centroidPosition, axis, radius, subMeshData); }

private:

	Vector3_d centroidPosition;
	AxisOrientationEnum axis;
	double radius;
};

class SphereBodyProxy : public ImmersedBoundaryProxy
{
public:

	// Check if sphere center point is inside the given box:
	bool isInsideBox(const SpaceBoundingBox& box) override
	{
		return	centerPosition.x < box.xMax && centerPosition.x > box.xMin
			&&	centerPosition.y < box.yMax && centerPosition.y > box.yMin
			&&	centerPosition.z < box.zMax && centerPosition.z > box.zMin ;
	}

	SphereBodyProxy(Vector3_d centerPosition, double radius)
	: centerPosition(centerPosition),
	  radius{radius}
	{}

	unique_ptr<ImmersedBoundary> createBoundary(const SubMeshDescriptor& subMeshData) const
	{ return std::make_unique<SphereBody>(centerPosition, radius, subMeshData); }

private:

	Vector3_d centerPosition;
	double radius;
};

typedef vector<unique_ptr<ImmersedBoundaryProxy>> ImmersedBoundaryProxyCollection;


#endif /* SRC_BOUNDARYPROXY_H_ */





