/*
 * BoundaryCondition.h
 *
 *  Created on: Oct 12, 2022
 *      Author: frederk
 */

#ifndef SRC_BOUNDARY_H_
#define SRC_BOUNDARY_H_

#include "includes_and_names.h"
#include "CustomExceptions.h"
#include "Array3D.h"
#include "ConfigSettings.h"
#include "Mesh.h"
#include "Node.h"


enum class AxisOrientationEnum
{
	x, y, z
};

enum class EdgeIndexEnum
{
	min, max
};

struct IndexBoundingBox
{
	uint iMin, iMax;
	uint jMin, jMax;
	uint kMin, kMax;

	IndexBoundingBox(uint iMax, uint jMax, uint kMax)
	: iMin{0}, iMax{iMax},
	  jMin{0}, jMax{jMax},
	  kMin{0}, kMax{kMax}
	  {}
};


// Base class for boundaries at the the edge of the Cartesian computational mesh.
// Suggested derived types: Inlet, outlet, walls, periodic, symmetry...
class MeshEdgeBoundary
{
public:
	MeshEdgeBoundary(AxisOrientationEnum normalAxis,
					 EdgeIndexEnum planeIndex);
	virtual ~MeshEdgeBoundary() = default;
	void identifyOwnedNodes(IndexBoundingBox& unclaimedNodes, Mesh& mesh);
	virtual void applyBoundaryCondition();
	const AxisOrientationEnum normalAxis;
	const EdgeIndexEnum planeIndex;
private:
	vector<uint> nodeIndices;
};

// Class to define inlet boundary condition:
class InletBoundary : public MeshEdgeBoundary
{
public:
	InletBoundary(AxisOrientationEnum normalAxis,
				  EdgeIndexEnum planeIndex,
				  double velocity);
	void applyBoundaryCondition() override;
private:
	double velocity;
};

// Class to define outlet boundary condition:
class OutletBoundary : public MeshEdgeBoundary
{
public:
	OutletBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex);
	void applyBoundaryCondition() override;
};

// Class to define a periodic boundary condition:
class PeriodicBoundary : public MeshEdgeBoundary
{
public:
	PeriodicBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex);
	void applyBoundaryCondition() override;
};

// Class to define a symmetry boundary condition:
class SymmetryBoundary : public MeshEdgeBoundary
{
public:
	SymmetryBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex);
	void applyBoundaryCondition() override;
};

// Base class for immersed boundaries.
// Suggested derived types: Adiabatic wall and iso-thermal wall
// or different shapes.
class ImmersedBoundary
{
public:
	ImmersedBoundary();
	virtual ~ImmersedBoundary() = default;
	virtual void identifyRelatedNodes(Mesh& mesh);
	virtual void applyBoundaryCondition();
private:
	vector<GhostNode> ghostNodes;
};

// Class to define boundary conditions at an immersed cylinder:
class CylinderBody : public ImmersedBoundary
{
public:
	CylinderBody(Vector3_d centroidPosition, AxisOrientationEnum axis, double radius);
	void identifyRelatedNodes(const ConfigSettings& params, Mesh& mesh) override;

	void applyBoundaryCondition() override;
private:
	Vector3_d centroidPosition;
	AxisOrientationEnum axis;
	double radius;

	IndexBoundingBox getCylinderBoundingBox(Mesh& mesh) const;
	void getSolidNodesInCylinder(const ConfigSettings& params, vector<uint>& solidNodeIndices, IndexBoundingBox indicesToCheck, Mesh& mesh);
	void findGhostNodes(const vector<uint>& solidNodeIndices, Mesh& mesh);
};

// Class to define boundary conditions at an immersed sphere:
class SphereBody : public ImmersedBoundary
{
public:
	SphereBody(Vector3_d centerPosition, double radius);
	void identifyRelatedNodes(Mesh& mesh) override;
	void applyBoundaryCondition() override;
private:
	Vector3_d centerPosition;
	double radius;
};

#endif /* SRC_BOUNDARY_H_ */




