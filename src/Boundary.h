/*
 * BoundaryCondition.h
 *
 *  Created on: Oct 12, 2022
 *      Author: frederk
 */

#ifndef SRC_BOUNDARY_H_
#define SRC_BOUNDARY_H_

struct IndexBoundingBox;
class MeshEdgeBoundary;
class ImmersedBoundary;

#include "includes_and_names.h"
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

typedef std::array<size_t, 8> Array8_u;
typedef std::array<double, 8> Array8_d;
typedef std::array<bool,   8> Array8_b;

struct IndexBoundingBox
{
	size_t iMin, iMax;
	size_t jMin, jMax;
	size_t kMin, kMax;

	IndexBoundingBox(size_t iMax, size_t jMax, size_t kMax)
	: iMin{0}, iMax{iMax},
	  jMin{0}, jMax{jMax},
	  kMin{0}, kMax{kMax}
	{}

	Array8_u asIndexList(const Mesh& mesh);
};

struct InterpolationValues
{
	Array8_d u;
	Array8_d v;
	Array8_d w;
	Array8_d p;
	Array8_d T;
};

struct InterpolationPositions
{
	Array8_d x;
	Array8_d y;
	Array8_d z;
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
	virtual void applyBoundaryCondition(double t, const ConfigSettings& params, Mesh& mesh) = 0;
	const AxisOrientationEnum normalAxis;
	const EdgeIndexEnum planeIndex;
protected:
	vector<size_t> nodeIndices;

	void getAdjacentIndices(size_t index1D, const Mesh& mesh, size_t& boundaryAdjacentIndex, size_t& nextToAdjacentIndex);
	size_t getPeriodicIndex(size_t index1D, const Mesh& mesh);
};

// Class to define inlet boundary condition:
class InletBoundary : public MeshEdgeBoundary
{
public:
	InletBoundary(AxisOrientationEnum normalAxis,
				  EdgeIndexEnum planeIndex,
				  double velocity);
	void applyBoundaryCondition(double t, const ConfigSettings& params, Mesh& mesh) override;
private:
	double velocity;
};

// Class to define outlet boundary condition:
class OutletBoundary : public MeshEdgeBoundary
{
public:
	OutletBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex);
	void applyBoundaryCondition(double t, const ConfigSettings& params, Mesh& mesh) override;
};

// Class to define a periodic boundary condition:
class PeriodicBoundary : public MeshEdgeBoundary
{
public:
	PeriodicBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex);
	void applyBoundaryCondition(double t, const ConfigSettings& params, Mesh& mesh) override;
};

// Class to define a symmetry boundary condition:
class SymmetryBoundary : public MeshEdgeBoundary
{
public:
	SymmetryBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex);
	void applyBoundaryCondition(double t, const ConfigSettings& params, Mesh& mesh) override;
};

// Base class for immersed boundaries.
// Suggested derived types: Adiabatic wall and iso-thermal wall
// or different shapes.
class ImmersedBoundary
{
public:
	ImmersedBoundary();
	virtual ~ImmersedBoundary() = default;
	virtual void identifyRelatedNodes(const ConfigSettings& params, Mesh& mesh) = 0;
	void applyBoundaryCondition(Mesh& mesh);
protected:
	double simplifiedInterpolation(const Array8_d& interpolationValues, const Vector3_u& lowerIndexNode, const Vector3_d& imagePointPosition, const Mesh& mesh);
	PrimitiveVariablesScalars simplifiedInterpolationAll(const InterpolationValues& interpolationValues, const Vector3_u& lowerIndexNode, const Vector3_d& imagePointPosition, const Mesh& mesh);
	PrimitiveVariablesScalars getGhostNodePrimitiveVariables(const PrimitiveVariablesScalars& imagePointPrimVars);
	vector<GhostNode> ghostNodes;
	std::map<size_t, GhostNode*> ghostNodeMap;
};

// Class to define boundary conditions at an immersed cylinder:
class CylinderBody : public ImmersedBoundary
{
public:
	CylinderBody(Vector3_d centroidPosition, AxisOrientationEnum axis, double radius);
	void identifyRelatedNodes(const ConfigSettings& params, Mesh& mesh) override;

private:
	Vector3_d centroidPosition;
	AxisOrientationEnum axis;
	double radius;

	IndexBoundingBox getCylinderBoundingBox(Mesh& mesh) const;
	void getSolidNodesInCylinder(const ConfigSettings& params, vector<size_t>& solidNodeIndices, IndexBoundingBox indicesToCheck, Mesh& mesh);
	void findGhostNodesWithFluidNeighbors(const vector<size_t>& solidNodeIndices, Mesh& mesh);
	void checkIfSurroundingShouldBeGhost(Mesh &mesh, vector<GhostNode>& newGhostNodes, const Vector3_u &surroundingNode);
	vector<GhostNode> setImagePointPositions(vector<GhostNode>& ghostNodesToProcess, Mesh &mesh);
};

// Class to define boundary conditions at an immersed sphere:
class SphereBody : public ImmersedBoundary
{
public:
	SphereBody(Vector3_d centerPosition, double radius);
	void identifyRelatedNodes(const ConfigSettings& params, Mesh& mesh) override;
private:
	Vector3_d centerPosition;
	double radius;
};

#endif /* SRC_BOUNDARY_H_ */




