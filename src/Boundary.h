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
#include <eigen3/Eigen/LU>


enum class AxisOrientationEnum
{
	x, y, z
};

enum class EdgeIndexEnum
{
	min, max
};

typedef std::array<size_t, 8> Array8_u;
typedef std::array<bool,   8> Array8_b;
typedef Eigen::Array<double, 8, 1> Vector8_d;
typedef Eigen::Matrix<double, 8, 8> Matrix8x8_d;

struct IndexBoundingBox
{
	size_t iMin, iMax;	// Indices in x-direction
	size_t jMin, jMax;	// Indices in y-direction
	size_t kMin, kMax;	// Indices in z-direction

	IndexBoundingBox(size_t iMax, size_t jMax, size_t kMax)
	: iMin{0}, iMax{iMax},
	  jMin{0}, jMax{jMax},
	  kMin{0}, kMax{kMax}
	{}

	IndexBoundingBox() = default;

	Array8_u asIndexList(const Mesh& mesh) const;
};

struct InterpolationValues
{
	Vector8_d u;	// Velocity component in x-direction
	Vector8_d v;	// Velocity component in y-direction
	Vector8_d w;	// Velocity component in z-direction
	Vector8_d p;	// Pressure
	Vector8_d T;	// Temperature
};

struct InterpolationPositions
{
	Vector8_d x;
	Vector8_d y;
	Vector8_d z;
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

	void getAdjacentIndices(size_t index1D,
							const Mesh& mesh,
							size_t& boundaryAdjacentIndex,
							size_t& nextToAdjacentIndex);

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

	void applyBoundaryCondition(Mesh& mesh, const ConfigSettings& params);

protected:

	vector<GhostNode> ghostNodes;
	std::map<size_t, GhostNode*> ghostNodeMap;

	void findGhostNodesWithFluidNeighbors(const vector<size_t>& solidNodeIndices, Mesh& mesh);

	void checkIfSurroundingShouldBeGhost(Mesh &mesh,
										 vector<GhostNode>& newGhostNodes,
										 const Vector3_u &surroundingNode);

	vector<GhostNode> setImagePointPositions(vector<GhostNode>& ghostNodesToProcess, Mesh &mesh);

private:
	virtual Vector3_d getNormalProbe(const Vector3_d& ghostNodePosition) = 0;

	double simplifiedInterpolation(const Vector8_d& interpolationValues,
								   const Vector3_u& lowerIndexNode,
								   const Vector3_d& imagePointPosition,
								   const Mesh& mesh);

	PrimitiveVariablesScalars simplifiedInterpolationAll(const InterpolationValues& interpolationValues,
														 const Vector3_u& lowerIndexNode,
														 const Vector3_d& imagePointPosition,
														 const Mesh& mesh);

	PrimitiveVariablesScalars getGhostNodePrimitiveVariables(const PrimitiveVariablesScalars& imagePointPrimVars);

	void populateVandermondeDirichlet(const InterpolationPositions& interpolationPoints,
									  Matrix8x8_d& vandermonde);

	void populateVandermondeNeumann(const InterpolationPositions& interpolationPoints,
									const Array8_b& ghostFlags,
									const vector<Vector3_d>& unitNormals,
									Matrix8x8_d& vandermonde);

	double trilinearInterpolation(const Vector8_d& interpolationValues,
													 const Vector3_d& imagePoint,
													 const Matrix8x8_d& vandermonde);

	PrimitiveVariablesScalars trilinearInterpolationAll(const InterpolationValues&,
														const Vector3_d& imagePoint,
														const Matrix8x8_d& vandermondeDirichlet,
														const Matrix8x8_d& vandermondeNeumann);

	void setInterpolationValuesFluidNode(uint counter, size_t surroundingNodeIndex1D,
										const Mesh &mesh,
										InterpolationValues& interpolationValues,		// OUTPUT
										InterpolationPositions& interpolationPositions);// OUTPUT

	void setInterpolationValuesGhostNode(
			uint counter,
			size_t surroundingNodeIndex1D,
			vector<Vector3_d>& unitNormals,					// <-
			InterpolationValues& interpolationValues,		// <- Output
			InterpolationPositions& interpolationPositions);// <-

	void setInterpolationValues(
			const IndexBoundingBox& surroundingNodes,
			const Mesh& mesh,
			InterpolationValues& interpolationValues,		// <-
			InterpolationPositions& interpolationPositions,	// <-
			Array8_b& ghostFlag,							// <- Output
			bool& allSurroundingAreFluid,					// <-
			vector<Vector3_d>& unitNormals);				// <-

	PrimitiveVariablesScalars interpolatePrimitiveVariables(
			const InterpolationValues& interpolationValues,
			const InterpolationPositions& interpolationPositions,
			bool allSurroundingAreFluid,
			const Array8_b& ghostFlag,
			const GhostNode& ghostNode,
			const IndexBoundingBox& surroundingNodes,
			const vector<Vector3_d>& unitNormals,
			const Mesh& mesh);
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

	Vector3_d getNormalProbe(const Vector3_d& ghostNodePosition) override;

	IndexBoundingBox getCylinderBoundingBox(Mesh& mesh) const;

	void getSolidNodesInCylinder(const ConfigSettings& params,
								 vector<size_t>& solidNodeIndices,
								 IndexBoundingBox indicesToCheck,
								 Mesh& mesh);
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

	Vector3_d getNormalProbe(const Vector3_d& ghostNodePosition) override;

	IndexBoundingBox getSphereBoundingBox(Mesh& mesh) const;

	void getSolidNodesInSphere(const ConfigSettings& params,
								 vector<size_t>& solidNodeIndices,
								 IndexBoundingBox indicesToCheck,
								 Mesh& mesh);
};

#endif /* SRC_BOUNDARY_H_ */




