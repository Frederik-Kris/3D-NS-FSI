/*
 * BoundaryCondition.h
 *
 *  Created on: Oct 12, 2022
 *      Author: frederk
 */

#ifndef SRC_BOUNDARY_H_
#define SRC_BOUNDARY_H_

#include "includes_and_names.h"
#include "Array3D.h"
#include "ConfigSettings.h"
#include "FlowVariableGroupStructs.h"
#include "Node.h"
#include <eigen3/Eigen/LU>
#include <map>


// The three spatial coordinate axes
enum class AxisOrientationEnum
{
	x, y, z
};

// Denotes if a boundary plane is at the lowest or highest coordinate value.
enum class EdgeIndexEnum
{
	min=-1, max=1
};

typedef Eigen::Array<double, 8, 1> Vector8_d;
typedef Eigen::Matrix<double, 8, 8> Matrix8x8_d;

// Package struct with data from a sub-mesh.
// Shortens argument lists, and lets Boundary not depend on Mesh
struct SubMeshDescriptor
{
	const Vector3_i& nNodes;
	const Vector3_d& spacings;
	const IndexBoundingBox& arrayLimits;
	const SpaceBoundingBox& boundingBox;
	Array3D<NodeTypeEnum>& nodeType;
	AllFlowVariablesArrayGroup& flowVariables;
	const int refinementLevel;
	const Vector3_i& regionID;

	SubMeshDescriptor(const Vector3_i& nNodes,
					  const Vector3_d& gridSpacings,
				   	  const IndexBoundingBox& arrayLimits,
					  const SpaceBoundingBox& boundingBox,
					  Array3D<NodeTypeEnum>& nodeTypeArray,
					  AllFlowVariablesArrayGroup& flowVariables,
					  int refinementLevel,
					  const Vector3_i& regionID)
	: nNodes{nNodes},
	  spacings{gridSpacings},
	  arrayLimits(arrayLimits),
	  boundingBox{boundingBox},
	  nodeType{nodeTypeArray},
	  flowVariables{flowVariables},
	  refinementLevel{refinementLevel},
	  regionID{regionID}
	{}
};

// Package struct with benchmark integral sizes, for cylinder and sphere test cases.
struct IntegralProperties
{
	double drag;
	double lift;
	vector<double> separationAngles;
};

// Package struct with the variables to interpolate in the image points.
// 8 values for each variables, since there are 8 nodes surrounding each image point.
struct InterpolationValues
{
	Vector8_d u;	// Velocity component in x-direction
	Vector8_d v;	// Velocity component in y-direction
	Vector8_d w;	// Velocity component in z-direction
	Vector8_d p;	// Pressure
	Vector8_d T;	// Temperature
};

// Package struct with the positions of interpolation points.
// 8 values for each coordinate, since there are 8 nodes surrounding each image point.
struct InterpolationPositions
{
	Vector8_d x;
	Vector8_d y;
	Vector8_d z;
};

struct SpecialTreatmentNodeInfo
{
	SpecialTreatmentNodeInfo(int ownNode,
							 const vector<int> borrowedNodes,
							 const AllFlowVariablesArrayGroup& subMeshToBorrow) :
		ownedNode{ownNode},
		borrowedNodes(borrowedNodes),
		borrowedSubMesh(subMeshToBorrow)
	{}

	int ownedNode;
	vector<int> borrowedNodes;
	const AllFlowVariablesArrayGroup& borrowedSubMesh;
};

struct RelatedNodesCollection
{
	IndexBoundingBox regular;
	vector<SpecialTreatmentNodeInfo> specialTreatment;
};

// Base class for boundaries at the the edge of the Cartesian computational mesh.
// Suggested derived types: Inlet, outlet, walls, periodic, symmetry...
class MeshEdgeBoundary
{
public:
	MeshEdgeBoundary(AxisOrientationEnum normalAxis,
					 EdgeIndexEnum planeIndex,
					 NodeTypeEnum ownedNodesType,
					 const SubMeshDescriptor& subMeshData);

	virtual ~MeshEdgeBoundary() = default;

	// Deep copy the boundary and return unique pointer to the copy.
	// Actual referenced object is not MeshEdgeBoundary, but correct derived class.
	virtual std::unique_ptr<MeshEdgeBoundary> getUniquePtrToCopy() = 0;

	// Get the number of extra node layers (zero or one) required by this BC.
	virtual int nExtraNodeLayer() = 0; // ← PURE virtual

	virtual void identifyRelatedNodes(IndexBoundingBox& unclaimedNodes,
									  const Array3D<SubMeshDescriptor>& neighborSubMeshes);

	/**
	 * Applies the boundary condition to the submesh.
	 *
	 * This is a pure virtual function that needs to be implemented in derived classes.
	 * It defines the behavior of applying the boundary condition to the submesh.
	 * The boundary conditions set the flow variables in its owned nodes.
	 *
	 * @param t The current time.
	 * @param params The configuration settings for the simulation.
	 */
	virtual void applyBoundaryCondition(double t, const ConfigSettings& params) = 0;


	Vector3_i getOutwardNormal() const;

	AxisOrientationEnum normalAxis;			// The axis that is normal to the boundary plane
	EdgeIndexEnum planeIndex;				// Denotes if the plane is at lowest or highest index side

protected:

	IndexBoundingBox peelOffBoundaryPlane(IndexBoundingBox&,
										  AxisOrientationEnum normalAxis,
										  EdgeIndexEnum planeIndex);

	void getAdjacentIndices(const Vector3_i& boundaryNode, 			// ← Input
			  	  	  	  	int& boundaryAdjacentIndex, 			// ← Output
							int& nextToAdjacentIndex);				// ← Output

	Vector3_i getPeriodicIndex(const Vector3_i& boundaryNode);

	NodeTypeEnum ownedNodesType;			// The node type that is assigned to the owned nodes
	RelatedNodesCollection relatedNodes;	// The nodes that the BC at this boundary uses
	SubMeshDescriptor subMeshData;			// Info and flow variables for the parent submesh of the boundary
};

// Class to define inlet boundary condition:
class InletBoundary : public MeshEdgeBoundary
{
public:
	InletBoundary(AxisOrientationEnum normalAxis,
				  EdgeIndexEnum planeIndex,
				  const SubMeshDescriptor& subMeshData,
				  double velocity);

	std::unique_ptr<MeshEdgeBoundary> getUniquePtrToCopy() override
	{ return std::make_unique<InletBoundary>(*this); }

	int nExtraNodeLayer() override { return 0; }

	void filterInletDensity(const Vector3_i& nMeshNodes,				// ← Input
	   	   	   	   	   	    const ConfigSettings& params,				// ← Input
							AllFlowVariablesArrayGroup& flowVariables);	// ← Output

	void applyBoundaryCondition(double t, const ConfigSettings& params)	override;

private:
	double velocity;	// The prescribed inlet velocity, that we reach after initial ramp-up.
};

// Class to define outlet boundary condition:
class OutletBoundary : public MeshEdgeBoundary
{
public:
	OutletBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex, const SubMeshDescriptor& subMeshData);

	std::unique_ptr<MeshEdgeBoundary> getUniquePtrToCopy() override
	{ return std::make_unique<OutletBoundary>(*this); }

	int nExtraNodeLayer() override { return 0; }

	void applyBoundaryCondition(double t, const ConfigSettings& params) override;
};

// Class to define a periodic boundary condition:
class PeriodicBoundary : public MeshEdgeBoundary
{
public:
	PeriodicBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex, const SubMeshDescriptor& subMeshData);

	std::unique_ptr<MeshEdgeBoundary> getUniquePtrToCopy() override
	{ return std::make_unique<PeriodicBoundary>(*this); }

	int nExtraNodeLayer() override { return 1; }

	void applyBoundaryCondition(double t, const ConfigSettings& params) override;
};

// Class to define a symmetry boundary condition:
class SymmetryBoundary : public MeshEdgeBoundary
{
public:
	SymmetryBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex, const SubMeshDescriptor& subMeshData);

	std::unique_ptr<MeshEdgeBoundary> getUniquePtrToCopy() override
	{ return std::make_unique<SymmetryBoundary>(*this); }

	int nExtraNodeLayer() override { return 1; }

	void applyBoundaryCondition(double t, const ConfigSettings& params)	override;
};

// Class to define a boundary condition where all variables are extrapolated:
class ExtrapolationBoundary : public MeshEdgeBoundary
{
public:
	ExtrapolationBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex, const SubMeshDescriptor& subMeshData);

	std::unique_ptr<MeshEdgeBoundary> getUniquePtrToCopy() override
	{ return std::make_unique<ExtrapolationBoundary>(*this); }

	int nExtraNodeLayer() override { return 1; }

	void applyBoundaryCondition(double t, const ConfigSettings& params)	override;
};

// Class to define boundary condition for a submesh's co-boundary with another submesh.
// Handles hanging nodes and fetching values from neighbor submesh etc.
class SubmeshInterfaceBoundary : public MeshEdgeBoundary
{
public:
	
	SubmeshInterfaceBoundary(AxisOrientationEnum normalAxis,
							 EdgeIndexEnum planeIndex,
							 const SubMeshDescriptor& subMeshData,
							 const SubMeshDescriptor& neighborSubMesh) :
		MeshEdgeBoundary(normalAxis, planeIndex, NodeTypeEnum::FluidGhost, subMeshData),
		neighborSubMesh(neighborSubMesh)
	{}

	std::unique_ptr<MeshEdgeBoundary> getUniquePtrToCopy() override
	{ return std::make_unique<SubmeshInterfaceBoundary>(*this); }

	int nExtraNodeLayer() override;

	void identifyRelatedNodes(IndexBoundingBox& unclaimedNodes,
							  const Array3D<SubMeshDescriptor>& neighborSubMeshes) override;

	void applyBoundaryCondition(double t, const ConfigSettings& params)	override;

protected:

	bool subMeshHasAppliedBC(const Vector3_i& thisRegionID, const Vector3_i& otherRegionID) const;

	Vector3_d indexInOtherSubMesh(const Vector3_i& node,
								  const SubMeshDescriptor& otherSubMesh,
								  const SubMeshDescriptor& thisSubMesh) const;

	Vector3_i whoOwnsThisNode(const Vector3_i& node,
							  const Vector3_i& regionID,
							  const Array3D<SubMeshDescriptor>& subMeshes) const;

	vector<int> getNodesAroundPoint(const Vector3_d& point, const IndexBoundingBox& arrayLimits) const;

	double getMeanNodeValue(vector<int> nodes, const Array3D<double>& flowVariable) const;

	SubMeshDescriptor neighborSubMesh; // Info about the adjacent SubMesh, including references to flow variables.
};

//////////////////////////////////////
// Mesh edge BCs above ↑
//////////////////////////////////////
// Immersed BCs below ↓
//////////////////////////////////////

// Base class for immersed boundaries.
// Suggested derived types: Adiabatic wall and iso-thermal wall
// or different shapes.
class ImmersedBoundary
{
public:

	ImmersedBoundary(const SubMeshDescriptor& subMeshData) :
		subMeshData(subMeshData)
	{}

	virtual ~ImmersedBoundary() = default;

	virtual std::unique_ptr<ImmersedBoundary> getUniquePtrToCopy() = 0;

	virtual void identifyRelatedNodes(const ConfigSettings& params,
   	   	    	  	  	  	  	  	  const Vector3_d& gridSpacing,
									  const Vector3_i& nMeshNodes,
									  const Vector3_d& meshOriginOffset,
									  Array3D<NodeTypeEnum>& nodeTypeArray	// ← OUTPUT
			  	  	  	  	  	  	  ) = 0; // ← PURE virtual

	void applyBoundaryCondition(const ConfigSettings& params);

	virtual IntegralProperties getIntegralProperties(const ConfigSettings& params) = 0;

protected:

	vector<GhostNode> ghostNodes;		// The solid ghost nodes adjacent to this surface
	std::map<int, int> ghostNodeMap;	// A codex from 1D index in the mesh, to index of corresponding entry in ghostNodes
	SubMeshDescriptor subMeshData;		// Info and flow variables for the submesh with this immersed surface.

	void findGhostNodesWithFluidNeighbors(const vector<int>& solidNodeIndices);

	void checkIfSurroundingShouldBeGhost(const Vector3_i &surroundingNode,
			   	   	   	   	   	   	     vector<GhostNode>& newGhostNodes);

	vector<GhostNode> setImagePointPositions(GhostNodeVectorIterator firstGhostToProcess);

	GhostNodeVectorIterator appendGhostNodes(const vector<GhostNode>& newGhostNodes);

	void populateVandermondeDirichlet(const InterpolationPositions& interpolationPoints,
									  Matrix8x8_d& vandermonde);

	void populateVandermondeNeumann(const InterpolationPositions& interpolationPoints,
									const Array8_b& ghostFlags,
									const vector<Vector3_d>& unitNormals,
									Matrix8x8_d& vandermonde);

	void setInterpolationValues(
			const IndexBoundingBox& surroundingNodes,
			InterpolationValues& interpolationValues,		// ↰
			InterpolationPositions& interpolationPositions,	// ↰
			Array8_b& ghostFlag,							// ← Output
			bool& allSurroundingAreFluid,					// ↲
			vector<Vector3_d>& unitNormals);				// ↲

	double trilinearInterpolation(const Vector8_d& interpolationValues,
													 const Vector3_d& imagePoint,
													 const Matrix8x8_d& vandermonde);

private:
	virtual Vector3_d getNormalProbe(const Vector3_d& ghostNodePosition) = 0;

	double simplifiedInterpolation(const Vector8_d& interpolationValues,
			 	 	 	 	 	   const Vector3_i& lowerIndexNode,
								   const Vector3_d& imagePointPosition);

	PrimitiveVariablesScalars simplifiedInterpolationAll(
			const InterpolationValues& interpolationValues,
			const Vector3_i& lowerIndexNode,
			const Vector3_d& imagePointPosition );

	PrimitiveVariablesScalars getGhostNodeBCVariables(const PrimitiveVariablesScalars& imagePointBCVars);

	PrimitiveVariablesScalars trilinearInterpolationAll(const InterpolationValues&,
														const Vector3_d& imagePoint,
														const Matrix8x8_d& vandermondeDirichlet,
														const Matrix8x8_d& vandermondeNeumann);

	void setInterpolationValuesFluidNode(int counter, int surroundingNodeIndex1D,
										 InterpolationValues& interpolationValues,			// ← OUTPUT
										 InterpolationPositions& interpolationPositions);	// ← OUTPUT

	void setInterpolationValuesGhostNode(
			int counter,
			int surroundingNodeIndex1D,
			vector<Vector3_d>& unitNormals,					// ←
			InterpolationValues& interpolationValues,		// ← Output
			InterpolationPositions& interpolationPositions);// ←

	PrimitiveVariablesScalars interpolateImagePointVariables(
			const InterpolationValues& interpolationValues,
			const InterpolationPositions& interpolationPositions,
			bool allSurroundingAreFluid,
			const Array8_b& ghostFlag,
			const GhostNode& ghostNode,
			const IndexBoundingBox& surroundingNodes,
			const vector<Vector3_d>& unitNormals );
};

// Class to define boundary conditions at an immersed cylinder:
class CylinderBody : public ImmersedBoundary
{
public:
	CylinderBody(Vector3_d centroidPosition, AxisOrientationEnum axis, double radius, const SubMeshDescriptor& subMeshData);

	std::unique_ptr<ImmersedBoundary> getUniquePtrToCopy() override
	{ return std::make_unique<CylinderBody>(*this); }

	void identifyRelatedNodes(const ConfigSettings& params,
	   	   	   	   	    	  const Vector3_d& gridSpacing,
							  const Vector3_i& nMeshNodes,
							  const Vector3_d& meshOriginOffset,
							  Array3D<NodeTypeEnum>& nodeTypeArray	// ← Output
							  ) override;

	IntegralProperties getIntegralProperties(const ConfigSettings& params) override;

private:
	Vector3_d centroidPosition;	// Only 2 coordinates used, depending on axis orientation
	AxisOrientationEnum axis;	// Orientation of centroid axis
	double radius;

	Vector3_d getNormalProbe(const Vector3_d& ghostNodePosition) override;

	IndexBoundingBox getCylinderBoundingBox(const Vector3_d& gridSpacing,
			  	  	  	  	  	  	  	  	const Vector3_i& nMeshNodes) const;

	void getSolidNodesInCylinder(const ConfigSettings& params,
			   	   	   	   	   	 	 	  const IndexBoundingBox& indicesToCheck,
										  vector<int>& solidNodeIndices // ← Output
			   	   	   	   	   	 	 	  );
};

// Class to define boundary conditions at an immersed sphere:
class SphereBody : public ImmersedBoundary
{
public:
	SphereBody(Vector3_d centerPosition, double radius, const SubMeshDescriptor& subMeshData);

	std::unique_ptr<ImmersedBoundary> getUniquePtrToCopy() override
	{ return std::make_unique<SphereBody>(*this); }

	void identifyRelatedNodes(const ConfigSettings& params,
   	   	    	  	  	  	  const Vector3_d& gridSpacing,
							  const Vector3_i& nMeshNodes,
							  const Vector3_d& meshOriginOffset,
							  Array3D<NodeTypeEnum>& nodeTypeArray	// ← Output
			  	  	  	  	  ) override;

	IntegralProperties getIntegralProperties(const ConfigSettings& params) override;

private:
	Vector3_d centerPosition;
	double radius;

	Vector3_d getNormalProbe(const Vector3_d& ghostNodePosition) override;

	IndexBoundingBox getSphereBoundingBox(const Vector3_d& gridSpacing, int filterNodeLayerWidth) const;

	void getSolidNodesInSphere(const ConfigSettings& params,
	   	   	   	   	   	   	   const IndexBoundingBox& indicesToCheck,
							   vector<int>& solidNodeIndices // ← Output
	   	   	   	   	   	   	   );
};

#endif /* SRC_BOUNDARY_H_ */




