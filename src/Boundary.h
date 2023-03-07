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
	min, max
};

typedef Eigen::Array<double, 8, 1> Vector8_d;
typedef Eigen::Matrix<double, 8, 8> Matrix8x8_d;

// Package struct with data from the mesh.
// Shortens argument lists, and lets Boundary not depend on Mesh
struct MeshDescriptor
{
	const Vector3_i& nNodes;
	const Vector3_d& spacings;
	const Vector3_d& originOffset;
	const Array3D_nodeType& nodeType;
	AllFlowVariablesArrayGroup& flowVariables;

	MeshDescriptor(const Vector3_i& nNodes,
				   const Vector3_d& gridSpacings,
				   const Vector3_d& originOffset,
				   Array3D_nodeType& nodeTypeArray,
				   AllFlowVariablesArrayGroup& flowVariables)
	: nNodes{nNodes},
	  spacings{gridSpacings},
	  originOffset{originOffset},
	  nodeType{nodeTypeArray},
	  flowVariables{flowVariables}
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


// Base class for boundaries at the the edge of the Cartesian computational mesh.
// Suggested derived types: Inlet, outlet, walls, periodic, symmetry...
class MeshEdgeBoundary
{
public:
	MeshEdgeBoundary(AxisOrientationEnum normalAxis,
					 EdgeIndexEnum planeIndex,
					 NodeTypeEnum ownedNodesType);

	virtual ~MeshEdgeBoundary() = default;

	void identifyOwnedNodes(IndexBoundingBox& unclaimedNodes,
							const Vector3_i& nMeshNodes,
							Array3D_nodeType& nodeTypeArray );

	virtual void applyBoundaryCondition(double t, const Vector3_i& nMeshNodes,		// <- Input
										const ConfigSettings& params,				// <- Input
										AllFlowVariablesArrayGroup& flowVariables)	// <- Output
										= 0; // <- PURE virtual

	const AxisOrientationEnum normalAxis;	// The axis that is normal to the boundary plane
	const EdgeIndexEnum planeIndex;			// Denotes if the plane is at lowest or highest index side
	const NodeTypeEnum ownedNodesType;		// The node type that is assigned to the owned nodes
	IndexBoundingBox ownedNodes;			// The nodes that the BC at this boundary can affect
protected:

	void getAdjacentIndices(int index1D, const Vector3_i& nMeshNodes,	// <- Input
			  	  	  	  	int& boundaryAdjacentIndex, 					// <- Output
							int& nextToAdjacentIndex);					// <- Output

	int getPeriodicIndex(int index1D, const Vector3_i& nMeshNodes);
};

// Class to define inlet boundary condition:
class InletBoundary : public MeshEdgeBoundary
{
public:
	InletBoundary(AxisOrientationEnum normalAxis,
				  EdgeIndexEnum planeIndex,
				  double velocity);

	void filterInletDensity(const Vector3_i& nMeshNodes,				// <- Input
	   	   	   	   	   	    const ConfigSettings& params,				// <- Input
							AllFlowVariablesArrayGroup& flowVariables);	// <- Output

	void applyBoundaryCondition(double t, const Vector3_i& nMeshNodes,		// <- Input
								const ConfigSettings& params,				// <- Input
								AllFlowVariablesArrayGroup& flowVariables)	// <- Output
								override;

private:
	double velocity;	// The prescribed inlet velocity, that we reach after initial ramp-up.
};

// Class to define outlet boundary condition:
class OutletBoundary : public MeshEdgeBoundary
{
public:
	OutletBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex);

	void applyBoundaryCondition(double t, const Vector3_i& nMeshNodes,		// <- Input
								const ConfigSettings& params,				// <- Input
								AllFlowVariablesArrayGroup& flowVariables)	// <- Output
								override;
};

// Class to define a periodic boundary condition:
class PeriodicBoundary : public MeshEdgeBoundary
{
public:
	PeriodicBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex);

	void applyBoundaryCondition(double t, const Vector3_i& nMeshNodes,		// <- Input
								const ConfigSettings& params,				// <- Input
								AllFlowVariablesArrayGroup& flowVariables)	// <- Output
								override;
};

// Class to define a symmetry boundary condition:
class SymmetryBoundary : public MeshEdgeBoundary
{
public:
	SymmetryBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex);

	void applyBoundaryCondition(double t, const Vector3_i& nMeshNodes,		// <- Input
								const ConfigSettings& params,				// <- Input
								AllFlowVariablesArrayGroup& flowVariables)	// <- Output
								override;
};

// Class to define a boundary condition where all variables are extrapolated:
class ExtrapolationBoundary : public MeshEdgeBoundary
{
public:
	ExtrapolationBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex);

	void applyBoundaryCondition(double t, const Vector3_i& nMeshNodes,		// <- Input
								const ConfigSettings& params,				// <- Input
								AllFlowVariablesArrayGroup& flowVariables)	// <- Output
								override;
};

// Base class for immersed boundaries.
// Suggested derived types: Adiabatic wall and iso-thermal wall
// or different shapes.
class ImmersedBoundary
{
public:
	ImmersedBoundary();

	virtual ~ImmersedBoundary() = default;

	virtual void identifyRelatedNodes(const ConfigSettings& params,
   	   	    	  	  	  	  	  	  const Vector3_d& gridSpacing,
									  const Vector3_i& nMeshNodes,
									  const Vector3_d& meshOriginOffset,
									  Array3D_nodeType& nodeTypeArray	// <- OUTPUT
			  	  	  	  	  	  	  ) = 0; // <- PURE virtual

	void applyBoundaryCondition(const ConfigSettings& params,
								MeshDescriptor& mesh // <- In/Output
			  	  	  	  	    );

	virtual IntegralProperties getIntegralProperties(const ConfigSettings& params, const MeshDescriptor& mesh) = 0;

protected:

	vector<GhostNode> ghostNodes;			// The solid ghost nodes adjacent to this surface
	std::map<int, int> ghostNodeMap;	// A codex from 1D index in the mesh, to index of corresponding entry in ghostNodes
	vector<int> filterNodes;				// The closest interior fluid nodes. Filtered when BC is applied.

	void findGhostNodesWithFluidNeighbors(const vector<int>& solidNodeIndices,
										  const Vector3_i& nMeshNodes,
										  Array3D_nodeType& nodeTypeArray);

	void checkIfSurroundingShouldBeGhost(const Vector3_i &surroundingNode,
			   	   	   	   	   	   	     vector<GhostNode>& newGhostNodes,
										 Array3D_nodeType& nodeTypeArray);

	vector<GhostNode> setImagePointPositions(GhostNodeVectorIterator firstGhostToProcess,
			   	   	   	   	   	   	   	   	 const Vector3_d& gridSpacing,
											 const Vector3_i& nMeshNodes,
											 const Vector3_d& meshOriginOffset,
											 Array3D_nodeType& nodeTypeArray);

	GhostNodeVectorIterator appendGhostNodes(const vector<GhostNode>& newGhostNodes, const Vector3_i& nMeshNodes);

	void populateVandermondeDirichlet(const InterpolationPositions& interpolationPoints,
									  Matrix8x8_d& vandermonde);

	void populateVandermondeNeumann(const InterpolationPositions& interpolationPoints,
									const Array8_b& ghostFlags,
									const vector<Vector3_d>& unitNormals,
									Matrix8x8_d& vandermonde);

	void setInterpolationValues(
			const IndexBoundingBox& surroundingNodes,
			const MeshDescriptor& mesh,
			InterpolationValues& interpolationValues,		// <-
			InterpolationPositions& interpolationPositions,	// <-
			Array8_b& ghostFlag,							// <- Output
			bool& allSurroundingAreFluid,					// <-
			vector<Vector3_d>& unitNormals);				// <-

	double trilinearInterpolation(const Vector8_d& interpolationValues,
													 const Vector3_d& imagePoint,
													 const Matrix8x8_d& vandermonde);

private:
	virtual Vector3_d getNormalProbe(const Vector3_d& ghostNodePosition) = 0;

	void filterClosestFluidNodes(const Vector3_i& nMeshNodes,
								 const AllFlowVariablesArrayGroup& flowVariables);

	double simplifiedInterpolation(const Vector8_d& interpolationValues,
			 	 	 	 	 	   const Vector3_i& lowerIndexNode,
								   const Vector3_d& imagePointPosition,
								   const Vector3_d& gridSpacing,
								   const Vector3_d& meshOriginOffset);

	PrimitiveVariablesScalars simplifiedInterpolationAll(
			const InterpolationValues& interpolationValues,
			const Vector3_i& lowerIndexNode,
			const Vector3_d& imagePointPosition,
			const Vector3_d& gridSpacing,
			const Vector3_d& meshOriginOffset );

	PrimitiveVariablesScalars getGhostNodeBCVariables(const PrimitiveVariablesScalars& imagePointBCVars);

	PrimitiveVariablesScalars trilinearInterpolationAll(const InterpolationValues&,
														const Vector3_d& imagePoint,
														const Matrix8x8_d& vandermondeDirichlet,
														const Matrix8x8_d& vandermondeNeumann);

	void setInterpolationValuesFluidNode(int counter, int surroundingNodeIndex1D,
			   	   	   	   	   	   	     const MeshDescriptor& mesh,
										 InterpolationValues& interpolationValues,			// <- OUTPUT
										 InterpolationPositions& interpolationPositions);	// <- OUTPUT

	void setInterpolationValuesGhostNode(
			int counter,
			int surroundingNodeIndex1D,
			vector<Vector3_d>& unitNormals,					// <-
			InterpolationValues& interpolationValues,		// <- Output
			InterpolationPositions& interpolationPositions);// <-

	PrimitiveVariablesScalars interpolateImagePointVariables(
			const InterpolationValues& interpolationValues,
			const InterpolationPositions& interpolationPositions,
			bool allSurroundingAreFluid,
			const Array8_b& ghostFlag,
			const GhostNode& ghostNode,
			const IndexBoundingBox& surroundingNodes,
			const vector<Vector3_d>& unitNormals,
			const MeshDescriptor& mesh);
};

// Class to define boundary conditions at an immersed cylinder:
class CylinderBody : public ImmersedBoundary
{
public:
	CylinderBody(Vector3_d centroidPosition, AxisOrientationEnum axis, double radius);

	void identifyRelatedNodes(const ConfigSettings& params,
	   	   	   	   	    	  const Vector3_d& gridSpacing,
							  const Vector3_i& nMeshNodes,
							  const Vector3_d& meshOriginOffset,
							  Array3D_nodeType& nodeTypeArray	// <- Output
							  ) override;

	IntegralProperties getIntegralProperties(const ConfigSettings& params, const MeshDescriptor& mesh) override;

private:
	Vector3_d centroidPosition;	// Only 2 coordinates used, depending on axis orientation
	AxisOrientationEnum axis;	// Orientation of centroid axis
	double radius;

	Vector3_d getNormalProbe(const Vector3_d& ghostNodePosition) override;

	IndexBoundingBox getCylinderBoundingBox(const Vector3_d& gridSpacing,
			  	  	  	  	  	  	  	  	const Vector3_i& nMeshNodes,
											int filterNodesLayerWidth) const;

	void getSolidAndFilterNodesInCylinder(const ConfigSettings& params,
			   	   	   	   	   	 	 	  const IndexBoundingBox& indicesToCheck,
										  const Vector3_d& gridSpacing,
										  const Vector3_i& nMeshNodes,
										  const Vector3_d& meshOriginOffset,
										  int nNodesFilterLayer,
										  vector<int>& solidNodeIndices, // <- Output
										  Array3D_nodeType& nodeTypeArray   // <- Output
			   	   	   	   	   	 	 	  );
};

// Class to define boundary conditions at an immersed sphere:
class SphereBody : public ImmersedBoundary
{
public:
	SphereBody(Vector3_d centerPosition, double radius);

	void identifyRelatedNodes(const ConfigSettings& params,
   	   	    	  	  	  	  const Vector3_d& gridSpacing,
							  const Vector3_i& nMeshNodes,
							  const Vector3_d& meshOriginOffset,
							  Array3D_nodeType& nodeTypeArray	// <- Output
			  	  	  	  	  ) override;

	IntegralProperties getIntegralProperties(const ConfigSettings& params, const MeshDescriptor& mesh) override;

private:
	Vector3_d centerPosition;
	double radius;

	Vector3_d getNormalProbe(const Vector3_d& ghostNodePosition) override;

	IndexBoundingBox getSphereBoundingBox(const Vector3_d& gridSpacing, int filterNodeLayerWidth) const;

	void getSolidAndFilterNodesInSphere(const ConfigSettings& params,
	   	   	   	   	   	   	   	   	    const IndexBoundingBox& indicesToCheck,
										const Vector3_d& gridSpacing,
										const Vector3_i& nMeshNodes,
										const Vector3_d& meshOriginOffset,
										int nNodesFilterLayer,
										vector<int>& solidNodeIndices, // <- Output
										Array3D_nodeType& nodeTypeArray   // <- Output
	   	   	   	   	   	   	   	   	    );
};

#endif /* SRC_BOUNDARY_H_ */




