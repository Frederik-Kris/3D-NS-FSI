/*
 * BoundaryCondition.h
 *
 *  Created on: Oct 12, 2022
 *      Author: frederk
 */

#ifndef SRC_BOUNDARY_H_
#define SRC_BOUNDARY_H_

class MeshEdgeBoundary;
class ImmersedBoundary;

#include "includes_and_names.h"
#include "Array3D.h"
#include "ConfigSettings.h"
#include "FlowVariableGroupStructs.h"
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

typedef Eigen::Array<double, 8, 1> Vector8_d;
typedef Eigen::Matrix<double, 8, 8> Matrix8x8_d;

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
					 EdgeIndexEnum planeIndex,
					 NodeTypeEnum ownedNodesType);

	virtual ~MeshEdgeBoundary() = default;

	void identifyOwnedNodes(IndexBoundingBox& unclaimedNodes,
							const Vector3_u& nMeshNodes,
							Array3D_nodeType& nodeTypeArray );

	virtual void applyBoundaryCondition(double t, const Vector3_u& nMeshNodes,		// <- Input
										const ConfigSettings& params,				// <- Input
										AllFlowVariablesArrayGroup& flowVariables)	// <- Output
										= 0; // <- PURE virtual

	const AxisOrientationEnum normalAxis;
	const EdgeIndexEnum planeIndex;
	const NodeTypeEnum ownedNodesType;
protected:
	vector<size_t> nodeIndices;

	void getAdjacentIndices(size_t index1D, const Vector3_u& nMeshNodes,	// <- Input
			  	  	  	  	size_t& boundaryAdjacentIndex, 					// <- Output
							size_t& nextToAdjacentIndex);					// <- Output

	size_t getPeriodicIndex(size_t index1D, const Vector3_u& nMeshNodes);
};

// Class to define inlet boundary condition:
class InletBoundary : public MeshEdgeBoundary
{
public:
	InletBoundary(AxisOrientationEnum normalAxis,
				  EdgeIndexEnum planeIndex,
				  double velocity);

	void applyBoundaryCondition(double t, const Vector3_u& nMeshNodes,		// <- Input
								const ConfigSettings& params,				// <- Input
								AllFlowVariablesArrayGroup& flowVariables)	// <- Output
								override;

private:
	double velocity;
};

// Class to define outlet boundary condition:
class OutletBoundary : public MeshEdgeBoundary
{
public:
	OutletBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex);

	void applyBoundaryCondition(double t, const Vector3_u& nMeshNodes,		// <- Input
								const ConfigSettings& params,				// <- Input
								AllFlowVariablesArrayGroup& flowVariables)	// <- Output
								override;
};

// Class to define a periodic boundary condition:
class PeriodicBoundary : public MeshEdgeBoundary
{
public:
	PeriodicBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex);

	void applyBoundaryCondition(double t, const Vector3_u& nMeshNodes,		// <- Input
								const ConfigSettings& params,				// <- Input
								AllFlowVariablesArrayGroup& flowVariables)	// <- Output
								override;
};

// Class to define a symmetry boundary condition:
class SymmetryBoundary : public MeshEdgeBoundary
{
public:
	SymmetryBoundary(AxisOrientationEnum normalAxis, EdgeIndexEnum planeIndex);

	void applyBoundaryCondition(double t, const Vector3_u& nMeshNodes,		// <- Input
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
									  const Vector3_u& nMeshNodes,
									  const Vector3_d& meshOriginOffset,
									  Array3D_nodeType& nodeTypeArray	// <- Output
			  	  	  	  	  	  	  ) = 0; // <- PURE virtual

	void applyBoundaryCondition(const Vector3_u& nMeshNodes,
			  	  	  	  	  	const Vector3_d& gridSpacing,
								const Vector3_d& meshOriginOffset,
								const ConfigSettings& params,
								const Array3D_nodeType& nodeTypeArray,
								AllFlowVariablesArrayGroup& flowVariables // <- Output
			  	  	  	  	    );

protected:

	vector<GhostNode> ghostNodes;
	std::map<size_t, GhostNode*> ghostNodeMap;

	void findGhostNodesWithFluidNeighbors(const vector<size_t>& solidNodeIndices,
										  const Vector3_u& nMeshNodes,
										  Array3D_nodeType& nodeTypeArray);

	void checkIfSurroundingShouldBeGhost(const Vector3_u &surroundingNode,
			   	   	   	   	   	   	     vector<GhostNode>& newGhostNodes,
										 Array3D_nodeType& nodeTypeArray);

	vector<GhostNode> setImagePointPositions(GhostNodeVectorIterator firstGhostToProcess,
			   	   	   	   	   	   	   	   	 const Vector3_d& gridSpacing,
											 const Vector3_u& nMeshNodes,
											 const Vector3_d& meshOriginOffset,
											 Array3D_nodeType& nodeTypeArray);

	GhostNodeVectorIterator appendGhostNodes(const vector<GhostNode>& newGhostNodes, const Vector3_u& nMeshNodes);

private:
	virtual Vector3_d getNormalProbe(const Vector3_d& ghostNodePosition) = 0;

	double simplifiedInterpolation(const Vector8_d& interpolationValues,
			 	 	 	 	 	   const Vector3_u& lowerIndexNode,
								   const Vector3_d& imagePointPosition,
								   const Vector3_d& gridSpacing,
								   const Vector3_d& meshOriginOffset);

	PrimitiveVariablesScalars simplifiedInterpolationAll(
			const InterpolationValues& interpolationValues,
			const Vector3_u& lowerIndexNode,
			const Vector3_d& imagePointPosition,
			const Vector3_d& gridSpacing,
			const Vector3_d& meshOriginOffset );

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
			   	   	   	   	   	   	     const Vector3_u& nMeshNodes, const Vector3_d& gridSpacing,
										 const Vector3_d& meshOriginOffset,
										 const AllFlowVariablesArrayGroup &flowVariables,
										 InterpolationValues& interpolationValues,			// <- OUTPUT
										 InterpolationPositions& interpolationPositions);	// <- OUTPUT

	void setInterpolationValuesGhostNode(
			uint counter,
			size_t surroundingNodeIndex1D,
			vector<Vector3_d>& unitNormals,					// <-
			InterpolationValues& interpolationValues,		// <- Output
			InterpolationPositions& interpolationPositions);// <-

	void setInterpolationValues(
			const IndexBoundingBox& surroundingNodes,
			const Vector3_u& nMeshNodes,
			const Vector3_d& gridSpacing,
			const Vector3_d& meshOriginOffset,
			const Array3D_nodeType& nodeTypeArray,
			const AllFlowVariablesArrayGroup& flowVariables,
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
			const Vector3_d& gridSpacing,
			const Vector3_d& meshOriginOffset);
};

// Class to define boundary conditions at an immersed cylinder:
class CylinderBody : public ImmersedBoundary
{
public:
	CylinderBody(Vector3_d centroidPosition, AxisOrientationEnum axis, double radius);

	void identifyRelatedNodes(const ConfigSettings& params,
	   	   	   	   	    	  const Vector3_d& gridSpacing,
							  const Vector3_u& nMeshNodes,
							  const Vector3_d& meshOriginOffset,
							  Array3D_nodeType& nodeTypeArray	// <- Output
							  ) override;

private:
	Vector3_d centroidPosition;
	AxisOrientationEnum axis;
	double radius;

	Vector3_d getNormalProbe(const Vector3_d& ghostNodePosition) override;

	IndexBoundingBox getCylinderBoundingBox(const Vector3_d& gridSpacing,
			  	  	  	  	  	  	  	  	const Vector3_u& nMeshNodes) const;

	void getSolidNodesInCylinder(const ConfigSettings& params,
			   	   	   	   	   	 const IndexBoundingBox& indicesToCheck,
								 const Vector3_d& gridSpacing,
								 const Vector3_u& nMeshNodes,
								 const Vector3_d& meshOriginOffset,
								 vector<size_t>& solidNodeIndices, // <- Output
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
							  const Vector3_u& nMeshNodes,
							  const Vector3_d& meshOriginOffset,
							  Array3D_nodeType& nodeTypeArray	// <- Output
			  	  	  	  	  ) override;

private:
	Vector3_d centerPosition;
	double radius;

	Vector3_d getNormalProbe(const Vector3_d& ghostNodePosition) override;

	IndexBoundingBox getSphereBoundingBox(const Vector3_d& gridSpacing) const;

	void getSolidNodesInSphere(const ConfigSettings& params,
	   	   	   	   	   	   	   const IndexBoundingBox& indicesToCheck,
							   const Vector3_d& gridSpacing,
							   const Vector3_u& nMeshNodes,
							   const Vector3_d& meshOriginOffset,
							   vector<size_t>& solidNodeIndices, // <- Output
							   Array3D_nodeType& nodeTypeArray   // <- Output
	   	   	   	   	   	   	   );
};

#endif /* SRC_BOUNDARY_H_ */




