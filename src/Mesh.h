/*
 * Mesh.h
 *
 *  Created on: Oct 7, 2022
 *      Author: frederk
 */

#ifndef SRC_MESH_H_
#define SRC_MESH_H_

class Mesh;

#include "Array3D.h"
#include "FlowVariableGroupStructs.h"
#include "includes_and_names.h"
#include "ConfigSettings.h"
#include "Boundary.h"

typedef vector<unique_ptr<MeshEdgeBoundary>> EdgeBoundaryCollection;
typedef vector<unique_ptr<ImmersedBoundary>> ImmersedBoundaryCollection;

class Mesh
{
public:

	Mesh(const ConfigSettings& params);

	void setupBoundaries(const ConfigSettings& params);

	void assertBoundaryConditionCompliance();

	void categorizeNodes(const ConfigSettings& params);

	void applyFilter_ifAppropriate(Array3D_d& variable_old, Array3D_d& variable_new,
									uint filterInterval, ulong timeLevel);

	ConservedVariablesScalars computeNorms_conservedVariables();

	void swapConservedVariables();

	Vector3_u getIndices3D(size_t index1D) const;

	size_t getIndex1D(size_t i, size_t j, size_t k) const;

	size_t getIndex1D(const Vector3_u& indices) const
	{ return getIndex1D(indices.i, indices.j, indices.k); }

	Vector3_d getNodePosition(size_t i, size_t j, size_t k) const;

	Vector3_d getNodePosition(const Vector3_u& indices) const
	{ return getNodePosition(indices.i, indices.j, indices.k); }

	IndexBoundingBox getSurroundingNodesBox(Vector3_d point) const;

	void applyAllBoundaryConditions(double t, const ConfigSettings& params);

	void setFlowVariablesAtNode(size_t index1D, const ConservedVariablesScalars&  conservedVariableScalars,
											  const PrimitiveVariablesScalars&  primitiveVariableScalars,
											  const TransportPropertiesScalars& transportPropertyScalars)
	{
		conservedVariables.rho  (index1D) = conservedVariableScalars.rho;
		conservedVariables.rho_u(index1D) = conservedVariableScalars.rho_u;
		conservedVariables.rho_v(index1D) = conservedVariableScalars.rho_v;
		conservedVariables.rho_w(index1D) = conservedVariableScalars.rho_w;
		conservedVariables.rho_E(index1D) = conservedVariableScalars.rho_E;
		primitiveVariables.u(index1D) = primitiveVariableScalars.u;
		primitiveVariables.v(index1D) = primitiveVariableScalars.v;
		primitiveVariables.w(index1D) = primitiveVariableScalars.w;
		primitiveVariables.p(index1D) = primitiveVariableScalars.p;
		primitiveVariables.T(index1D) = primitiveVariableScalars.T;
		transportProperties.mu   (index1D) = transportPropertyScalars.mu;
		transportProperties.kappa(index1D) = transportPropertyScalars.kappa;
	}

	void setFlowVariablesAtNode(Vector3_u node, const ConservedVariablesScalars&  conservedVariableScalars,
												const PrimitiveVariablesScalars&  primitiveVariableScalars,
												const TransportPropertiesScalars& transportPropertyScalars)
	{
		conservedVariables.rho  (node) = conservedVariableScalars.rho;
		conservedVariables.rho_u(node) = conservedVariableScalars.rho_u;
		conservedVariables.rho_v(node) = conservedVariableScalars.rho_v;
		conservedVariables.rho_w(node) = conservedVariableScalars.rho_w;
		conservedVariables.rho_E(node) = conservedVariableScalars.rho_E;
		primitiveVariables.u(node) = primitiveVariableScalars.u;
		primitiveVariables.v(node) = primitiveVariableScalars.v;
		primitiveVariables.w(node) = primitiveVariableScalars.w;
		primitiveVariables.p(node) = primitiveVariableScalars.p;
		primitiveVariables.T(node) = primitiveVariableScalars.T;
		transportProperties.mu   (node) = transportPropertyScalars.mu;
		transportProperties.kappa(node) = transportPropertyScalars.kappa;
	}

	const size_t NI, NJ, NK;	// Mesh size. Number of nodes in x,y,z directions
	const size_t nNodesTotal;	// Total number of nodes in the mesh
	const double dx, dy, dz;	// Grid spacing in x-, y- and z-direction
	ConservedVariablesArrayGroup conservedVariables;	// Mass density, momentum & total energy
	ConservedVariablesArrayGroup conservedVariablesOld;	// Previous time level. Or temporary storage, when needed.
	PrimitiveVariablesArrayGroup primitiveVariables;	// Velocity, pressure & Temperature
	TransportPropertiesArrayGroup transportProperties;	// Viscosity and thermal conductivity
	RK4slopesArrayGroup RK4slopes;						// 4 slopes for each conserved variable
	vector<size_t> activeNodeIndices;	// Indices to active fluid nodes
	Array3D_nodeType nodeType;			// Type/category of each node (ghost, fluid, etc.)
private:
	EdgeBoundaryCollection edgeBoundaries;			// Boundaries at the edges of the Cartesian mesh
	ImmersedBoundaryCollection immersedBoundaries;	// Boundaries at immersed bodies

	double getNormOfChange(const Array3D_d& oldValue, const Array3D_d& newValue);
};

#endif /* SRC_MESH_H_ */
