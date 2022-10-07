/*
 * Mesh.cpp
 *
 *  Created on: Oct 7, 2022
 *      Author: frederk
 */

#include "Mesh.h"

Mesh::Mesh(uint nMeshNodesX, uint nMeshNodesY, uint nMeshNodesZ,
		double domainLengthX, double domainLengthY, double domainLengthZ ) :
NI{nMeshNodesX}, NJ{nMeshNodesY}, NK{nMeshNodesZ},
rho     (NI, NJ, NK),
rho_u   (NI, NJ, NK),
rho_v   (NI, NJ, NK),
rho_w   (NI, NJ, NK),
E       (NI, NJ, NK),
u       (NI, NJ, NK),
v       (NI, NJ, NK),
w       (NI, NJ, NK),
p       (NI, NJ, NK),
T       (NI, NJ, NK),
mu      (NI, NJ, NK),
kappa   (NI, NJ, NK),
k1_rho  (NI, NJ, NK),
k2_rho  (NI, NJ, NK),
k3_rho  (NI, NJ, NK),
k4_rho  (NI, NJ, NK),
k1_rho_u(NI, NJ, NK),
k2_rho_u(NI, NJ, NK),
k3_rho_u(NI, NJ, NK),
k4_rho_u(NI, NJ, NK),
k1_rho_v(NI, NJ, NK),
k2_rho_v(NI, NJ, NK),
k3_rho_v(NI, NJ, NK),
k4_rho_v(NI, NJ, NK),
k1_rho_w(NI, NJ, NK),
k2_rho_w(NI, NJ, NK),
k3_rho_w(NI, NJ, NK),
k4_rho_w(NI, NJ, NK),
k1_E    (NI, NJ, NK),
k2_E    (NI, NJ, NK),
k3_E    (NI, NJ, NK),
k4_E    (NI, NJ, NK),
interm_rho  (NI, NJ, NK),
interm_rho_u(NI, NJ, NK),
interm_rho_v(NI, NJ, NK),
interm_rho_w(NI, NJ, NK),
interm_E    (NI, NJ, NK)
{
	setGridSpacings(domainLengthX, domainLengthY, domainLengthZ);
}


// Calculate the space between nodes in the grid, based on domain size and no. of nodes.
void Mesh::setGridSpacings(double domainLengthX,
							 double domainLengthY,
							 double domainLengthZ )
{
	dx = domainLengthX / (NI - 1);
	dy = domainLengthY / (NJ - 1);
	dz = domainLengthZ / (NK - 1);
	cout << "Grid spacings set: dx = " << dx << " , dy = " << dy << " , dz = " << dz << endl;
}

// Filters one variable field, i.e., one solution array, 'filterVariable' and stores the filtered
// result in 'variableTemporaryStorage'. Then, the arrays are swapped by move-semantics.
// Only filters if the modulo of time level plus one, by the filter interval is zero.
// This causes the first filtering to happen as late as possible.
void Mesh::applyFilter_ifAppropriate(Array3D_d& filterVariable, Array3D_d& variableTemporaryStorage,
									uint filterInterval, uint timeLevel)
{
	if(filterInterval > 0)
		if( (timeLevel+1) % filterInterval == 0 )
		{
			uint iMax{NI-1}, jMax{NJ-1}, kMax{NK-1};
			// Copy boundary nodes:
			for(uint i{0}; i<=iMax; ++i)
				for(uint j{0}; j<=jMax; ++j)
				{
					variableTemporaryStorage(i,j,0   ) = filterVariable(i,j,0   );
					variableTemporaryStorage(i,j,kMax) = filterVariable(i,j,kMax);
				}
			for(uint i{0}; i<=iMax; ++i)
				for(uint k{1}; k<=kMax-1; ++k)
				{
					variableTemporaryStorage(i,0,   k) = filterVariable(i,0,   k);
					variableTemporaryStorage(i,jMax,k) = filterVariable(i,jMax,k);
				}
			for(uint j{1}; j<=jMax-1; ++j)
				for(uint k{1}; k<=kMax-1; ++k)
				{
					variableTemporaryStorage(0,   j,k) = filterVariable(0,   j,k);
					variableTemporaryStorage(iMax,j,k) = filterVariable(iMax,j,k);
				}
			// Apply filter to interior nodes:
			for(uint i{1}; i<=iMax-1; ++i)
				for(uint j{1}; j<=jMax-1; ++j)
					for(uint k{1}; k<=kMax-1; ++k)
					{
						variableTemporaryStorage(i,j,k) = 1./2.  *   filterVariable(i,j,k)
														+ 1./12. * ( filterVariable(i+1,j,k) + filterVariable(i-1,j,k)
																   + filterVariable(i,j+1,k) + filterVariable(i,j-1,k)
																   + filterVariable(i,j,k+1) + filterVariable(i,j,k-1) );
					}
			filterVariable.dataSwap(variableTemporaryStorage);	// Swap the arrays using move-sematics (super-fast)
		}
}

// Compute the norm of change in the conserved variables, and store in the history vectors.
void Mesh::computeNorms_conservedVariables()
{
	normHistory_rho  .push_back( getNormOfChange(rho  , interm_rho  ) );
	normHistory_rho_u.push_back( getNormOfChange(rho_u, interm_rho_u) );
	normHistory_rho_v.push_back( getNormOfChange(rho_v, interm_rho_v) );
	normHistory_rho_w.push_back( getNormOfChange(rho_w, interm_rho_w) );
	normHistory_E    .push_back( getNormOfChange(E    , interm_E    ) );
}

// Swap the contents of all the arrays of conserved variables and the intermediate arrays, by move-semantics.
// This operation is super fast and needs no extra copy. Only the ownership of the data is changed.
void Mesh::swapConservedVariables()
{
	rho  .dataSwap(interm_rho  );
	rho_u.dataSwap(interm_rho_u);
	rho_v.dataSwap(interm_rho_v);
	rho_w.dataSwap(interm_rho_w);
	E    .dataSwap(interm_E    );
}

// Compute the 2-norm of the difference between two arrys. Intended to monitor the change between two consecutive time levels.
// E.g. to check convergence of solution.
double Mesh::getNormOfChange(const Array3D_d& oldValue, const Array3D_d& newValue)
{
	double sumOfSquaredChanges = 0;
	uint numberOfMeshNodes = NI * NJ * NK;
	for(uint i=0; i<numberOfMeshNodes; ++i)
		sumOfSquaredChanges += pow(oldValue(i)-newValue(i), 2);
	double normOfChange = sqrt( dx*dy*dz * sumOfSquaredChanges );
	return normOfChange;
}





