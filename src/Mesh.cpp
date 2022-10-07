/*
 * Mesh.cpp
 *
 *  Created on: Oct 7, 2022
 *      Author: frederk
 */

#include "Mesh.h"

Mesh::Mesh()
{

}

// Calculate the space between nodes in the grid, based on domain size and no. of nodes.
void Solver::setGridSpacings()
{
	dx = params.L_x / (params.NI - 1);
	dy = params.L_y / (params.NJ - 1);
	dz = params.L_z / (params.NK - 1);
	cout << "Grid spacings set: dx = " << dx << " , dy = " << dy << " , dz = " << dz << endl;
}

// Filters one variable field, i.e., one solution array, 'filterVariable' and stores the filtered
// result in 'variableTemporaryStorage'. Then, the arrays are swapped by move-semantics.
// Only filters if the modulo of time level plus one, by the filter interval is zero.
// This causes the first filtering to happen as late as possible.
void Solver::applyFilter_ifAppropriate(Array3D_d& filterVariable, Array3D_d& variableTemporaryStorage)
{
	if(params.filterInterval > 0)
		if( (timeLevel+1) % params.filterInterval == 0 )
		{
			uint iMax{params.NI-1}, jMax{params.NJ-1}, kMax{params.NK-1};
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
void Solver::computeNorms_conservedVariables()
{
	normHistory_rho  .push_back( getNormOfChange(rho  , interm_rho  ) );
	normHistory_rho_u.push_back( getNormOfChange(rho_u, interm_rho_u) );
	normHistory_rho_v.push_back( getNormOfChange(rho_v, interm_rho_v) );
	normHistory_rho_w.push_back( getNormOfChange(rho_w, interm_rho_w) );
	normHistory_E    .push_back( getNormOfChange(E    , interm_E    ) );
}

// Swap the contents of all the arrays of conserved variables and the intermediate arrays, by move-semantics.
// This operation is super fast and needs no extra copy. Only the ownership of the data is changed.
void Solver::swapConservedVariables()
{
	rho  .dataSwap(interm_rho  );
	rho_u.dataSwap(interm_rho_u);
	rho_v.dataSwap(interm_rho_v);
	rho_w.dataSwap(interm_rho_w);
	E    .dataSwap(interm_E    );
}

// Compute the 2-norm of the difference between two arrys. Intended to monitor the change between two consecutive time levels.
// E.g. to check convergence of solution.
double Solver::getNormOfChange(const Array3D_d& oldValue, const Array3D_d& newValue)
{
	double sumOfSquaredChanges = 0;
	uint numberOfMeshNodes = params.NI * params.NJ * params.NK;
	for(uint i=0; i<numberOfMeshNodes; ++i)
		sumOfSquaredChanges += pow(oldValue(i)-newValue(i), 2);
	double normOfChange = sqrt( dx*dy*dz * sumOfSquaredChanges );
	return normOfChange;
}





