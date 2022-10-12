/*
 * BoundaryCondition.h
 *
 *  Created on: Oct 12, 2022
 *      Author: frederk
 */

#ifndef SRC_BOUNDARY_H_
#define SRC_BOUNDARY_H_

#include "includes_and_names.h"

enum class DomainBoundaryNormalAxisEnum
{
	x, y, z
};

enum class BoundaryConditionTypeEnum
{
	inlet, outlet, adiabaticWall, isothermalWall, periodic, symmetry
};

class MeshEdgeBoundary
{
public:
	MeshEdgeBoundary(BoundaryConditionTypeEnum bcType,
				   DomainBoundaryNormalAxisEnum normalAxis,
				   uint planeIndex) :
					   normalAxis{normalAxis},
					   planeIndex{planeIndex}
	{}
	virtual ~MeshEdgeBoundary();
	virtual void applyBoundaryCondition();
	const DomainBoundaryNormalAxisEnum normalAxis;
	const uint planeIndex;
};

class ImmersedBoundary
{
public:
	ImmersedBoundary(BoundaryConditionTypeEnum bcType)
	{}
	virtual void applyBoundaryCondition();
	virtual ~ImmersedBoundary();
};

#endif /* SRC_BOUNDARY_H_ */




