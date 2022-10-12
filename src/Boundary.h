/*
 * BoundaryCondition.h
 *
 *  Created on: Oct 12, 2022
 *      Author: frederk
 */

#ifndef SRC_BOUNDARY_H_
#define SRC_BOUNDARY_H_

#include "includes_and_names.h"

enum class BoundaryTypeEnum
{
	domainEdge, immersed
};

enum class DomainBoundaryNormalAxisEnum
{
	x, y, z
};

enum class BoundaryConditionTypeEnum
{
	inlet, outlet, adiabaticWall, isothermalWall, symmetry
};

class Boundary
{
protected:
	Boundary(BoundaryTypeEnum boundaryType, BoundaryConditionTypeEnum bcType) :
		boundaryType{boundaryType}, bcType{bcType} {}
public:
	virtual ~Boundary();
	virtual void applyBoundaryCondition();
	const BoundaryTypeEnum boundaryType;
	BoundaryConditionTypeEnum bcType;
};

class DomainBoundary : public Boundary
{
public:
	DomainBoundary(BoundaryTypeEnum boundaryType,
				   BoundaryConditionTypeEnum bcType,
				   DomainBoundaryNormalAxisEnum normalAxis,
				   uint planeIndex) :
					   Boundary(boundaryType, bcType),
					   normalAxis{normalAxis},
					   planeIndex{planeIndex}
	{}
	const DomainBoundaryNormalAxisEnum normalAxis;
	const uint planeIndex;
};

class ImmersedBoundary : public Boundary
{

};

#endif /* SRC_BOUNDARY_H_ */
