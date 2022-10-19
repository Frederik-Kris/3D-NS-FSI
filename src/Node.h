/*
 * Node.h
 *
 *  Created on: Oct 19, 2022
 *      Author: frederk
 */

#ifndef SRC_NODE_H_
#define SRC_NODE_H_

#include "includes_and_names.h"

struct GhostNode
{
	Vector3_u indices;
	Vector3_d bodyInterceptPoint;
	Vector3_d imagePoint;

	GhostNode(Vector3_u indices) : indices(indices) {}
};


#endif /* SRC_NODE_H_ */
