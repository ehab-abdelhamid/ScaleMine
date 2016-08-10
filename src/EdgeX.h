/**
 * Copyright 2016 Ehab Abdelhamid, Ibrahim Abdelaziz

EdgeX.h

This file is part of ScaleMine.

ScaleMine is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

ScaleMine is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with Grami.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * representing graph edge
 */
#ifndef EDGEX_H_
#define EDGEX_H_

#include "NodeX.h"

class EdgeX
{
private:
	double label;	//edge label
	NodeX* otherNode;	//edges are originating from one node to the other node, this referes to the other node

public:
	EdgeX(double label, NodeX* otherNode);
	double getLabel() {return label;}
	NodeX* getOtherNode(){return otherNode;}
};

#endif /* EDGEX_H_ */
