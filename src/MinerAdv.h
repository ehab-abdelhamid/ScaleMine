/**
 * Copyright 2016 Ehab Abdelhamid, Ibrahim Abdelaziz

MinerAdv.h

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

#ifndef MINERADV_H_
#define MINERADV_H_

#include "Miner.h"

class MinerAdv : public Miner
{
public:
	static long maxNumCandids;

	unsigned long numOfVisitedNodes;
	unsigned long numIterations;

	int nThreads;

	void startMining(string ,int ,int ,int ,int );
	bool popACandidateApprox(vector<CLMap*>& , map<int, Pattern*>& , int destination);
	void sendACandidateApprox(string , Pattern* , map<int, Pattern*>& , int );

	void setInputGraph(GraphX* g) { this->graph = g; }
	GraphX* getInputGraph() { return graph; }
	void printTotalExpectedTime();
};

#endif /* MINERADV_H_ */
