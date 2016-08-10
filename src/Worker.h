/**
 * Copyright 2016 Ehab Abdelhamid, Ibrahim Abdelaziz

Worker.h

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
 * represents a compute node
 */
#ifndef WORKER_H_
#define WORKER_H_

#include "Pattern.h";


#define BUFFERSIZE 2000

class Worker
{
protected:
	GraphX* graph;
	int frequency;

public:
	int machineID;
	int nThreads;
	void start(int machineID, int nThreads);
	static void *processing(Worker* worker, int id);
	bool isFrequent(Pattern* candidate, int subTaskNum, int subTaskMax, int* mniTable, tr1::unordered_map<int, tr1::unordered_set<int>* >* , unsigned long maxIters );
	int getApproxFreq(Pattern* candidate, int& invalidCol, int& getApproxFreq, bool&, unsigned long& , unsigned long& , unsigned long& );
	int getFrequency() { return frequency; }
};



#endif /* WORKER_H_ */
