/**
 * Copyright 2016 Ehab Abdelhamid, Ibrahim Abdelaziz

GraMiCounter.h

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

#ifndef GRAMICOUNTER_H_
#define GRAMICOUNTER_H_

#include <tr1/unordered_set>
#include <vector>
#include "GraphX.h"
#include "Pattern.h"

class pair_with_edge
{
public:
	int id1;
	int id2;
	double edgeLabel;
	int minDomainSize;
};

class GraMiCounter
{
private:
	static void setDomainsOrder(tr1::unordered_map<int, tr1::unordered_set<int>*>& domains_values, int* order, Pattern* pattern);
public:
	static int numByPassedNodes;
	static bool useAC3;
	static int numSideEffectNodes;
	static int numPostponedNodes;

	static int isFrequent(GraphX*, Pattern*, int support, double approximate);
	static int isFrequent_adv(GraphX*, Pattern*, int support, int subTaskNum, int subTaskMax, int* mniTable, tr1::unordered_map<int, tr1::unordered_set<int>* >* postponedNodes = 0, unsigned long maxIters = Settings::postponeNodesAfterIterations);
	static int isFrequent_approx(GraphX*, Pattern*, int support, int& invalidCol, int& predictedValids, bool allowEarlyBreak, bool& exact, unsigned long& numVisiteNodes, unsigned long& numIterations, vector<unsigned long>& listOfNumOfIters, unsigned long& predictedTime);
	static void AC_3(GraphX*, tr1::unordered_map<int, tr1::unordered_set<int>*>& , Pattern* , int support);
	static void AC_3(GraphX*, tr1::unordered_map<int, tr1::unordered_set<int>*>& , GraphX* , int support, int invalidCol = -1);
	static bool refine(GraphX*, tr1::unordered_map<int, tr1::unordered_set<int>*>& , tr1::unordered_set<int>* , tr1::unordered_set<int>* , pair_with_edge* , int );
	static bool isItAcyclic(GraphX& );
};

#endif /* GRAMICOUNTER_H_ */
