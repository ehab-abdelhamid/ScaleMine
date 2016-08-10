/**
 * Copyright 2016 Ehab Abdelhamid, Ibrahim Abdelaziz

Pattern.h

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

#ifndef PATTERN_H_
#define PATTERN_H_

#include <set>
#include<list>
#include<tr1/unordered_set>
#include "GraphX.h"
#include "EdgeX.h"

using namespace std;

/**
 * represents a primary graph of a specific subgraph
 * primary graphs are used for join-based subgraph extension
 */
class PrimaryGraph
{
private:
	GraphX* graph;
	int srcNodeID;//source from where the edge was deleted
	EdgeX* edge;//the edge that was deleted
	bool b;//false until a value is set for this primary graph

public:
	PrimaryGraph() {b=false;}
	~PrimaryGraph() {delete graph;}
	void setValues(GraphX* graph, int srcNodeID, EdgeX* edge) {this->graph = graph;this->srcNodeID = srcNodeID; this->edge = edge; b = true;}
	GraphX* getGraph() { return graph; }
	int getSrcNodeID() { return srcNodeID; }
	EdgeX* getEdge() { return edge; }
	//string getCanLabel() { return CL; }
	bool isTheSameWith(PrimaryGraph* );
	bool isTheSameWith(GraphX* otherPG);
	bool isSet() { return b; }
};

/**
 * represents a frequent subgraph
 */
class Pattern
{
private:
	int ID;
	int frequency; //initialized to -1 to indicate that it is not set
	GraphX* graph = 0;
	bool graphCopied;
	vector<tr1::unordered_set<int>*> occurences;
	unsigned long predictedTime;
	int invalidCol = -1;
	int predictedValids = -1;
	int predictedFreq = -1;

	int subtasking = 1;//how many tasks (workers) will be given to finish this pattern
	int subtaskingFixed = 1;
	void createEffeGraphsList();
	void freeEffeGraphs();
	int* tempMNITable;
	int remainingSubtasks;
	bool resultExact;//set to exact if the approximate frequency function can return exact result
	int* postponedNodes_mniTable;
	unsigned long maxIters;

public:
	bool selected = false;
	int static st_counter;
	static int maxPatternID;
	bool postponeExpensiveNodes;
	bool inTheProcessOfPrediction = false;
	Pattern* predictedPattern = 0;

	Pattern(GraphX* ,bool copyGraph = true);
	Pattern(Pattern* );
	void init();
	~Pattern();
	int getFrequency();
	int getID();
	void setFrequency(int );
	void addNode(int nodeID, int patternNodeID);
	GraphX* getGraph() {return graph; }
	string toString();
	void Combine(Pattern* ,int addToID=0);//the two patterns should have a similar graph
	vector<tr1::unordered_set<int>*>* getOccurences() {return &occurences;}
	void invalidateFrequency() {frequency=-1;}
	void extend(int srcID, int destID, double destLabel, double edgeLabel);//apply one edge extension
	bool doINeedToCountInThisPartition(int );
	void setPartitionFrequency(int freq, int partID);
	int getSize();//get pattern size as the number of edges

	list<PrimaryGraph*> primaryGraphs;
	set<std::pair<PrimaryGraph*, PrimaryGraph* > > getJoiningPG(Pattern* pattern);
	void generatePrimaryGraphs();
	bool hasUniqueLabels();
	void makeIDNegative();
	void setPredictedTime(unsigned long pt) {this->predictedTime = pt;}
	void setInvalidCol(int invC, int pValids);
	void borrowTimeInfor(Pattern*, int nWorkers);
	unsigned long getPredictedTime() {return this->predictedTime;}
	int getInvalidCol() { return invalidCol; }
	int getPredictedValids() { return predictedValids; }
	void setSubtasking(int st, int nWorkers);
	int getSubtaskingValue() { return this->subtasking; }
	int getSubtaskingValueFixed() { return this->subtaskingFixed; }
	void setSubtaskTaken() { this->subtasking--; }
	bool subtasksFinished() { return subtasking==0; }
	int subTaskDone(int* mniTable, int* postponed_mniTable, int support);
	void resetMNI();
	void setResultExact() { resultExact = true; }
	bool isResultExact() { return resultExact; }
	void setPredictedFreq(int pf) { this->predictedFreq = pf; }
	int getPredictedFreq() { return this->predictedFreq; }
	void setMaxIters(unsigned long mni) { maxIters = mni; }
	unsigned long getMaxIters() { return maxIters; }
};

#endif /* PATTERN_H_ */
