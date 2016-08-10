/**
 * Copyright 2016 Ehab Abdelhamid, Ibrahim Abdelaziz

GraphX.h

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
 * A representation of a graph, consisting of a set of nodes, and each node has a set of edges
 */
#ifndef GRAPHX_H_
#define GRAPHX_H_

#include <tr1/unordered_map>
#include <tr1/unordered_set>
#include<set>
#include"NodeX.h"
#include "Settings.h"

using namespace std;

class GraphX
{
private:
	int id;
	int type;//0-undirected, 1-directed, default value is 0. IMPORTANT: current algorithms do not support directed graphs
	tr1::unordered_map<int, NodeX*> nodes;
	int numOfEdges;
	tr1::unordered_map<double, set<int>* > nodesByLabel;
	string CL;//canonical label
	long freq = -1;//

public:
	//initializers
	GraphX();
	GraphX(int, int );
	void init(int, int );
	GraphX(GraphX* );
	~GraphX();

	//load graph data
	bool parseData(istream& data ,tr1::unordered_map<string, void* >& );
	bool loadFromFile(string fileName ,tr1::unordered_map<string, void* >& );
	bool loadFromString(string data ,tr1::unordered_map<string, void* >& );

	//add new node
	NodeX* AddNode(int id, double label);

	//adding an edge
	void addEdge(int , int , double);
	void addEdge(NodeX* , NodeX* , double);
	void addEdge(int , int , double ,tr1::unordered_map<string, void* >& );

	//remove an edge
	void removeEdge(int , int );

	//remove a node and discard its incident edges
	void removeNode_IgnoreEdges(int );

	//iteratot functions
	tr1::unordered_map<int,NodeX*>::const_iterator getNodesIterator();
	tr1::unordered_map<int,NodeX*>::const_iterator getNodesEndIterator();

	int getID() {return id;}
	int getType() {return type;}
	double getEdgeLabel(int srcNode, int destNode);
	int getNumOfNodes();
	NodeX* getNodeWithID(int nodeID);
	set<int>* getNodesByLabel(double );
	friend ostream& operator<<(ostream& os, const GraphX& g);
	int getNumOfEdges() {return numOfEdges;}

	//check if graph is connected
	bool isConnected();

	//this subgraph isomorphism code is correct
	void isIsomorphic(GraphX* graph, vector<map<int, int>* >& , tr1::unordered_map<int, tr1::unordered_set<int>*>& ,int restrictedDomainID = -1, int restrictedNodeID = -1, bool pruneByDomainValues = true, tr1::unordered_map<int, tr1::unordered_set<int>* >* postponedNodes = 0, unsigned long maxIters = Settings::postponeNodesAfterIterations);

	//graph isomorphism function
	bool isTheSameWith(GraphX* otherG);

	//get a unique canonical label of this graph
	string getCanonicalLabel();

	void resetCL() {CL.clear();}
	unsigned long numIterations;
	bool CanSatisfyNodeLabels(GraphX* otherG);
	void setFrequency(long freq) { this->freq = freq; }
};

ostream& operator<<(ostream& , const GraphX& );

class Set_Iterator
{
public:
	int static st_SI;
	Set_Iterator() {st_SI++;}
	~Set_Iterator() {st_SI--;}
	set<int>* se;//list of candidate nodes in the graph to match with
	set<int>::iterator it;//iterator over elements in se
	bool isIterEnd() { return it==se->end(); }
};

#endif /* GRAPHX_H_ */
