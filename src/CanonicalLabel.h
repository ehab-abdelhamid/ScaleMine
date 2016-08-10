/**
 * Copyright 2016 Ehab Abdelhamid, Ibrahim Abdelaziz

CanonicalLabel.h

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
 * This class is responsible for generating a unique graph canonical label
 */
#ifndef CANONICALLABEL_H_
#define CANONICALLABEL_H_

#include<iostream>
#include<map>
#include<vector>
#include<string>
#include"GraphX.h"

using namespace std;

bool CL_tester();

class NodeWithInfo
{
public:
	NodeX* node;
	string nl;//neighbors list
	int partID;

	NodeWithInfo(NodeX* node);
	~NodeWithInfo()
	{
		NodeWithInfo::st_counter--;
	}
	string toString();
	int compareTo(NodeWithInfo& nInfo);
	bool operator < (const NodeWithInfo& ) const;

	friend bool NodeWithInfo_ascending (NodeWithInfo* a, NodeWithInfo* b);
	friend bool NodeWithInfo_descending (NodeWithInfo* a, NodeWithInfo* b);

	static int st_counter;
};

bool NodeWithInfo_ascending (NodeWithInfo* a, NodeWithInfo* b);
bool NodeWithInfo_descending (NodeWithInfo* a, NodeWithInfo* b);

class CL_Partition
{
private:
	int partID;
	vector<NodeWithInfo*> nodes;
	vector<vector<NodeWithInfo* >* > combinations;

	//the below used for sorting, every partition must have all nodes with the same values for the below variables
	int degree;
	string nodeLabel;
	string nl;
	void combinations_fn(vector<NodeWithInfo* >* notused, bool);

public:
	static int st_counter;
	int counter;
	bool operator< (const CL_Partition& str) const;
	void setSortingValues();
	CL_Partition(int partID);
	void addNode(NodeWithInfo* nInfo);
	vector<NodeWithInfo*>::const_iterator getNodesEnum() const;
	vector<NodeWithInfo*>::const_iterator getNodesEnd() const {return nodes.end();}
	int getNumNodes();
	void clearNodes();
	void setID(int partID);
	int getID();
	map<string, CL_Partition*>* getNewParts();
	string toString();
	vector<vector<NodeWithInfo* >* >* getCombinations();
	void sortNodes();

	void CL_Part_Destructor();

	friend bool ascending (CL_Partition* a,CL_Partition* b);
	friend bool descending (CL_Partition* a,CL_Partition* b);

	string getPrefix();
};

bool ascending (CL_Partition* a,CL_Partition* b);
bool descending (CL_Partition* a,CL_Partition* b);

class PartID_label
{
public:
	int static st_counter;
	PartID_label() {PartID_label::st_counter++;}
	~PartID_label() {PartID_label::st_counter--;}
	int partID;
	string label;

	string toString();
};

class CanonicalLabel
{
private:
	static const bool enablePrint = false;

	static bool sortPartitions(vector<CL_Partition* >& parts);
	static char* generateCanLabel(vector<NodeWithInfo* >* nodes, GraphX* graph, double**, bool oneLineOnly = false);

public:
	static long elapsed1;
	static long elapsed2;
	static char sigBuffer[10000];

	static string generate(GraphX* graph);
	static string generateNeighborsList(NodeWithInfo* nInfo, map<int , NodeWithInfo* >& allNodes);
};

#endif /* CANONICALLABEL_H_ */
