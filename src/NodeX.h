/**
 * Copyright 2016 Ehab Abdelhamid, Ibrahim Abdelaziz

NodeX.h

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

#ifndef NODEX_H_
#define NODEX_H_

#include <map>
#include<vector>
#include<ostream>
#include <tr1/unordered_map>

using namespace std;

class NodeX
{
private:
	int id;
	double label;
	tr1::unordered_map<int, void*> edges;
	tr1::unordered_map<int, void*> revEdges;

public:
	NodeX(int id, double value);
	~NodeX();
	void addEdge(NodeX* , double edgeLabel, int graphType);
	void removeEdge(NodeX* , int graphType);
	friend ostream& operator<<(ostream& os, const NodeX& n);
	int getID() {return id;}
	double getLabel() {return label;}
	tr1::unordered_map<int, void*>::iterator getEdgesIterator() {return edges.begin();}
	tr1::unordered_map<int, void*>::iterator getEdgesEndIterator() {return edges.end();}
	int getEdgesSize() {return edges.size(); }
	void* getEdgeForDestNode(int destNodeID );
	bool isItConnectedWithNodeID(int nodeID);
	bool isItConnectedWithNodeID(int nodeID, double label);
	bool isNeighborhoodConsistent(NodeX* );
};

#endif /* NODEX_H_ */
