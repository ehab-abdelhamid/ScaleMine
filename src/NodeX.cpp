/**
 * Copyright 2016 Ehab Abdelhamid, Ibrahim Abdelaziz

NodeX.cpp

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
 * Represent a graph node
 */

#include<iostream>
#include <stdlib.h>
#include<tr1/unordered_map>
#include "EdgeX.h"
#include "NodeX.h"

NodeX::NodeX(int id, double label)
{
	this->id = id;
	this->label = label;
}

NodeX::~NodeX()
{
	for( tr1::unordered_map<int, void*>::const_iterator ii=edges.begin(); ii!=edges.end(); ++ii)
    {
    	EdgeX* edge = (EdgeX*)((*ii).second);
    	delete edge;
    }
	edges.clear();

	for( tr1::unordered_map<int, void*>::const_iterator ii=revEdges.begin(); ii!=revEdges.end(); ++ii)
	{
	  	EdgeX* edge = (EdgeX*)((*ii).second);
	   	delete edge;
	}
	revEdges.clear();
}

/**
 * Add an edge to this node
 * Parameters: the other node, and the edge label
 */
void NodeX::addEdge(NodeX* otherNode, double edgeLabel, int graphType)
{
	if(edges.find(otherNode->getID())!=edges.end())
		return;
	EdgeX* edge = new EdgeX(edgeLabel, otherNode);

	edges[otherNode->getID()] = edge;

	if(graphType==1)
	{
		EdgeX* edge = new EdgeX(edgeLabel, this);
		otherNode->revEdges[this->getID()] = edge;
	}
}

void NodeX::removeEdge(NodeX* otherNode, int graphType)
{
	tr1::unordered_map<int, void*>::iterator iter = edges.find(otherNode->getID());
	if(iter!=edges.end())
		delete iter->second;
	edges.erase(otherNode->getID());

	if(graphType==1)
	{
		//do somethingh here
		tr1::unordered_map<int, void*>::iterator iter = otherNode->revEdges.find(this->getID());
		if(iter!=revEdges.end())
			delete iter->second;
		revEdges.erase(this->getID());
	}
}

void* NodeX::getEdgeForDestNode(int destNodeID )
{
	tr1::unordered_map<int, void*>::iterator temp = edges.find(destNodeID);
	if(temp==edges.end())
		return NULL;
	else
		return (*temp).second;
}

bool NodeX::isItConnectedWithNodeID(int nodeID)
{
	//check node connectivity
	tr1::unordered_map<int, void*>::iterator iter = edges.find(nodeID);
	if(iter==edges.end())
		return false;

	return true;
}

bool NodeX::isItConnectedWithNodeID(int nodeID, double label)
{
	//check node connectivity
	tr1::unordered_map<int, void*>::iterator iter = edges.find(nodeID);
	if(iter==edges.end())
		return false;

	//check edge label
	if(((EdgeX*)iter->second)->getLabel()!=label)
		return false;

	return true;
}

/**
 * check whether the 'node' parameter in neighborhood consistent with 'this' node
 */
bool NodeX::isNeighborhoodConsistent(NodeX* node)
{
	tr1::unordered_map<double, int> labels;
	//populate labels of this node
	for(tr1::unordered_map<int, void*>::iterator iter = edges.begin();iter!=edges.end();iter++)
	{
		double otherNodeLabel = ((EdgeX*)iter->second)->getLabel();
		tr1::unordered_map<double, int>::iterator tempIter = labels.find(otherNodeLabel);
		if(tempIter==labels.end())
			labels.insert(std::pair<double, int>(otherNodeLabel, 1));
		else
		{
			int currentCount = tempIter->second;
			labels.erase(otherNodeLabel);
			labels.insert(std::pair<double, int>(otherNodeLabel, currentCount+1));
		}
	}

	//check labels against this's labels
	for(tr1::unordered_map<int, void*>::iterator iter = node->getEdgesIterator();iter!=node->getEdgesEndIterator();iter++)
	{
		double otherNodeLabel = ((EdgeX*)iter->second)->getLabel();
		tr1::unordered_map<double, int>::iterator tempIter = labels.find(otherNodeLabel);
		if(tempIter==labels.end())
			return false;
		int currentCount = tempIter->second;
		labels.erase(otherNodeLabel);
		if(currentCount>1)
			labels.insert(std::pair<double, int>(otherNodeLabel, currentCount-1));
	}

	return true;
}

ostream& operator<<(ostream& os, const NodeX& n)
{
    os << n.id << '[' << n.label << "]" <<endl;
    for( tr1::unordered_map<int, void*>::const_iterator ii=n.edges.begin(); ii!=n.edges.end(); ++ii)
    {
    	EdgeX* edge = (EdgeX*)((*ii).second);
    	os<<"--"<<edge->getLabel()<<"-->"<<edge->getOtherNode()->getID()<<endl;
    }
    return os;
}
