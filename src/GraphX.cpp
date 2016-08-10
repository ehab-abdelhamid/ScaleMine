/**
 * Copyright 2016 Ehab Abdelhamid, Ibrahim Abdelaziz

GraphX.cpp

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

#include <tr1/unordered_set>
#include <tr1/unordered_map>
#include<iostream>
#include<fstream>
#include<sstream>
#include <stdlib.h>
#include <unistd.h>
#include <limits>
#include "GraphX.h"
#include "EdgeX.h"
#include "utils.h"
#include "Pattern.h"
#include "GraMiCounter.h"
#include "Settings.h"

bool checkCurrentStatus(GraphX* graph1, GraphX* graph2, vector<Set_Iterator* >& currentS, int* nodesOrder, tr1::unordered_map<int, int>& selectedQueryMap)
{
	NodeX* dataGraphNode = graph2->getNodeWithID(*(currentS[currentS.size()-1]->it));
	NodeX* queryGraphNode = graph1->getNodeWithID(nodesOrder[currentS.size()-1]);

	for(tr1::unordered_map<int, void*>::iterator iter = queryGraphNode->getEdgesIterator();iter!=queryGraphNode->getEdgesEndIterator();++iter)
	{
		int otherNodeID = iter->first;
		EdgeX* edge = (EdgeX*)(iter->second);

		if(selectedQueryMap.find(otherNodeID)==selectedQueryMap.end())
			continue;

		//get the data node that the current data node should be connected with
		int dataNodeID = *(currentS[selectedQueryMap.find(otherNodeID)->second]->it);

		//check the current graph edges, whether it has a connection to dataNodeID or not
		if(!dataGraphNode->isItConnectedWithNodeID(dataNodeID, edge->getLabel()))
		{
			//cout<<"false"<<endl;
			return false;
		}
	}

	return true;
}

void insertIntoEdgeFreq(NodeX* src, int destNodeID, double destNodeLabel, double edgeLabel, tr1::unordered_map<string, void* >& edgeToFreq, bool addSrcOnly)
{
	string key;
	if(src->getLabel()>destNodeLabel)
	{
		stringstream sstmGF;
		if(src->getLabel()==destNodeLabel)
			sstmGF << src->getLabel()<<destNodeLabel<<","<<edgeLabel<<",";
		else
			sstmGF << src->getLabel()<<","<<destNodeLabel<<","<<edgeLabel<<",";
		key = sstmGF.str();
	}
	else
	{
		stringstream sstmGF;
		if(destNodeLabel==src->getLabel())
			sstmGF <<destNodeLabel<<src->getLabel()<<","<<edgeLabel<<",";
		else
			sstmGF <<destNodeLabel<<","<<src->getLabel()<<","<<edgeLabel<<",";
		key = sstmGF.str();
	}

	tr1::unordered_map<string, void* >::iterator iter = edgeToFreq.find(key);
	Pattern* pattern;
	if(iter==edgeToFreq.end())
	{
		GraphX* graph = new GraphX();
		NodeX* node1 = graph->AddNode(0, src->getLabel());
		NodeX* node2 = graph->AddNode(1, destNodeLabel);
		graph->addEdge(node1, node2, edgeLabel);
		pattern = new Pattern(graph);
		delete graph;

		edgeToFreq[key] = pattern;
	}
	else
	{
		pattern = (Pattern*)((*iter).second);
	}

	if(pattern->getGraph()->getNodeWithID(0)->getLabel()==src->getLabel())
	{
		pattern->addNode(src->getID(), 0);
		if(!addSrcOnly) pattern->addNode(destNodeID, 1);
	}
	else
	{
		pattern->addNode(src->getID(), 1);
		if(!addSrcOnly) pattern->addNode(destNodeID, 0);
	}
}

/**
 * Constructor
 */
GraphX::GraphX()
{
	init(0,0);
}

GraphX::GraphX(int id, int type)
{
	init(id, type);
}

void GraphX::init(int id, int type)
{
	this->id = id;
	this->type = type;
	if(type==1)
	{
		cout<<"IMPORTANT! current set of algorithms are not tested on directed graphs!";
		exit(0);
	}
	numOfEdges = 0;
	this->freq = -1;
}

GraphX::GraphX(GraphX* graph)
{
	this->id = graph->getID();
	this->type = graph->getType();
	numOfEdges = graph->getNumOfEdges();

	//copy nodes
	for(tr1::unordered_map<int, NodeX*>::const_iterator iter = graph->getNodesIterator();iter!=graph->getNodesEndIterator();++iter)
	{
		NodeX* oldNode = iter->second;
		this->AddNode(oldNode->getID(), oldNode->getLabel());
	}

	//copy edges
	for(tr1::unordered_map<int, NodeX*>::const_iterator iter = graph->getNodesIterator();iter!=graph->getNodesEndIterator();++iter)
	{
		NodeX* oldNode = iter->second;
		NodeX* node = nodes.find(oldNode->getID())->second;

		for(tr1::unordered_map<int, void*>::iterator iter2 = oldNode->getEdgesIterator();iter2!=oldNode->getEdgesEndIterator();++iter2)
		{
			EdgeX* edge = (EdgeX*)(iter2->second);
			node->addEdge(nodes.find(edge->getOtherNode()->getID())->second, edge->getLabel(), this->type);
		}
	}
}

/**
 * Add a node to the graph
 * Parameters are: node id, and node label
 */
NodeX* GraphX::AddNode(int id, double label)
{
	CL.clear();

	tr1::unordered_map<int, NodeX*>::iterator temp = nodes.find(id);
	if(temp!=nodes.end())
		return temp->second;

	NodeX* node = new NodeX(id, label);

	nodes.insert(std::pair<int, NodeX*>(id, node));

	//add to the 'nodes by label' map
	tr1::unordered_map<double, set<int>* >::iterator iter = nodesByLabel.find(label);
	if(iter==nodesByLabel.end())
	{
		nodesByLabel.insert(std::pair<double, set<int>*>(label, new set<int>()));
		iter = nodesByLabel.find(label);
	}
	iter->second->insert(id);

	return node;
}

bool GraphX::CanSatisfyNodeLabels(GraphX* otherG)
{
	for(tr1::unordered_map<double, set<int>* >::iterator iter = otherG->nodesByLabel.begin();iter!=otherG->nodesByLabel.end();iter++)
	{
		double label = (*iter).first;
		int count = (*iter).second->size();
		if(count==0) continue;

		tr1::unordered_map<double, set<int>* >::iterator thisIter = nodesByLabel.find(label);
		if(thisIter==nodesByLabel.end())
			return false;
		int thisCount = (*iter).second->size();
		if(thisCount<count)
			return false;
	}

	return true;
}

/**
 * Add an edge between two nodes in the graph, the nodes must exist before adding the edge
 * Parameters: the source node ID, the destination node ID, and the edge label
 * For undirected graphs, one more edge will be added in the reverse direction
 */
void GraphX::addEdge(NodeX* srcNode, NodeX* destNode, double edgeLabel)
{
	CL.clear();

	srcNode->addEdge(destNode, edgeLabel, this->type);

	if(this->type==0)
	{
		destNode->addEdge(srcNode, edgeLabel, this->type);
	}

	numOfEdges++;
}

/**
 * Add an edge to this graph given source node, destination node, and an edge label
 */
void GraphX::addEdge(int srcID, int destID, double edgeLabel)
{
	CL.clear();

	NodeX* srcNode;
	NodeX* destNode;

	tr1::unordered_map<int, NodeX*>::iterator iter = nodes.find(srcID);
	if(iter!=nodes.end())
		srcNode = nodes[srcID];
	else
		return;

	iter = nodes.find(destID);
	if(iter!=nodes.end())
		destNode = nodes[destID];
	else
		return;

	srcNode->addEdge(destNode, edgeLabel, this->type);

	if(this->type==0)
	{
		destNode->addEdge(srcNode, edgeLabel, this->type);
	}

	numOfEdges++;

	return;
}

/**
 * remove an edge by its incident node ids
 */
void GraphX::removeEdge(int id1, int id2)
{
	CL.clear();

	nodes[id1]->removeEdge(nodes[id2], this->type);
	if(nodes[id1]->getEdgesSize()==0)
		removeNode_IgnoreEdges(id1);

	if(this->type==0)
	{
		nodes[id2]->removeEdge(nodes[id1], this->type);
		if(nodes[id2]->getEdgesSize()==0)
			removeNode_IgnoreEdges(id2);
	}

	numOfEdges--;
}

/**
 * remove a node from the graph, ignoring edges (as a prepost all edges connecting to this node should be already removed)
 */
void GraphX::removeNode_IgnoreEdges(int nodeID)
{
	CL.clear();

	nodesByLabel.find(getNodeWithID(nodeID)->getLabel())->second->erase(nodeID);
	delete nodes.find(nodeID)->second;
	nodes.erase(nodeID);
}

void GraphX::addEdge(int srcID, int destID, double edgeLabel, tr1::unordered_map<string, void* >& edgeToFreq)
{
	if(CL.length()>0)
		CL.clear();

	NodeX* srcNode;
	NodeX* destNode;

	tr1::unordered_map<int, NodeX*>::iterator iter = nodes.find(srcID);
	if(iter!=nodes.end())
		srcNode = iter->second;
	else
		return;

	iter = nodes.find(destID);
	if(iter!=nodes.end())
		destNode = iter->second;
	else
		return;

	this->addEdge(srcNode, destNode, edgeLabel);

	insertIntoEdgeFreq(srcNode, destID, destNode->getLabel(), edgeLabel, edgeToFreq, false);

	if(this->type==0)
	{
		insertIntoEdgeFreq(destNode, srcID, srcNode->getLabel(), edgeLabel, edgeToFreq, false);
	}
}

/**
 * Load a graph file that has .lg format
 * return true if loading is done correctly, otherwise false
 */
bool GraphX::loadFromFile(string fileName, tr1::unordered_map<string, void* >& edgeToFreq)
{
	CL.clear();

	cout<<"Loading graph from file: "<<fileName<<endl;
	ifstream file (fileName.c_str(), ios::in);
	if(!file)
	{
		cout << "While opening a file an error is encountered" << endl;
		return false;
    }

	if(!parseData(file, edgeToFreq))
		return false;

	file.close();

	return true;
}

/**
 * load a graph from the given string. string should follow the .lg format
 */
bool GraphX::loadFromString(string data, tr1::unordered_map<string, void* >& edgeToFreq)
{
	CL.clear();

	istringstream str(data);

	bool b = parseData(str, edgeToFreq);

	//destruct data in the 'edgeToFreq' structure
	for(tr1::unordered_map<string, void* >::iterator iter = edgeToFreq.begin();iter!=edgeToFreq.end();iter++)
	{
		delete ((Pattern*)iter->second);
	}
	edgeToFreq.clear();

	if(!b)
		return false;

	return true;
}

/**
 * the graph loader parser
 */
bool GraphX::parseData(istream& data, tr1::unordered_map<string, void* >& edgeToFreq)
{
	//read the first line
	char temp_ch;
	data>>temp_ch;data>>temp_ch;data>>temp_ch;

	int numEdgesLoaded = 0;

	bool firstEdgeMet = false;
	while (!data.eof())
	{
		char ch;
		ch = '\0';
		data>>ch;

		//to add nodes
		if(ch=='v')
		{
			int id;
			double label;
			data>>id;
			data>>label;

			this->AddNode(id, label);
		}
		else if(ch=='e')//to add edges
		{
			if(!firstEdgeMet)
			{
				firstEdgeMet = true;

				if(freq>-1)
				{
					for(tr1::unordered_map<double, set<int>* >::iterator iter = this->nodesByLabel.begin();iter!=nodesByLabel.end();iter++)
					{
						if(iter->second->size()<freq)
						{
							tr1::unordered_map<double, set<int>* >::iterator iter2 = this->nodesByLabel.find(iter->first);
							if(iter2==this->nodesByLabel.end())
							{
								cout<<"Error: isufnm44"<<endl;
								exit(0);
							}
							set<int>* nodesToRemove = iter2->second;
							for(set<int>::iterator iter1 = nodesToRemove->begin();iter1!=nodesToRemove->end();)
							{
								int nodeID = *iter1;
								iter1++;
								this->removeNode_IgnoreEdges(nodeID);
							}
						}
					}
				}
			}

			int id1;
			int id2;
			double label;
			data>>id1;
			data>>id2;
			data>>label;
			this->addEdge(id1, id2, label, edgeToFreq);
			if(numEdgesLoaded%1000000==0 && numEdgesLoaded>0)
			{
				char name[100];
				gethostname(name, 99);
			}
			numEdgesLoaded++;
		}
	}

	return true;
}

NodeX* GraphX::getNodeWithID(int nodeID)
{
	tr1::unordered_map<int, NodeX*>::iterator iter = nodes.find(nodeID);
	if(iter==nodes.end())
		return NULL;
	else
		return iter->second;
}

set<int>* GraphX::getNodesByLabel(double label)
{
	tr1::unordered_map<double, set<int>* >::iterator iter = nodesByLabel.find(label);
	if(iter==nodesByLabel.end())
		return NULL;
	return iter->second;
}

tr1::unordered_map<int,NodeX*>::const_iterator GraphX::getNodesIterator()
{
	return nodes.begin();
}

tr1::unordered_map<int,NodeX*>::const_iterator GraphX::getNodesEndIterator()
{
	return nodes.end();
}

/**
 * get edge label from the scr node to the destination node
 * if either the src node is not found, or the detination node is not found from the src node, then return 0.0001
 */
double GraphX::getEdgeLabel(int srcNodeID, int destNodeID)
{
	NodeX* srcNode;
	tr1::unordered_map<int, NodeX* >::iterator temp = nodes.find(srcNodeID);
	if(temp==nodes.end())
			return 0.0001;
	srcNode = (*temp).second;

	EdgeX* edge = (EdgeX*)srcNode->getEdgeForDestNode(destNodeID);
	if(edge==NULL)
		return 0.0001;

	return edge->getLabel();
}

int GraphX::getNumOfNodes()
{
	return this->nodes.size();
}

/**
 * return treu if the graph is connected
 */
bool GraphX::isConnected()
{
	if(nodes.size()<2)
		return true;

	//start from any node
	NodeX* node = nodes.begin()->second;
	map<int,NodeX*> visited;
	map<int,NodeX*> toVisit;

	toVisit.insert(std::pair<int, NodeX*>(node->getID(), node));
	while(toVisit.size()>0)
	{
		//1- pop a node for the to be visited list, 2- remove it from toVisit, and 3- add it to the visited list
		node = toVisit.begin()->second;//1
		toVisit.erase(toVisit.begin());//2
		visited.insert(std::pair<int, NodeX*>(node->getID(), node));//3
		//add its neighbors
		for(tr1::unordered_map<int, void*>::iterator iter = node->getEdgesIterator();iter!=node->getEdgesEndIterator();++iter)
		{
			int id = iter->first;
			if(visited.find(id)!=visited.end())
				continue;
			EdgeX* edge = (EdgeX*)iter->second;
			toVisit.insert(std::pair<int, NodeX*>(id, edge->getOtherNode()));
		}
	}

	if(visited.size()==nodes.size())
		return true;
	else
		return false;

}

/**
 * return true if the two graphs are the same
 */
bool GraphX::isTheSameWith(GraphX* otherG)
{
	if(this->getNumOfNodes()!=otherG->getNumOfNodes() ||
			this->getNumOfEdges()!=otherG->getNumOfEdges())
	{
		return false;
	}

	//compare using CL
	string thisCL = this->getCanonicalLabel();
	string otherCL = otherG->getCanonicalLabel();
	if(thisCL.compare(otherCL)!=0)
		return false;
	else
	{
		bool a = false;
		if(thisCL.at(0)=='X') a = true;

		if(!a)
			return true;
	}

	vector<map<int, int>* > result;
	tr1::unordered_map<int, tr1::unordered_set<int>*> domains_values;

	this->isIsomorphic(otherG, result, domains_values, -1, -1, false);

	bool b;

	if(result.size()>0)
	{
		b = true;
	}
	else
	{
		b = false;
	}

	for(vector<map<int, int>* >::iterator iter1 = result.begin();iter1!=result.end();iter1++)
	{
		delete (*iter1);
	}
	result.clear();

	return b;
}

/**
 * check whether two graphs are isomorphic or not
 * if isomorphic, then 'results' will have mapping between nodes from this graph to nodes in the given graph
 * 'this' is the query graph, while graph i the data graph, this means I can query this (smaller) in the data (Bigger)
 *  [*] this function does not consider edge label! FIX IT!!!!
 */
void GraphX::isIsomorphic(GraphX* graph, vector<map<int, int>* >& results, tr1::unordered_map<int, tr1::unordered_set<int>*>& domains_values, int restrictedDomainID, int restrictedNodeID, bool pruneByDomainValues, tr1::unordered_map<int, tr1::unordered_set<int>* >* postponedNodes, unsigned long maxIters)
{
	//a variable for counting
	numIterations = 0;

	//populate the order in which query graph will be traversed
	tr1::unordered_set<int> checked;
	int nodesOrder[getNumOfNodes()];//the result of this part

	//get the nodes order
	vector<int> toCheck;

	//check for the restricted node, if it exists let it be the first to check
	if(restrictedDomainID==-1)
		toCheck.push_back(nodes.begin()->second->getID());
	else
		toCheck.push_back(restrictedDomainID);

	int count = 0;
	while(toCheck.size()>0)
	{
		int current = toCheck.front();
		toCheck.erase(toCheck.begin());
		if(checked.find(current)!=checked.end())
			continue;
		nodesOrder[count] = current;
		checked.insert(current);
		count++;
		NodeX* currNode = nodes[current];

		int start = toCheck.size();

		for(tr1::unordered_map<int, void*>::iterator iter = currNode->getEdgesIterator();iter!=currNode->getEdgesEndIterator();++iter)
		{
			int otherID = iter->first;
			if(checked.find(otherID)==checked.end())
			{
				//put them in order based on the degree
				vector<int>::iterator iter1 = toCheck.begin();
				for(int i=0;i<start;i++)
					iter1++;
				for(;iter1!=toCheck.end();++iter1)
				{
					if(nodes[otherID]->getEdgesSize()<nodes[(*iter1)]->getEdgesSize())
						break;
				}
				toCheck.insert(iter1, otherID);
			}
		}
	}

	vector<Set_Iterator* > selected;
	tr1::unordered_set<int> selectedDataMap;//only for fast check for existence
	tr1::unordered_map<int, int> selectedQueryMap;//map value is query node id, value is its order
	//add it to the selected list
	Set_Iterator* isi = new Set_Iterator();

	//if the restricted node ID exists, then limit its domain to the restricted node
	if(restrictedDomainID==-1)
	{
		isi->se = graph->getNodesByLabel(this->getNodeWithID(nodesOrder[0])->getLabel());
		if(isi->se==NULL)
		{
			delete isi;
			return;
		}
	}
	else
	{
		isi->se = new set<int>();
		isi->se->insert(restrictedNodeID);
	}
	isi->it = isi->se->begin();
	selected.push_back(isi);
	selectedDataMap.insert(*(isi->it));
	selectedQueryMap.insert(std::pair<int, int>(nodesOrder[selected.size()-1], selected.size()-1));

	while(true)
	{
		numIterations++;
		if(postponedNodes!=0)
		{
			if(numIterations>maxIters)
			{
				tr1::unordered_map<int, tr1::unordered_set<int>* >::iterator iter = postponedNodes->find(restrictedDomainID);
				tr1::unordered_set<int>* nodesList;

				if(iter==postponedNodes->end())
				{
					nodesList = new tr1::unordered_set<int>();
					postponedNodes->insert(std::pair<int, tr1::unordered_set<int>* >(restrictedDomainID, nodesList));
				}
				else
				{
					nodesList = iter->second;
				}
				nodesList->insert(restrictedNodeID);

				if(Settings::debugMSG)
				{
					cout<<"isIsomorphic function exceeded the allowed processing limit! NodeID: "<<restrictedNodeID<<", DomainID: "<<restrictedDomainID<<endl;
					cout<<"maxIters = "<<maxIters<<", numIterations = "<<numIterations<<endl;
				}

				{
					vector<Set_Iterator* >::iterator iter = selected.begin();
					if(iter!=selected.end())
					{
						if(restrictedDomainID!=-1)
							delete (*iter)->se;
						delete (*iter);
						iter++;
						for(;iter!=selected.end();iter++)
						{
							(*iter)->se->clear();
							delete (*iter)->se;
							delete (*iter);
						}
					}
				}

				return;
			}
		}

		//take care of the counting part
		Set_Iterator* currentISI = selected.back();//.at(selected.size()-1);

		while(currentISI->isIterEnd())
		{
			numIterations++;

			//if we finished the domain of the first node
			if(selected.size()==1)
			{
				if(restrictedDomainID!=-1)
				{
					delete (*(selected.begin()))->se;
				}
				delete (*(selected.begin()));
				selected.clear();
				return;
			}

			//clear the data node associated with the last selected
			//selectedDataMap.erase(*(selected.back()->it));//commented on 27/10/2015
			//remove the last selected
			Set_Iterator* tempToDel = selected.back();
			selected.pop_back();
			delete tempToDel->se;
			delete tempToDel;
			//delete the query node associated with the last selected
			selectedQueryMap.erase(nodesOrder[selected.size()]);

			currentISI = selected.back();//.at(selected.size()-1);

			selectedDataMap.erase(*(currentISI->it));
			currentISI->it++;
			if(!currentISI->isIterEnd())//only the if line is added on 27/10/2015
				selectedDataMap.insert(*(currentISI->it));
		}
		//check current status
		NodeX* queryGraphNode = getNodeWithID(nodesOrder[selected.size()-1]);
		numIterations+=queryGraphNode->getEdgesSize();

		bool b = checkCurrentStatus(this, graph, selected, nodesOrder, selectedQueryMap);
		//if valid
		if(b)
		{
			//I found a solution
			if(selected.size()==this->getNumOfNodes())
			{
				map<int, int>* m = new map<int,int>();
				int c = 0;
				for(vector<Set_Iterator* >::iterator i = selected.begin();i!=selected.end();++i,c++)
				{
					m->insert(std::pair<int, int>(nodesOrder[c], *((*i)->it)));
				}
				results.push_back(m);

				if(restrictedDomainID!=-1)
				{
					//delete elements in selected
					for(vector<Set_Iterator* >::iterator iter = selected.begin();iter!=selected.end();iter++)
					{
						(*iter)->se->clear();
						delete (*iter)->se;
						delete (*iter);
					}
					return;
				}

				selectedDataMap.erase(*(selected.back()->it));
				selected.back()->it++;
				if(!selected.back()->isIterEnd())
					selectedDataMap.insert(*(selected.back()->it));
			}
			else//no solution found yet!
			{
				Set_Iterator* si = new Set_Iterator();
				si->se = new set<int>();
				//get extension from current selected (last si)
				//get the index of a query node connected to the query node to be selected for the check
				NodeX* temp = NULL;
				NodeX* tempN1 = this->getNodeWithID(nodesOrder[selected.size()]);
				vector<Set_Iterator* >::iterator tempIT = selected.begin();

				for(int i=0;i<selected.size();i++, tempIT++)
				{
					if(tempN1->isItConnectedWithNodeID(nodesOrder[i]))
					{
						temp = graph->getNodeWithID((*(*tempIT)->it));
						break;
					}
				}
				for(tr1::unordered_map<int, void*>::iterator iter = temp->getEdgesIterator();iter!=temp->getEdgesEndIterator();++iter)
				{
					numIterations++;

					if(numIterations>maxIters)
						break;

					int otherNodeID = iter->first;

					NodeX* otherNode = graph->getNodeWithID(otherNodeID);

					//check for the AC_3ed domain, to check for node occurrence, if not then no need to check it
					if(pruneByDomainValues && domains_values.size()>0)
					{
						tr1::unordered_set<int>* tempD = domains_values[nodesOrder[selected.size()]];

						if(tempD->find(otherNodeID)==tempD->end())
						{
							continue;
						}
					}

					//check for node label
					if(otherNode->getLabel()!=tempN1->getLabel())
					{
						continue;
					}

					if(selectedDataMap.find(otherNodeID)==selectedDataMap.end())
					{
						bool b = true;
						//check for other edges
						for(tr1::unordered_map<int, void*>::iterator iter1 = tempN1->getEdgesIterator();iter1!=tempN1->getEdgesEndIterator();++iter1)
						{
							numIterations++;

							if(numIterations>maxIters)
								break;

							EdgeX* edge = (EdgeX*)iter1->second;
							tr1::unordered_map<int, int>::iterator itr = selectedQueryMap.find(edge->getOtherNode()->getID());
							if(itr == selectedQueryMap.end())
							{
								continue;
							}
							int tempOrder = itr->second;
							double edgeLabel = edge->getLabel();
							if(!otherNode->isItConnectedWithNodeID((*(selected[tempOrder]->it)), edgeLabel))
							{
								b = false;
								break;
							}
						}

						if(b)
						{
							si->se->insert(otherNodeID);
						}
					}
					else
					{
					}
				}

				//I discovered a bug when testing the efficient graph, this fix should work now, otherwise the old code is below
				if(si->se->size()>0)
				{
					si->it = si->se->begin();
					selected.push_back(si);
					selectedDataMap.insert(*(si->it));
					selectedQueryMap.insert(std::pair<int, int>(nodesOrder[selected.size()-1], selected.size()-1));
				}
				else
				{
					si->it = si->se->end();
					delete si->se;
					delete si;
					//
					selectedDataMap.erase(*(selected.back()->it));
					selected.back()->it++;
					if(!selected.back()->isIterEnd())
						selectedDataMap.insert(*(selected.back()->it));
				}
			}
		}
		else
		{
			selectedDataMap.erase(*(selected.back()->it));
			selected.back()->it++;
			if(!selected.back()->isIterEnd())
				selectedDataMap.insert(*(selected.back()->it));
		}
	}

	//delete elements in selected
	vector<Set_Iterator* >::iterator iter = selected.begin();
	if(iter!=selected.end())
	{
		if(restrictedDomainID!=-1)
			delete (*iter)->se;
		delete (*iter);
		iter++;
		for(;iter!=selected.end();iter++)
		{
			(*iter)->se->clear();
			delete (*iter)->se;
			delete (*iter);
		}
	}
	selected.clear();
}

/**
 * get the canonical label of this graph
 */
string GraphX::getCanonicalLabel()
{
	if(CL.length()>0)
	{
		return CL;
	}

	CL = CanonicalLabel::generate(this);

	return CL;
}

ostream& operator<<(ostream& os, const GraphX& g)
{
	os<<"# t 1\n";
	//output the nodes
	for(tr1::unordered_map<int,NodeX*>::const_iterator ii=g.nodes.begin(); ii!=g.nodes.end(); ++ii)
	{
		NodeX* node = ii->second;
		os<<"v "<<node->getID()<<" "<<node->getLabel()<<"\n";
	}

	//output the edges
	tr1::unordered_set<string> savedEdges;//list to keep track of already saved edges
	for(tr1::unordered_map<int,NodeX*>::const_iterator ii=g.nodes.begin(); ii!=g.nodes.end(); ++ii)
	{
		NodeX* node = ii->second;
		for(tr1::unordered_map<int, void*>::iterator iter1 = node->getEdgesIterator();iter1!=node->getEdgesEndIterator();iter1++)
		{
			EdgeX* edge = (EdgeX*)iter1->second;

			//check whether it has been added before or not
			string sig = intToString(node->getID())+"_"+doubleToString(edge->getLabel())+"_"+intToString(edge->getOtherNode()->getID());
			if(node->getID()<edge->getOtherNode()->getID())
				sig = intToString(edge->getOtherNode()->getID())+"_"+doubleToString(edge->getLabel())+"_"+intToString(node->getID());
			if(savedEdges.find(sig)==savedEdges.end())
			{
				savedEdges.insert(sig);
				os<<"e "<<node->getID()<<" "<<edge->getOtherNode()->getID()<<" "<<edge->getLabel()<<"\n";
			}
		}
	}

    return os;
}

//Destructor
GraphX::~GraphX()
{
	for(tr1::unordered_map<int,NodeX*>::iterator ii=nodes.begin(); ii!=nodes.end(); ++ii)
	{
		delete ii->second;
	}
	nodes.clear();

	for( tr1::unordered_map<double, set<int>* >::iterator ii=nodesByLabel.begin(); ii!=nodesByLabel.end(); ++ii)
	{
		delete ii->second;
	}
	nodesByLabel.clear();
}
