/**
 * Copyright 2016 Ehab Abdelhamid, Ibrahim Abdelaziz

Pattern.cpp

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
 * represents a frequent subgraph
 */
#include<iostream>
#include<limits>
#include <stdlib.h>
#include "Pattern.h"
#include "utils.h"
#include "Settings.h"

int Pattern::maxPatternID = 0;

bool PrimaryGraph::isTheSameWith(PrimaryGraph* otherPG)
{
	if(this->graph->getNumOfNodes()!=otherPG->graph->getNumOfNodes() ||
				this->graph->getNumOfEdges()!=otherPG->graph->getNumOfEdges())
		return false;

	vector<map<int, int>* > result;
	tr1::unordered_map<int, tr1::unordered_set<int>*> domains_values;

	this->getGraph()->isIsomorphic(otherPG->getGraph(), result, domains_values);

	bool b;

	if(result.size()>0)
		b = true;
	else
		b = false;

	for(vector<map<int, int>* >::iterator iter1 = result.begin();iter1!=result.end();iter1++)
	{
		delete (*iter1);
	}
	result.clear();

	return b;
}

bool PrimaryGraph::isTheSameWith(GraphX* otherPG)
{
	if(this->graph->getNumOfNodes()!=otherPG->getNumOfNodes() ||
					this->graph->getNumOfEdges()!=otherPG->getNumOfEdges())
		return false;

	vector<map<int, int>* > result;
	tr1::unordered_map<int, tr1::unordered_set<int>*> domains_values;

	this->getGraph()->isIsomorphic(otherPG, result, domains_values);

	bool b;

	if(result.size()>0)
		b = true;
	else
		b = false;

	for(vector<map<int, int>* >::iterator iter1 = result.begin();iter1!=result.end();iter1++)
	{
		delete (*iter1);
	}
	result.clear();

	return b;
}

Pattern::Pattern(GraphX* graph, bool copyGraph)
{
	this->graphCopied = copyGraph;
	this->resultExact = false;
	this->ID = Pattern::maxPatternID;
	Pattern::maxPatternID++;
	this->frequency = -1;
	if(copyGraph)
		this->graph = new GraphX(graph);
	else
		this->graph = graph;

	this->maxIters = Settings::postponeNodesAfterIterations;

	init();
}

Pattern::Pattern(Pattern* pattern)
{
	this->ID = Pattern::maxPatternID;
	Pattern::maxPatternID++;
	this->frequency = pattern->getFrequency();
	//get a similar graph
	this->graph = new GraphX(pattern->getGraph());
	this->graphCopied = true;
	this->predictedTime = pattern->getPredictedTime();
	this->setInvalidCol(pattern->getInvalidCol(), pattern->getPredictedValids());

	this->resultExact = pattern->isResultExact();
	if(resultExact)
		this->frequency = pattern->getFrequency();

	init();

	this->Combine(pattern);
}

void Pattern::resetMNI()
{
	if(subtaskingFixed>0)
	{
		for(int i=0;i<this->graph->getNumOfNodes();i++)
		{
			tempMNITable[i] = 0;
			postponedNodes_mniTable[i] = 0;
		}
	}
}

int Pattern::getID()
{
	return ID;
}

void Pattern::init()
{
	st_counter++;
	int numNodes = graph->getNumOfNodes();
	for(int i=0;i<numNodes;i++)
		occurences.push_back(new tr1::unordered_set<int>());
	predictedTime = 0;

	tempMNITable = 0;
	remainingSubtasks = 0;

	if(Settings::postponeExpensiveNodes)
		postponeExpensiveNodes = true;
	else
		postponeExpensiveNodes = false;

	tempMNITable = new int[this->graph->getNumOfNodes()];
	postponedNodes_mniTable = new int[this->graph->getNumOfNodes()];

	predictedPattern = 0;
}

void Pattern::setInvalidCol(int invC, int pValids)
{
	invalidCol = invC;
	predictedValids = pValids;
	if(invC==-1)
		postponeExpensiveNodes = false;
}

void Pattern::generatePrimaryGraphs()
{
	if(this->primaryGraphs.size()>0)
	{
		if(Settings::debugMSG)
		{
			cout<<"Already generated, return!"<<endl;
		}
		return;
	}

	long long start1 = getmsofday();

	//the set of already removed edges signatures, node1ID+"_"+node2ID if node1ID>node2ID, other wise it is node2ID+"_"+node1ID
	set<string> alreadyRemoved;
	//go over all nodes, try to remove an edge one at at time
	for(tr1::unordered_map<int,NodeX*>::const_iterator iter = graph->getNodesIterator();iter!=graph->getNodesEndIterator();++iter)
	{
		NodeX* node1 = iter->second;
		//if this node has one edge, then the other node (in this edge) should be treated as the source node
		if(node1->getEdgesSize()==1)
			continue;
		for(tr1::unordered_map<int, void*>::iterator iter = node1->getEdgesIterator();iter!=node1->getEdgesEndIterator();++iter)
		{
			EdgeX* edge = ((EdgeX*)iter->second);
			NodeX* node2 = edge->getOtherNode();

			//check whether this edge has been removed before
			string node1ID = intToString(node1->getID());
			string node2ID = intToString(node2->getID());
			string sig;
			if(node1->getID()>node2->getID())
				sig = node1ID+"_"+node2ID;
			else
				sig = node2ID+"_"+node1ID;
			if(alreadyRemoved.find(sig)==alreadyRemoved.end())
			{
				alreadyRemoved.insert(sig);
				GraphX* rGraph = new GraphX(graph);
				rGraph->removeEdge(node1->getID(), node2->getID());
				if(!rGraph->isConnected())
				{
					delete rGraph;
					continue;
				}

				//check if the given primary graph is already generated before
				bool exists = false;
				for(list<PrimaryGraph* >::iterator iter = primaryGraphs.begin();iter!=primaryGraphs.end();iter++)
				{
					PrimaryGraph* pgTemp = (*iter);
					if(pgTemp->isTheSameWith(rGraph))
					{
						exists = true;
						break;
					}
				}
				if(!exists)
				{
					PrimaryGraph* pg = new PrimaryGraph();
					pg->setValues(rGraph, node1->getID(), edge);
					this->primaryGraphs.push_back(pg);
				}
				else
				{
					delete rGraph;
				}
			}
		}
	}
}

/**
 * get the joining primarygraphs of the two patterns
 */
set<std::pair<PrimaryGraph*, PrimaryGraph* > > Pattern::getJoiningPG(Pattern* pattern)
{
	set<std::pair<PrimaryGraph*, PrimaryGraph* > > llist;

	for(list<PrimaryGraph*>::iterator ii1=primaryGraphs.begin(); ii1!=primaryGraphs.end(); ++ii1)
	{
		PrimaryGraph* pg1 = (*ii1);
		for( list<PrimaryGraph* >::iterator ii2=pattern->primaryGraphs.begin(); ii2!=pattern->primaryGraphs.end(); ++ii2)
		{
			PrimaryGraph* pg2 = (*ii2);
			if(pg1->isTheSameWith(pg2))
			{
				llist.insert(std::pair<PrimaryGraph*, PrimaryGraph* >(pg1, pg2));
			}
		}
	}

	return llist;
}

void Pattern::setFrequency(int newFreq)
{
	frequency = newFreq;
}

int Pattern::getFrequency()
{
	if(frequency>-1)
		return frequency;

	if(graph->getNumOfNodes()==0)
		return -1;

	int min = occurences[0]->size();
	for(int i=1;i<graph->getNumOfNodes();i++)
	{
		if(min>occurences[i]->size())
			min = occurences[i]->size();
	}

	return min;
}

void Pattern::borrowTimeInfor(Pattern* otherPattern, int nWorkers)
{
	//in case we already set subtaskingfixed, then there is no need to borrow information, because it was set earlier
	if(this->subtaskingFixed>1)
		return;

	if(Settings::debugMSG)
	{
		cout<<"Pattern#"<<this->getID()<<" is borrowing time information from Pattern#"<<otherPattern->getID()<<"this.Subtasking = "<<subtaskingFixed<<", other.subtasking = "<<otherPattern->getSubtaskingValueFixed()<<endl;
		cout<<"Borrowed predicted time is: "<<otherPattern->getPredictedTime()<<endl;
	}

	this->predictedTime = otherPattern->getPredictedTime();
	this->setSubtasking(otherPattern->getSubtaskingValueFixed(), nWorkers);
	this->setInvalidCol(otherPattern->getInvalidCol(), otherPattern->getPredictedValids());

	this->resultExact = otherPattern->isResultExact();
	if(resultExact)
		this->frequency = otherPattern->getFrequency();

	this->setMaxIters(otherPattern->getMaxIters());
}

void Pattern::addNode(int nodeID, int patternNodeID)
{
	occurences[patternNodeID]->insert(nodeID);
}

string Pattern::toString()
{
	string str = intToString(getFrequency())+"\n";
	for(int i=0;i<graph->getNumOfNodes();i++)
	{
		str = str + intToString(i) + " with label: "+doubleToString(graph->getNodeWithID(i)->getLabel())+"\nNodes list:\n";

		for(tr1::unordered_set<int>::iterator iter = occurences[i]->begin();iter!=occurences[i]->end();++iter)
		{
			str = str +","+intToString(*iter);
		}
		str = str+"\n";
	}
	return str;
}

void Pattern::Combine(Pattern* otherP, int addToID)
{
	invalidateFrequency();
	for(int i=0;i<graph->getNumOfNodes();i++)
	{
		for(tr1::unordered_set<int>::iterator iter = otherP->getOccurences()->at(i)->begin();iter!=otherP->getOccurences()->at(i)->end();++iter)
		{
			occurences[i]->insert((*iter)+addToID);
		}
	}
}

void Pattern::extend(int srcID, int destID, double destLabel, double edgeLabel)
{

	if(graph->getNodeWithID(destID)==NULL)
		graph->AddNode(destID, destLabel);
	graph->addEdge(srcID, destID, edgeLabel);

	int numNodes = graph->getNumOfNodes();
	while(occurences.size()<numNodes)
		occurences.push_back(new tr1::unordered_set<int>());

	//delete primary graphs
	for( list<PrimaryGraph* >::iterator ii=primaryGraphs.begin(); ii!=primaryGraphs.end(); ++ii)
	{
		delete (*ii);
	}

	primaryGraphs.clear();
}

int Pattern::getSize()
{
	return graph->getNumOfEdges();
}

bool Pattern::hasUniqueLabels()
{
	tr1::unordered_set<double> labels;
	for(tr1::unordered_map<int, NodeX*>::const_iterator iter = graph->getNodesIterator();iter!=graph->getNodesEndIterator();iter++)
	{
		double label = iter->second->getLabel();
		if(labels.find(label)==labels.end())
			labels.insert(label);
		else
			return false;
	}
	return true;
}

void Pattern::setSubtasking(int st, int nWorkers)
{
	//if we do not do task division
	if(!Settings::divideBigTasks)
	{
		st = 1;
	}

	if(Settings::fixedNumSubtasks>-1)
	{
		if(subtaskingFixed>1)
			return;
		else
			st = Settings::fixedNumSubtasks;
	}

	if(st>nWorkers)
		st = nWorkers;

	if(Settings::debugMSG)
		cout<<"setSubTasking. st = "<<st<<endl;
	this->subtasking = st;
	this->subtaskingFixed = st;
	remainingSubtasks = st;
	if(st>0)
	{
		if(tempMNITable==0)
			tempMNITable = new int[this->graph->getNumOfNodes()];
		if(postponedNodes_mniTable==0)
			postponedNodes_mniTable = new int[this->graph->getNumOfNodes()];
		for(int i=0;i<this->graph->getNumOfNodes();i++)
		{
			tempMNITable[i] = 0;
			postponedNodes_mniTable[i] = 0;
		}
	}
}

/**
 * returns -1 only if no more subtasks are remaining
 * return -2 if we need to re-run but not using the postponed nodes option
 * return -3 if we need to re-run using
 * otherwise return the support
 */
int Pattern::subTaskDone(int* mniTable, int* ppn_mniTable, int support)
{
	this->remainingSubtasks--;

	if(Settings::debugMSG)
	{
		cout<<"Master side: for candidate#"<<this->ID<<", remainingSubtasks="<<remainingSubtasks<<endl;
		cout<<"MNI table BEFORE adding new values:"<<endl;
	}

	for(int i=0;i<this->graph->getNumOfNodes();i++)
	{
		if(Settings::debugMSG)
			cout<<tempMNITable[i]<<", ";
		this->tempMNITable[i] += mniTable[i];
	}

	if(Settings::debugMSG)
		cout<<endl;

	if(ppn_mniTable!=0)
	{
		if(Settings::debugMSG)
		{
			cout<<"Postponed MNI table BEFORE adding new values:"<<endl;
		}

		for(int i=0;i<this->graph->getNumOfNodes();i++)
		{
			if(Settings::debugMSG)
				cout<<this->postponedNodes_mniTable[i]<<", ";
			this->postponedNodes_mniTable[i] += ppn_mniTable[i];
		}

		if(Settings::debugMSG)
			cout<<endl;
	}

	if(Settings::debugMSG)
	{
		cout<<"MNI table AFTER adding new values:"<<endl;
		for(int i=0;i<this->graph->getNumOfNodes();i++)
		{
			cout<<tempMNITable[i]<<", ";
		}
		cout<<endl;
	}

	if(ppn_mniTable!=0)
	{
		if(Settings::debugMSG)
		{
			cout<<"Postponed MNI table AFTER adding new values:"<<endl;
			for(int i=0;i<this->graph->getNumOfNodes();i++)
			{
				cout<<postponedNodes_mniTable[i]<<", ";
			}
			cout<<endl;
		}
	}

	if(this->remainingSubtasks==0)
	{
		if(this->getInvalidCol()==-1)
		{
			int min1 = std::numeric_limits<int>::max();
			for(int i=0;i<this->graph->getNumOfNodes();i++)
			{
				int temp = tempMNITable[i];
				if(temp<min1)
					min1 = temp;
			}

			int min2 = std::numeric_limits<int>::max();
			for(int i=0;i<this->graph->getNumOfNodes();i++)
			{
				int temp = tempMNITable[i]+postponedNodes_mniTable[i];
				if(temp<min2)
					min2 = temp;
			}

			if(min1==-1)
			{
				cout<<"Error 4563: -1 should not be returned";
				exit(0);
			}
			if(min2<min1)
			{
				cout<<"Error 4163: min2 should always be min1 or more!";
				exit(0);
			}

			//decide whether we need to rerun or not
			if(min1!=min2)
			{
				if(min1>=support)
				{
					return -2;
				}

				if(min2>=support)
				{
					return -2;
				}
			}

			return min1;
		}
		else
		{
			if(Settings::debugMSG)
				cout<<"Invalid column is used ... mni col size = "<<tempMNITable[this->getInvalidCol()]<<"+postponedNodes_mniTable: "<<postponedNodes_mniTable[this->getInvalidCol()]<<endl;

			if(tempMNITable[this->getInvalidCol()]>=support)// && postponedNodes_mniTable[this->getInvalidCol()]!=0)
			{
				if(Settings::debugMSG)
					cout<<"return frequency to repeat processing [1]"<<endl;//cout<<"return -2 [1]"<<endl;
				return tempMNITable[this->getInvalidCol()];
			}
			if((tempMNITable[this->getInvalidCol()]+postponedNodes_mniTable[this->getInvalidCol()])>=support && postponedNodes_mniTable[this->getInvalidCol()]!=0)
			{
				if(Settings::debugMSG)
					cout<<"return -2 [2]: "<<(tempMNITable[this->getInvalidCol()]+postponedNodes_mniTable[this->getInvalidCol()])<<">= support = "<<support<<endl;
				return -2;
			}

			return tempMNITable[this->getInvalidCol()];
		}
	}
	return -1;
}

void Pattern::makeIDNegative()
{
	if(ID>0)
		ID = ID*-1;
}

Pattern::~Pattern()
{
	st_counter--;
	if(this->graphCopied)
		delete graph;
	for(int i=0;i<occurences.size();i++)
		delete occurences.at(i);

	for( list<PrimaryGraph* >::iterator ii=primaryGraphs.begin(); ii!=primaryGraphs.end(); ++ii)
	{
		delete (*ii);
	}

	primaryGraphs.clear();
	if(tempMNITable!=0)
		delete[] tempMNITable;

	if(postponedNodes_mniTable!=0)
		delete[] postponedNodes_mniTable;
}
