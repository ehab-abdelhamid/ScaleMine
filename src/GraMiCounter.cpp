/**
 * Copyright 2016 Ehab Abdelhamid, Ibrahim Abdelaziz

GraMiCounter.cpp

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
 * responsible for computing the support of a given subgraph in an input graph given a support threshold
 * the exact version implements the algorithms in: "GRAMI: Frequent Subgraph and Pattern Mining in a Single Large Graph", VLDB 2014
 */
#include <tr1/unordered_map>
#include <tr1/unordered_set>
#include <limits>
#include <math.h>
#include <cstdlib>
#include<sys/time.h>
#include "GraMiCounter.h"
#include "utils.h"
#include "Settings.h"

int GraMiCounter::numByPassedNodes;
int GraMiCounter::numSideEffectNodes;
bool GraMiCounter::useAC3;
int Settings::samplesBucket_accur = 40;

//utility functions
long long getmsofday_3()
{
   struct timeval tv;
   struct timezone tz;
   gettimeofday(&tv, &tz);
   return (long long)tv.tv_sec*1000 + tv.tv_usec/1000;
}

string getSig(int a, int b, double el)
{
	char ch[20];

	string sig;
	if(a<b)
		sig = intToString(a)+"_"+intToString(b)+"_"+doubleToString(el);
	else
		sig = intToString(b)+"_"+intToString(a)+"_"+doubleToString(el);
	return sig;
}

//function to delete results space
void deleteResults(vector<map<int, int>* >& result)
{
	for(vector<map<int, int>* >::iterator iter1 = result.begin();iter1!=result.end();iter1++)
	{
		delete (*iter1);
	}
	result.clear();
}

/**
 * GraMi based approximate counter, for the 'pattern' in the 'graph' given the 'support' minimum threshold.
 * the 'approximate' parameter indicates whther we allow approximation or not, if approximate==-1 then it is not approximate
 */
int GraMiCounter::isFrequent(GraphX* graph, Pattern* pattern, int support, double approximate)
{
	int freq = 0;

	//create domains
	tr1::unordered_map<int, tr1::unordered_set<int>*> domains_values;
	GraphX* pg = pattern->getGraph();
	for(tr1::unordered_map<int, NodeX*>::const_iterator iter = pg->getNodesIterator();iter!=pg->getNodesEndIterator();++iter)
	{
		int varNodeID = iter->first;
		domains_values.insert(std::pair<int, tr1::unordered_set<int>*>(varNodeID, new tr1::unordered_set<int>()));
	}

	// [*] apply saved substructure optimization

	//insert domains values (considering node, neighborhood count, and degree consistencies)
	for(tr1::unordered_map<int, tr1::unordered_set<int>*>::iterator iter = domains_values.begin();iter!=domains_values.end();iter++)
	{
		//get the pattern node
		NodeX* pNode = pattern->getGraph()->getNodeWithID(iter->first);
		tr1::unordered_set<int>* currentDomain = domains_values[pNode->getID()];

		//check node consistency
		set<int>* nodes_SameLabel = graph->getNodesByLabel(pNode->getLabel());

		for(set<int>::iterator iter1 = nodes_SameLabel->begin();iter1!=nodes_SameLabel->end();iter1++)
		{
			NodeX* dNode = graph->getNodeWithID((*iter1));
			if(dNode->getEdgesSize()>=pNode->getEdgesSize() &&	//node degree consistency
					true)	// [*] apply neighborhood consistency
			{
				currentDomain->insert(dNode->getID());
			}
		}
		if(currentDomain->size()<support)
			return 0;
	}

	int patternSize = pattern->getGraph()->getNumOfNodes();
	long long start = getmsofday_3();

	//apply arc consistency constraint
	AC_3(graph, domains_values, pattern, support);

	long long end = getmsofday_3();
	long long elapsed = end - start;
	//check for domains lengths
	for(int i=0;i<patternSize;i++)
	{
		if(domains_values.find(i)->second->size()<support)
		{
			for(tr1::unordered_map<int, tr1::unordered_set<int>*>::iterator iter = domains_values.begin();iter!=domains_values.end();iter++)
				delete iter->second;

			return 0;
		}
	}

	// [*] set domains order
	int* ordering  = new int[patternSize];

	for(int i=0;i<patternSize;i++)
	{
		ordering[i] = i;
	}

	//create solutions data structure
	map<int, set<int>*> domains_solutions;
	for(int i=0;i<patternSize;i++)
	{
		domains_solutions.insert(std::pair<int, set<int>* >(i, new set<int>()));
	}

	for(int i=0;i<patternSize;i++)
	{
		int domainID = ordering[i];

		// [*] apply automorphism

		//go over elements in the current domain, check if a solution exists
		tr1::unordered_set<int>* currentDomain = domains_values.find(domainID)->second;

		if(Settings::debugMSG)
			cout<<"old currentDomain size:"<<currentDomain->size()<<endl;

		if(approximate!=-1)
		{
			map<int ,int> idMap;
			int c = 0;
			for(tr1::unordered_set<int>::iterator iter = currentDomain->begin();iter!=currentDomain->end();iter++)
			{
				idMap.insert(std::pair<int, int>(c, *iter));
				c++;
			}
			delete currentDomain;
			domains_values.erase(domainID);
			currentDomain = new tr1::unordered_set<int>();
			domains_values.insert(std::pair<int, tr1::unordered_set<int>* >(domainID, currentDomain));

			while(currentDomain->size()<(idMap.size()*approximate))
			{
				int r = rand()%idMap.size();
				currentDomain->insert(idMap.at(r));
			}
		}

		if(Settings::debugMSG)
			cout<<"new currentDomain size:"<<currentDomain->size()<<endl;

		int counter = 0;
		for(tr1::unordered_set<int>::iterator iter = currentDomain->begin();iter!=currentDomain->end();++iter)
		{
			//cout<<"DID:"<<domainID<<", c:"<<counter<<endl;
			counter++;
			int nodeID = (*iter);
			bool b = false;
			vector<map<int, int>* > result;
			//check if this node has been passed previously as a correct solution or not
			if(domains_solutions.find(domainID)->second->find(nodeID)!=domains_solutions.find(domainID)->second->end())
			{
				b = true;
			}
			else
			{
				pattern->getGraph()->isIsomorphic(graph, result, domains_values, domainID, nodeID);
				if(result.size()>0)
					b = true;
				else
				{
					//in case the remaining + existing solutions can not satisfy the support, break
					if(approximate==-1)
					{
						if(currentDomain->size()-counter+domains_solutions.find(domainID)->second->size()<support)
						{
							deleteResults(result);
							break;
						}
					}
					else
					{
						if(currentDomain->size()-counter+domains_solutions.find(domainID)->second->size()<(support*approximate))
						{
							deleteResults(result);
							break;
						}
					}
				}
			}
			if(b)
			{
				if(approximate==-1)
				{
					//there is a solution for this node, add all valid node values to the solutions domain
					for(vector<map<int, int>* >::iterator iter1 = result.begin();iter1!=result.end();iter1++)
					{
						map<int, int>* currentMapping = (*iter1);
						for(map<int, int>::iterator iter2 = currentMapping->begin();iter2!=currentMapping->end();++iter2)
						{
							int dID = iter2->first;
							int nID = iter2->second;
							domains_solutions.find(dID)->second->insert(nID);
						}
					}
				}
				else
				{
					domains_solutions.find(domainID)->second->insert(nodeID);
				}
			}
			if(approximate==-1)
			{
				if(domains_solutions.find(domainID)->second->size()>=support)
				{
					deleteResults(result);
					break;
				}
			}
			else
			{
				if(domains_solutions.find(domainID)->second->size()/approximate>=support)
				{
					deleteResults(result);
					break;
				}
			}

			deleteResults(result);
		}

		//in case the solution can not satisfy the support, break
		if(approximate==-1)
		{
			if(domains_solutions.find(domainID)->second->size()<support)
				break;
		}
		else
		{
			int tf = domains_solutions.find(domainID)->second->size()/approximate;
			if(tf<support)
			{
				break;
			}
		}
	}

	//delete the domains
	for(tr1::unordered_map<int, tr1::unordered_set<int>*>::iterator iter = domains_values.begin();iter!=domains_values.end();iter++)
	{
		delete iter->second;
	}

	//get MNI frequency using the domains solutions
	if(approximate==-1)
	{
		int min = std::numeric_limits<int>::max();
		for(map<int, set<int>*>::iterator iter = domains_solutions.begin();iter!=domains_solutions.end();iter++)
		{
			int temp = iter->second->size();
			if(temp<min)
				min = temp;
		}
		freq = min;
	}
	else //in case of approximation
	{
		int min = std::numeric_limits<int>::max();
		for(map<int, set<int>*>::iterator iter = domains_solutions.begin();iter!=domains_solutions.end();iter++)
		{
			int domainID = iter->first;
			int temp = iter->second->size();
			//temp = domains_values.find(domainID)->second->size();
			temp = temp / approximate;
			if(temp<min)
				min = temp;
		}
		freq = min;
	}

	//delete domains solutions
	for(map<int, set<int>*>::iterator iter = domains_solutions.begin();iter!=domains_solutions.end();iter++)
	{
		delete iter->second;
	}

	//delete the ordering
	delete ordering;

	return freq;
}

int GraMiCounter::isFrequent_adv(GraphX* graph, Pattern* pattern, int support, int subTaskNum, int subTaskMax, int* mniTable, tr1::unordered_map<int, tr1::unordered_set<int>* >* postponedNodes, unsigned long maxIters)
{
	int freq = 0;

	for(int i=0;i<pattern->getGraph()->getNumOfNodes();i++)
	{
		mniTable[i] = 0;
	}

	//create domains
	tr1::unordered_map<int, tr1::unordered_set<int>*> domains_values;
	GraphX* pg = pattern->getGraph();
	for(tr1::unordered_map<int, NodeX*>::const_iterator iter = pg->getNodesIterator();iter!=pg->getNodesEndIterator();++iter)
	{
		int varNodeID = iter->first;

		tr1::unordered_set<int>* currentDomain = new tr1::unordered_set<int>();
		domains_values.insert(std::pair<int, tr1::unordered_set<int>*>(varNodeID, currentDomain));

		//insert domains values (considering node, neighborhood count, and degree consistencies)
		NodeX* pNode = pattern->getGraph()->getNodeWithID(varNodeID);
		int pNodeEdgesSize = pNode->getEdgesSize();
		//check node consistency
		set<int>* nodes_SameLabel = graph->getNodesByLabel(pNode->getLabel());

		currentDomain->insert(nodes_SameLabel->begin(), nodes_SameLabel->end());


		//case where we distribute the task on different workers and there is an invalid column
		//for this case, we limit the domain (inavlidCol domain) to nodes with IDs divisible by taskNumber
		if(varNodeID==pattern->getInvalidCol() && subTaskNum>-1)
		{
			for(tr1::unordered_set<int>::iterator iter1 = currentDomain->begin();iter1!=currentDomain->end();)
			{
				if(graph->getNodeWithID((*iter1))->getEdgesSize()>=pNodeEdgesSize &&	//node degree consistency
						true)
				{
					//for subtasks when the invalidCol is set and nodeID belongs to current subtask, remove the node
					if(((*iter1)%subTaskMax)!=subTaskNum)
					{
						iter1 = currentDomain->erase((iter1));
					}
					else
						iter1++;
				}
				else
				{
					iter1 = currentDomain->erase((iter1));

				}
			}
		}
		else
		{
			for(tr1::unordered_set<int>::iterator iter1 = currentDomain->begin();iter1!=currentDomain->end();)
			{
				if(graph->getNodeWithID((*iter1))->getEdgesSize()>=pNodeEdgesSize &&	//node degree consistency
						true)
				{
					iter1++;
				}
				else
				{

					iter1 = currentDomain->erase((iter1));
					if(currentDomain->size()<support)
					{
						for(tr1::unordered_map<int, tr1::unordered_set<int>*>::iterator iter1 = domains_values.begin();iter1!=domains_values.end();iter1++)
						{
							iter1->second->clear();
							delete iter1->second;
						}

						if(Settings::debugMSG)
							cout<<"Cannot create domain: "<<varNodeID<<" with enough nodes"<<endl;

						return 0;
					}
				}
			}
		}
	}

	bool useAC3 = GraMiCounter::useAC3;
	int numEdges = pattern->getGraph()->getNumOfEdges();

	long totalNumCandidNodes = 0;
	for(tr1::unordered_map<int, tr1::unordered_set<int>*>::iterator iter = domains_values.begin();iter!=domains_values.end();iter++)
	{
		totalNumCandidNodes+=iter->second->size();
	}

	if(subTaskNum>-1)
	{
		useAC3 = false;
	}
	else
	{
		useAC3 = true;
	}

	int patternSize = pattern->getGraph()->getNumOfNodes();

	long long start = getmsofday_3();

	if(useAC3)
		AC_3(graph, domains_values, pattern, support);

	long long end = getmsofday_3();
	long long elapsed = end - start;

	if(Settings::debugMSG)
		if(useAC3)
			cout<<"AC_3 took "<<(elapsed/1000)<<" sec and "<<(elapsed%1000)<<" ms"<<endl;
		else
			cout<<"AC_3 is not used."<<endl;

	if(subTaskNum>-1 && pattern->getInvalidCol()>-1)
	{
		if(Settings::debugMSG)
		{
			for(tr1::unordered_map<int, tr1::unordered_set<int>*>::iterator iter = domains_values.begin();iter!=domains_values.end();iter++)
			{
				cout<<"domain size"<<iter->second->size()<<" [invCol is set]"<<endl;
			}
		}
	}
	else
	{
		//check for domains lengths
		for(tr1::unordered_map<int, tr1::unordered_set<int>*>::iterator iter = domains_values.begin();iter!=domains_values.end();iter++)
		{
			if(Settings::debugMSG)
				cout<<"domain size"<<iter->second->size()<<endl;

			if(iter->second->size()<support)
			{
				for(iter = domains_values.begin();iter!=domains_values.end();iter++)
				{
					iter->second->clear();
					delete iter->second;
				}

				if(Settings::debugMSG)
					cout<<"Proved infrequent by AC3[1] for candidate graph with id = "<<pattern->getID()<<endl;

				return 0;
			}
		}
	}

	//check for unique labels optimization, if applicable then return true (as we already checked in the previous  step for infrequentness)
	//this check is only done for not divided tasks. You can optimize for this
	if(subTaskNum==-1 && pattern->hasUniqueLabels() && GraMiCounter::isItAcyclic(*(pattern->getGraph())))
	{
		if(Settings::debugMSG)
			cout<<"BINGO1!"<<endl;

		int freq = 0;
		int min = std::numeric_limits<int>::max();
		for(tr1::unordered_map<int, tr1::unordered_set<int>*>::iterator iter1 = domains_values.begin();iter1!=domains_values.end();iter1++)
		{
			int temp = iter1->second->size();
			iter1->second->clear();
			delete iter1->second;//new
			if(temp<min)
				min = temp;
			mniTable[iter1->first] = temp;
		}
		freq = min;

		return freq;
	}

	// [*] set domains order
	int* ordering  = new int[patternSize];
	setDomainsOrder(domains_values, ordering, pattern);

	//create solutions data structure
	map<int, set<int>*> domains_solutions;
	map<int, set<int>*> side_solutions;
	for(int i=0;i<patternSize;i++)
	{
		domains_solutions.insert(std::pair<int, set<int>* >(i, new set<int>()));
		side_solutions.insert(std::pair<int, set<int>* >(i, new set<int>()));
	}

	unsigned long oldMaxIters_ = maxIters;

	for(int i=0;i<patternSize;i++)
	{
		int postponedNodesCounter = 0;
		maxIters = oldMaxIters_;

		if(i>0 && pattern->getInvalidCol()>-1)
		{
			mniTable[pattern->getInvalidCol()] = domains_solutions.find(pattern->getInvalidCol())->second->size();

			for(tr1::unordered_map<int, tr1::unordered_set<int>*>::iterator iter1 = domains_values.begin();iter1!=domains_values.end();iter1++)
			{
				iter1->second->clear();
				delete iter1->second;//new
			}

			//delete domains solutions
			int domainsSolutionsSize = domains_solutions.find(pattern->getInvalidCol())->second->size();
			for(map<int, set<int>*>::iterator iter = domains_solutions.begin();iter!=domains_solutions.end();iter++)
			{
				iter->second->clear();
				delete iter->second;
			}

			//delete side solutions
			for(map<int, set<int>*>::iterator iter = side_solutions.begin();iter!=side_solutions.end();iter++)
			{
				iter->second->clear();
				delete iter->second;
			}

			delete[] ordering;

			return domainsSolutionsSize;
		}
		int domainID = ordering[i];

		//go over elements in the current domain, check if a solution exists
		tr1::unordered_set<int>* currentDomain = domains_values.find(domainID)->second;

		//get number of nodes that can be checked for this subtask
		int numNodesToCheck = 0;

		if(subTaskNum>-1)
		{
			for(tr1::unordered_set<int>::iterator iter = currentDomain->begin();iter!=currentDomain->end();++iter)
			{
				int nodeID = (*iter);
				if((nodeID%subTaskMax)!=subTaskNum)
					continue;
				numNodesToCheck++;
			}
		}
		else
		{
			numNodesToCheck = currentDomain->size();
		}

		double nodesShare = 1.0/subTaskMax;

		if(nodesShare>1)
		{
			cout<<"NodesShare = "<<nodesShare<<endl<<flush;
			exit(0);
		}

		std::tr1::unordered_set<int> placesSet;
		int counter = 0;

		for(tr1::unordered_set<int>::iterator iter = currentDomain->begin();iter!=currentDomain->end();++iter)
		{
			int nodeID;
			nodeID = (*iter);

			if(subTaskNum>-1 && (nodeID%subTaskMax)!=subTaskNum)
			{
				continue;
			}

			counter++;

			bool b = false;
			vector<map<int, int>* > result;

			//check if this node has been passed previously as a correct solution or not
			if(side_solutions.find(domainID)->second->find(nodeID)!=side_solutions.find(domainID)->second->end())
			{
				b = true;
				GraMiCounter::numSideEffectNodes++;
			}
			else
			{
				int MLDecision = -1;
				int proposedMLDecision = -1;

				{
					unsigned long startTime = getmsofday();

					tr1::unordered_map<int, tr1::unordered_set<int>* >::iterator iterr;
					int temp = 0;

					if(postponedNodes!=0)
					{
						iterr = postponedNodes->find(domainID);

						if(iterr!=postponedNodes->end())
						{
							temp = iterr->second->size();
						}
					}

					pattern->getGraph()->isIsomorphic(graph, result, domains_values, domainID, nodeID, true, postponedNodes, maxIters);

					unsigned long elapsed = getmsofday() - startTime;
					if(elapsed>100000)
					{
						cout<<"checking: "<<counter<<"/"<<currentDomain->size()<<". NodeID = "<<nodeID<<", with #edges = "<<graph->getNodeWithID(nodeID)->getEdgesSize()<<endl;
						cout<<"Took "<<elapsed<<"ms"<<endl;
						cout<<"Max Iter = "<<maxIters<<", num of iterations = "<<pattern->getGraph()->numIterations<<endl;
					}

					if(postponedNodes!=0 && iterr!=postponedNodes->end())
					{
						postponedNodesCounter+=(iterr->second->size()-temp);

						if(postponedNodesCounter%100==0)
						{
							postponedNodesCounter = 1;

							unsigned long oldMaxIters = maxIters;
							maxIters = maxIters * 2;
							if(maxIters<oldMaxIters)
							{
								maxIters = oldMaxIters;
								if(Settings::debugMSG)
								{
									cout<<"Trying to increasing maxIters failed!"<<endl;
								}
							}
							else
							{
								if(Settings::debugMSG)
								{
									cout<<"Increasing maxIters from "<<oldMaxIters<<" to "<<maxIters<<endl;
								}
							}

						}
					}
				}

				if(result.size()>0)
				{
					b = true;
				}
				else
				{
					//in case the remaining + existing solutions can not satisfy the support, break
					if(Settings::smartBreak && subTaskNum==-1)
					{
						long maximumNumPossibleValids = currentDomain->size()-counter+domains_solutions.find(domainID)->second->size();
						if(postponedNodes!=0)
						{
							tr1::unordered_map<int, tr1::unordered_set<int>* >::iterator iter = postponedNodes->find(domainID);
							if(iter!=postponedNodes->end())
								maximumNumPossibleValids+=iter->second->size();
						}
						if(maximumNumPossibleValids<support)
						{
							deleteResults(result);
							if(Settings::debugMSG)
								cout<<"Counter:"<<counter<<", Ramaining ["<<(currentDomain->size()-counter)<<"] solutions ["<<domains_solutions.find(domainID)->second->size()<<"] Can not satisfy the required support"<<endl;
							break;
						}
					}
					else
					{
						if(pattern->getInvalidCol()!=-1)
						{

							int remaining = numNodesToCheck-counter;

							double progress = counter/((double)numNodesToCheck);
							long numValidsAndPostponed = domains_solutions.find(domainID)->second->size();
							if(postponedNodes!=0)
							{
								tr1::unordered_map<int, tr1::unordered_set<int>* >::iterator iter = postponedNodes->find(domainID);
								if(iter!=postponedNodes->end())
									numValidsAndPostponed+=iter->second->size();
							}
							double lhs = numValidsAndPostponed-(pattern->getPredictedValids()*progress*nodesShare)+remaining;
							double rhs = ((support-pattern->getPredictedValids())*Settings::validsPerTaskMargin)*nodesShare;

							if(Settings::smartBreak && lhs<rhs && remaining>(numNodesToCheck*0.05))
							{
								deleteResults(result);

								if(Settings::debugMSG)
									cout<<"predicted num of valids allow us to break. Remaining/Total = "<<remaining<<"/"<<numNodesToCheck<<" rhs = "<<rhs<<", lhs = "<<lhs<<". Valid = "<<domains_solutions.find(domainID)->second->size()<<". Predicted valids = "<<pattern->getPredictedValids()<<"Progress = "<<progress<<", nodesShare = "<<nodesShare<<endl;

								GraMiCounter::numByPassedNodes+=remaining;

								for(;iter!=currentDomain->end();++iter)
								{
									int nodeID;
									nodeID = (*iter);

									if((nodeID%subTaskMax)!=subTaskNum)
									{
										continue;
									}

									domains_solutions.find(domainID)->second->insert(nodeID);
								}

								break;
							}
						}
					}
				}
			}

			if(b)
			{
				//there is a solution for this node, add all valid node values to the solutions domain
				for(vector<map<int, int>* >::iterator iter1 = result.begin();iter1!=result.end();iter1++)
				{
					map<int, int>* currentMapping = (*iter1);
					//this is done only for not divided tasks
					if(subTaskNum==-1)
					{
						for(map<int, int>::iterator iter2 = currentMapping->begin();iter2!=currentMapping->end();++iter2)
						{
							int dID = iter2->first;
							int nID = iter2->second;
							domains_solutions.find(dID)->second->insert(nID);
							side_solutions.find(dID)->second->insert(nID);
						}
					}
					else
					{
						for(map<int, int>::iterator iter2 = currentMapping->begin();iter2!=currentMapping->end();++iter2)
						{
							int dID = iter2->first;
							int nID = iter2->second;

							if(subTaskNum>-1 && (nID%subTaskMax)!=subTaskNum)
							{
								continue;
							}

							domains_solutions.find(dID)->second->insert(nID);
							side_solutions.find(dID)->second->insert(nID);
						}
					}
				}
			}

			if(!Settings::fullCount)
			{
				if(domains_solutions.find(domainID)->second->size()>=support)
				{
					deleteResults(result);

					if(Settings::debugMSG)
						cout<<"The domain solution size exceeds the required support: "<<domains_solutions.find(domainID)->second->size()<<endl;
					break;
				}
			}

			deleteResults(result);
		}

		//in case the solution can not satisfy the support, break

		long numValidsAndPostponed = domains_solutions.find(domainID)->second->size();
		if(postponedNodes!=0)
		{
			tr1::unordered_map<int, tr1::unordered_set<int>* >::iterator iter = postponedNodes->find(domainID);
			if(iter!=postponedNodes->end())
				numValidsAndPostponed+=iter->second->size();
		}

		if(subTaskNum==-1 && numValidsAndPostponed<support)
		{
			if(Settings::debugMSG)
				cout<<"Current domain: "<<domainID<<" has solutions with size < support"<<endl<<flush;
			break;
		}
	}

	//delete the domains
	for(tr1::unordered_map<int, tr1::unordered_set<int>*>::iterator iter = domains_values.begin();iter!=domains_values.end();iter++)
	{
		iter->second->clear();
		delete iter->second;
	}

	//get MNI frequency using the domains solutions
	int min = std::numeric_limits<int>::max();
	for(map<int, set<int>*>::iterator iter = domains_solutions.begin();iter!=domains_solutions.end();iter++)
	{
		int temp = iter->second->size();
		if(temp<min)
			min = temp;
		mniTable[iter->first] = temp;
	}
	freq = min;

	//delete domains solutions
	for(map<int, set<int>*>::iterator iter = domains_solutions.begin();iter!=domains_solutions.end();iter++)
	{
		iter->second->clear();
		delete iter->second;
	}

	//delete side solutions
	for(map<int, set<int>*>::iterator iter = side_solutions.begin();iter!=side_solutions.end();iter++)
	{
		iter->second->clear();
		delete iter->second;
	}

	//delete the ordering
	delete[] ordering;

	return freq;
}

/**
 * This version is utilized for building an approximate search space
 * invalidCol: column that is expected to be invalid
 * predictedValids: number of predicted valid nodes in the invalid column
 * exact returns true if the resulted frequency is a result of optimization steps, which is exact
 */
int GraMiCounter::isFrequent_approx(GraphX* graph, Pattern* pattern, int support, int& invalidCol, int& predictedValids, bool allowEarlyBreak, bool& exact, unsigned long& numVisiteNodes, unsigned long& numIterations, vector<unsigned long>& listOfNumOfIters, unsigned long& predictedTime)
{
	exact = false;
	numVisiteNodes = 0;
	numIterations = 0;

	//no invalid column yet! this is a column discovered being invalid while processing that column
	invalidCol = -1;
	predictedValids = -1;

	if(Settings::debugMSG)
		cout<<"Target support = "<<support<<endl;
	int freq = 0;
	predictedTime = 0;

	long long start = getmsofday_3();
	//create domains
	tr1::unordered_map<int, tr1::unordered_set<int>*> domains_values;

	GraphX* pg = pattern->getGraph();

	//decide whether to use AC3 or not
	bool useAC3 = GraMiCounter::useAC3;

	if(pattern->hasUniqueLabels() && GraMiCounter::isItAcyclic(*(pattern->getGraph())))
	{
		useAC3 = true;
	}
	else
	{
		long totalNumCandidNodes = 0;
		for(tr1::unordered_map<int, NodeX*>::const_iterator iter = pg->getNodesIterator();iter!=pg->getNodesEndIterator();++iter)
		{
			int varNodeID = iter->first;
			double domainLabel = pg->getNodeWithID(varNodeID)->getLabel();
			int numCandidsForDomain = graph->getNodesByLabel(domainLabel)->size();
			totalNumCandidNodes+=numCandidsForDomain;
		}
		if(pattern->getGraph()->getNumOfEdges()==0)
		{
			cout<<"skdsnfnff"<<endl<<flush;
			exit(0);
		}

		int numEdges = pattern->getGraph()->getNumOfEdges();

		if(graph->getNumOfEdges()/graph->getNumOfNodes()>10)
		{
			useAC3 = false;
		}
		else
		{
			useAC3 = true;
		}
	}

	//end of decide use AC3
	if(Settings::debugMSG)
		cout<<"useAC3 = "<<useAC3<<endl;


	for(tr1::unordered_map<int, NodeX*>::const_iterator iter = pg->getNodesIterator();iter!=pg->getNodesEndIterator();++iter)
	{
		int varNodeID = iter->first;
		tr1::unordered_set<int>* currentDomain = new tr1::unordered_set<int>();
		domains_values.insert(std::pair<int, tr1::unordered_set<int>*>(varNodeID, currentDomain));

		if(useAC3)//GraMiCounter::useAC3)
		{
			//insert domains values (considering node, neighborhood count, and degree consistencies)
			NodeX* pNode = pattern->getGraph()->getNodeWithID(varNodeID);
			int pNodeEdgesSize = pNode->getEdgesSize();
			//check node consistency
			set<int>* nodes_SameLabel = graph->getNodesByLabel(pNode->getLabel());

			for(set<int>::iterator iter1 = nodes_SameLabel->begin();iter1!=nodes_SameLabel->end();iter1++)
			{
				if(graph->getNodeWithID((*iter1))->getEdgesSize()>=pNodeEdgesSize &&	//node degree consistency
						//dNode->isNeighborhoodConsistent(pNode))
						true)	// [*] apply neighborhood consistency
				{
					currentDomain->insert((*iter1));
				}
				else
				{
					//currentDomain->insert((*iter1));
				}
			}
			if(currentDomain->size()<support)
			{
				for(tr1::unordered_map<int, tr1::unordered_set<int>*>::iterator iter1 = domains_values.begin();iter1!=domains_values.end();iter1++)
					delete iter1->second;
				exact = true;
				return 0;
			}
		}
	}

	long long end = getmsofday_3();
	long long elapsed = end - start;

	predictedTime+=elapsed;
	int patternSize = pattern->getGraph()->getNumOfNodes();

	//apply arc consistency constraint
	if(useAC3)//GraMiCounter::useAC3)
	{
		start = getmsofday_3();

		AC_3(graph, domains_values, pattern, support);

		end = getmsofday_3();
		elapsed = end - start;

		if(Settings::debugMSG)
			cout<<"AC_3 took "<<(elapsed/1000)<<" sec and "<<(elapsed%1000)<<" ms"<<endl;

		//check for domains lengths
		for(tr1::unordered_map<int, tr1::unordered_set<int>*>::iterator iter = domains_values.begin();iter!=domains_values.end();iter++)
		{
			if(Settings::debugMSG)
				cout<<"domain size: "<<iter->second->size()<<endl;
			if(iter->second->size()<support)
			{
				for(iter = domains_values.begin();iter!=domains_values.end();iter++)
					delete iter->second;

				if(Settings::debugMSG)
					cout<<"domain size is less than support"<<endl;
				exact = true;
				return 0;
			}
		}
	}

	//check for unique labels optimization, if applicable then return true (as we already checked in the previous  step for infrequentness)
	if(useAC3 && pattern->hasUniqueLabels() && GraMiCounter::isItAcyclic(*(pattern->getGraph())))//GraMiCounter::useAC3
	{
		if(Settings::debugMSG)
			cout<<"BINGO2!"<<endl;
		exact = true;

		int min = std::numeric_limits<int>::max();
		for(tr1::unordered_map<int, tr1::unordered_set<int>*>::iterator iter = domains_values.begin();iter!=domains_values.end();iter++)
		{
			int temp = iter->second->size();
			delete iter->second;//new
			if(temp<min)
				min = temp;
		}

		return min;
	}

	// [*] set domains order
	int* ordering  = new int[patternSize];
	setDomainsOrder(domains_values, ordering, pattern);

	double* approximator;
	int* fullSize;
	if(Settings::useSearchSpacePrediction)
	{
		approximator = new double[domains_values.size()];
		fullSize = new int[domains_values.size()];
	}

	//create solutions data structure
	map<int, set<int>*> domains_solutions;
	map<int, set<int>*> side_solutions;
	for(int i=0;i<patternSize;i++)
	{
		domains_solutions.insert(std::pair<int, set<int>* >(i, new set<int>()));
		side_solutions.insert(std::pair<int, set<int>* >(i, new set<int>()));
		if(Settings::useSearchSpacePrediction)
		{
			approximator[i] = -1;
			fullSize[i] = -1;
		}
	}

	start = getmsofday_3();

	for(int i=0;i<patternSize;i++)
	{
		long long startPredict = getmsofday_3();

		int domainID = ordering[i];
		//cout<<"Checking domain ID: "<<domainID<<endl;

		//go over elements in the current domain, check if a solution exists
		tr1::unordered_set<int>* currentDomain = domains_values.find(domainID)->second;

		//insert domain values (considering node, neighborhood count, and degree consistencies)
		long long start = getmsofday_3();
		NodeX* pNode = pattern->getGraph()->getNodeWithID(domainID);

		if(!useAC3)//GraMiCounter::useAC3)
		{
			int pNodeEdgesSize = pNode->getEdgesSize();
			//check node consistency
			set<int>* nodes_SameLabel = graph->getNodesByLabel(pNode->getLabel());

			for(set<int>::iterator iter1 = nodes_SameLabel->begin();iter1!=nodes_SameLabel->end();iter1++)
			{
				if(graph->getNodeWithID((*iter1))->getEdgesSize()>=pNodeEdgesSize &&	//node degree consistency
						//dNode->isNeighborhoodConsistent(pNode))
						true)	// [*] apply neighborhood consistency
				{
					currentDomain->insert((*iter1));
				}
				else
				{
					//currentDomain->insert((*iter1));
				}
			}
			if(currentDomain->size()<support)
			{
				for(tr1::unordered_map<int, tr1::unordered_set<int>*>::iterator iter1 = domains_values.begin();iter1!=domains_values.end();iter1++)
					delete iter1->second;

				predictedTime+=(getmsofday_3()-startPredict);

				delete[] ordering;
				delete[] approximator;
				delete[] fullSize;

				return 0;
			}
		}

		long long end = getmsofday_3();
		long long elapsed = end - start;

		if(Settings::debugMSG)
			cout<<"Creating domain["<<domainID<<"] took "<<(elapsed/1000)<<" sec and "<<(elapsed%1000)<<" ms"<<endl;
		//end of filling domain with nodes

		//std::tr1::unordered_set<int> placesSet;
		int counter = 0;

		int* idList;
		if(Settings::useSearchSpacePrediction)
		{
			idList = new int[Settings::maxNumberOfSamples];
			//std::tr1::unordered_set<int> placesSet;
			int counter = 0;
			tr1::unordered_set<int>::iterator iter = currentDomain->begin();

			if(currentDomain->size()==0)
			{
				cout<<"err783743"<<endl<<flush;
				exit(0);
			}

			while(counter<Settings::maxNumberOfSamples)
			{
				int r = rand()%currentDomain->size();
				{
					for(int j=0;j<r;j++)
						iter++;
					idList[counter] = (*iter);
					iter = currentDomain->begin();
					counter++;
				}
			}

		}


		double targetAvg = ((double)support)/currentDomain->size();
		targetAvg = targetAvg*Settings::samplesBucket;//added for the new sample scheme
		int numSampleAvgs = Settings::maxNumberOfSamples/Settings::samplesBucket;
		double* sampleAvg = new double[numSampleAvgs];
		for(int i=0;i<numSampleAvgs;i++)
			sampleAvg[i] = -1;

		counter=0;

		int validNodes = 0;
		int invalidNodes = 0;

		int totalValidNodes = 0;
		int totalInvalidNodes = 0;

		predictedTime+=(getmsofday_3()-startPredict);

		startPredict = getmsofday_3();

		double lastSampleMean;
		double lastSampleSD;

		for(tr1::unordered_set<int>::iterator iter = currentDomain->begin();iter!=currentDomain->end();++iter)
		{
			int nodeID;

			if(Settings::useSearchSpacePrediction)
			{
				if(counter>0 && (counter%Settings::samplesBucket)==0)
				{
					srand (time(NULL));

					int numOfSampleAvgs = counter/Settings::samplesBucket;
					if(numOfSampleAvgs==0)
					{
						cout<<"kjkdswqqq3333"<<endl<<flush;
						exit(0);
					}
					sampleAvg[numOfSampleAvgs-1] = ((double)validNodes)/Settings::samplesBucket;//old calculation

					if(Settings::debugMSG)
						cout<<counter<<" setting sample avg for "<<(numOfSampleAvgs-1)<<" to "<<sampleAvg[numOfSampleAvgs-1]<<". valid nodes = "<<validNodes<<", invalidNodes = "<<invalidNodes<<endl;
					validNodes = 0;
					invalidNodes = 0;

					//check if the sampling we used is good enough to take a decision
					if(allowEarlyBreak && counter>=Settings::minNumberOfSamples)
					{
						double sampleSD = 0;
						double sampleMean = 0;
						//vew calculations based on the binomial distribution
						{
						//estimated distribution mean
						double p = ((double)totalValidNodes)/(totalValidNodes+totalInvalidNodes);

						//calculating sample mean and starndard deviation
						sampleMean = p*Settings::samplesBucket;
						sampleSD = sqrt(Settings::samplesBucket*p*(1-p));
						sampleSD = sampleSD / sqrt(Settings::samplesBucket);
						}

						lastSampleMean = sampleMean;
						lastSampleSD = sampleSD;

						//compute confidence interval
						//for 99% z-value is 2.58
						double z = 2.58;
						double lowInterval = sampleMean - (z*sampleSD);
						double highInterval = sampleMean + (z*sampleSD);

						if(Settings::debugMSG)
							cout<<"sampleMean = "<<sampleMean<<", sampleSD = "<<sampleSD<<", lowInterval = "<<lowInterval<<", highInterval = "<<highInterval<<endl;
						//now check if the required threshold ratio is beyond the interval, if yes then break (we are confident enough of whether it is frequent or not frequent)
						if(targetAvg<lowInterval || targetAvg>highInterval)
						{
							if(Settings::debugMSG)
								cout<<"Early break ... good: sampleMean = "<<sampleMean<<", targetAvg = "<<targetAvg<<endl;
							break;
						}
					}
				}

				if(counter>=Settings::maxNumberOfSamples)
				{
					break;
				}

				nodeID = idList[counter];
			}
			else
			{
				nodeID = (*iter);
			}

			counter++;

			bool b = false;
			vector<map<int, int>* > result;
			//check if this node has been passed previously as a correct solution or not
			if(side_solutions.find(domainID)->second->find(nodeID)!=side_solutions.find(domainID)->second->end())
			{
				b = true;
			}
			else
			{
				{
					tr1::unordered_map<int, tr1::unordered_set<int>* >* ppn = 0;
					if(Settings::postponeExpensiveNodes)
						ppn = new tr1::unordered_map<int, tr1::unordered_set<int>* >();
					if(useAC3)//GraMiCounter::useAC3)
						pattern->getGraph()->isIsomorphic(graph, result, domains_values, domainID, nodeID);
					else
						pattern->getGraph()->isIsomorphic(graph, result, domains_values, domainID, nodeID, false);

					if(Settings::postponeExpensiveNodes)
						delete ppn;

					numVisiteNodes++;
					numIterations += pattern->getGraph()->numIterations;
					//put them in order
					vector<unsigned long>::iterator i=listOfNumOfIters.begin();
					unsigned long value = pattern->getGraph()->numIterations;
					for(;i!=listOfNumOfIters.end();i++)
						if(value<(*i))
							break;
					listOfNumOfIters.insert(i, value);
				}

				if(result.size()>0)
				{
					b = true;
				}
				else
				{
					//in case the remaining + existing solutions can not satisfy the support, break
					if(!Settings::useSearchSpacePrediction && currentDomain->size()-counter+domains_solutions.find(domainID)->second->size()<support)
					{
						deleteResults(result);
						break;
					}
				}
			}

			if(b)
			{
				//there is a solution for this node, add all valid node values to the solutions domain
				for(vector<map<int, int>* >::iterator iter1 = result.begin();iter1!=result.end();iter1++)
				{
					map<int, int>* currentMapping = (*iter1);
					for(map<int, int>::iterator iter2 = currentMapping->begin();iter2!=currentMapping->end();++iter2)
					{
						int dID = iter2->first;
						int nID = iter2->second;
						domains_solutions.find(dID)->second->insert(nID);
						side_solutions.find(dID)->second->insert(nID);
					}
				}
				validNodes++;
				totalValidNodes++;
			}
			else
			{
				invalidNodes++;
				totalInvalidNodes++;
			}

			if(!Settings::fullCount && !Settings::useSearchSpacePrediction)
			{
				if(domains_solutions.find(domainID)->second->size()>=support)
				{
					deleteResults(result);

					if(Settings::debugMSG)
						cout<<"The domain solution size exceeds the required support: "<<domains_solutions.find(domainID)->second->size()<<endl;
					break;
				}
			}

			deleteResults(result);
		}

		if(counter==0)
		{
			cout<<"errretyru78"<<endl<<flush;
			exit(0);
		}
		predictedTime+=((getmsofday_3()-startPredict)*(((double)currentDomain->size())/counter));

		//compute approximator
		approximator[domainID] = 0;
		int numOfSampleAvgs = counter/Settings::samplesBucket;
		for(int ii=0;ii<numOfSampleAvgs;ii++)
		{
			approximator[domainID]+=sampleAvg[ii];
		}

		if(numOfSampleAvgs==0)
		{
			cout<<"popoyuytgv"<<endl<<flush;
			exit(0);
		}
		approximator[domainID] = approximator[domainID] / numOfSampleAvgs;

		if(Settings::debugMSG)
			cout<<"approximator["<<domainID<<"] = "<<approximator[domainID]<<endl;

		if(Settings::useSearchSpacePrediction)
			fullSize[domainID] = currentDomain->size();
		delete[] sampleAvg;

		if(Settings::debugMSG)
			cout<<"SampleAvg freed"<<endl;

		if(Settings::useSearchSpacePrediction)
		{
			delete[] idList;
		}

		if(Settings::debugMSG)
			cout<<"idList freed"<<endl;

		//in case the solution can not satisfy the support, break
		if(Settings::useSearchSpacePrediction)
		{
			//option 1: break right away
			double validNodesPredicted = approximator[domainID]*((double)fullSize[domainID]);
			if(validNodesPredicted<support)
			{
				if(Settings::debugMSG)
					cout<<"["<<patternSize<<"]"<<flush;

				//the used one:
				//if a column found to be invalid, use it and report the time, so that if this column really took so much time, then it is divided
				invalidCol = domainID;
				predictedValids = validNodesPredicted;

				if(Settings::debugMSG)
					cout<<"Breaking from domain "<<domainID<<" as [approximator[domainID]*((double)fullSize[i])<support]"<<validNodesPredicted<<endl;
				break;
			}
		}
		else
		{
			if(domains_solutions.find(domainID)->second->size()<support)
				break;
		}
	}

	end = getmsofday_3();
	elapsed = end - start;

	if(Settings::debugMSG)
		cout<<endl<<"Searching for valids took "<<(elapsed/1000)<<" sec and "<<(elapsed%1000)<<" ms"<<endl;

	//delete the domains
	for(tr1::unordered_map<int, tr1::unordered_set<int>*>::iterator iter = domains_values.begin();iter!=domains_values.end();iter++)
	{
		delete iter->second;
	}

	//get MNI frequency using the domains solutions
	if(!Settings::useSearchSpacePrediction)
	{
		int min = std::numeric_limits<int>::max();
		for(map<int, set<int>*>::iterator iter = domains_solutions.begin();iter!=domains_solutions.end();iter++)
		{
			int temp = iter->second->size();
			if(temp<min)
				min = temp;
		}
		freq = min;
	}
	else
	{
		int min = std::numeric_limits<int>::max();
		for(int i = 0;i<domains_solutions.size();i++)
		{
			int temp;
			if(approximator[i]>-1)
				temp = approximator[i]*((double)fullSize[i]);
			else
				temp = 0;
			if(temp<min)
				min = temp;

			if(Settings::debugMSG)
				cout<<"Domain ["<<i<<"]: approximator = "<<approximator[i]<<", fullSize = "<<fullSize[i]<<endl;
		}
		freq = min;
		if(Settings::debugMSG)
				cout<<"FrequencyCD = "<<freq<<endl;
	}

	delete[] approximator;
	delete[] fullSize;

	//delete domains solutions
	for(map<int, set<int>*>::iterator iter = domains_solutions.begin();iter!=domains_solutions.end();iter++)
	{
		delete iter->second;
	}

	//delete side solutions
	for(map<int, set<int>*>::iterator iter = side_solutions.begin();iter!=side_solutions.end();iter++)
	{
		delete iter->second;
	}

	//delete the ordering
	delete[] ordering;

//	cout<<"Q10"<<endl<<flush;

	return freq;
}

void GraMiCounter::AC_3(GraphX* graph, tr1::unordered_map<int, tr1::unordered_set<int>*>& domains_values, Pattern* pattern, int support)
{
	GraMiCounter::AC_3(graph, domains_values, pattern->getGraph(), support, pattern->getInvalidCol());
}

void GraMiCounter::AC_3(GraphX* graph, tr1::unordered_map<int, tr1::unordered_set<int>*>& domains_values, GraphX* pGraph, int support, int invalidCol)
{
	map<string,pair_with_edge*> arcs;

	//create pairs of connection from the pattern
	for(tr1::unordered_map<int, NodeX*>::const_iterator iter = pGraph->getNodesIterator();iter!=pGraph->getNodesEndIterator();++iter)
	{
		NodeX* pNode = iter->second;
		for(tr1::unordered_map<int, void*>::iterator iter1 = pNode->getEdgesIterator();iter1!=pNode->getEdgesEndIterator();iter1++)
		{
			EdgeX* edge = (EdgeX*)(iter1->second);
			int id1 = pNode->getID();
			int id2 = edge->getOtherNode()->getID();
			string sig = getSig(id1, id2, edge->getLabel());

			if(arcs.find(sig)!=arcs.end())
				continue;

			pair_with_edge* pwe = new pair_with_edge();
			pwe->id1 = id1;
			pwe->id2 = id2;
			pwe->edgeLabel = edge->getLabel();
			pwe->minDomainSize = domains_values.find(pwe->id1)->second->size();
			if(pwe->minDomainSize<domains_values.find(pwe->id2)->second->size())
				pwe->minDomainSize=domains_values.find(pwe->id2)->second->size();

			map<string,pair_with_edge*>::iterator iter = arcs.begin();
			for(;iter!=arcs.end();iter++)
			{
				pair_with_edge* temp_pwe = iter->second;
				//cout<<iter->first<<endl<<flush;
				if(temp_pwe->minDomainSize>pwe->minDomainSize)
					break;
			}
			if(arcs.find(sig)!=arcs.end())//new 2April
				delete pwe;
			else
				arcs.insert(iter, std::pair<string, pair_with_edge*>(sig, pwe));
		}
	}

	//order arcs to start with invalidcolumn if exists

	while(arcs.size()>0)
	{
		pair_with_edge* pwe = 0;
		string sig;

		if(invalidCol!=-1)
		{
			map<string,pair_with_edge*>::iterator iter = arcs.begin();
			for(;iter!=arcs.end();iter++)
			{
				pair_with_edge* temp_pwe = iter->second;
				if(temp_pwe->id1==invalidCol || temp_pwe->id2==invalidCol)
				{
					pwe = iter->second;
					sig = iter->first;
					break;
				}
			}
		}

		if(pwe==0)
		{
			pwe = arcs.begin()->second;
			sig = arcs.begin()->first;
		}

		//this lookup is repeated in the refine function, think of a way not to repeat it
		tr1::unordered_set<int>* D1 = domains_values.find(pwe->id1)->second;
		tr1::unordered_set<int>* D2 = domains_values.find(pwe->id2)->second;
		int old_size1 = D1->size();
		int old_size2 = D2->size();

		if(refine(graph, domains_values, D1, D2, pwe, support))
		{
			map<string,pair_with_edge*>::iterator iter = arcs.begin();
			for(;iter!=arcs.end();iter++)
			{
				pair_with_edge* temp_pwe = iter->second;
				delete temp_pwe;
			}
			arcs.clear();
			return;
		}

		//add affected arcs
		NodeX* pNode = pGraph->getNodeWithID(pwe->id1);
		if(old_size1!=D1->size() && pNode->getEdgesSize()>1)
		{
			int id1 = pwe->id1;
			for(tr1::unordered_map<int, void*>::iterator iter1 = pNode->getEdgesIterator();iter1!=pNode->getEdgesEndIterator();iter1++)
			{
				EdgeX* edge = (EdgeX*)(iter1->second);

				int id2 = edge->getOtherNode()->getID();
				if(id2==pwe->id2)
					continue;
				string sig = getSig(id1, id2, pwe->edgeLabel);

				if(arcs.find(sig)!=arcs.end())
					;//delete pwe;
				else
				{
					//create pwe
					pair_with_edge* pwe = new pair_with_edge();
					pwe->id1 = id1;
					pwe->id2 = id2;
					pwe->edgeLabel = edge->getLabel();
					pwe->minDomainSize = domains_values.find(pwe->id1)->second->size();
					if(pwe->minDomainSize<domains_values.find(pwe->id2)->second->size())
					pwe->minDomainSize=domains_values.find(pwe->id2)->second->size();

					arcs.insert(std::pair<string, pair_with_edge*>(sig, pwe));
				}
			}
		}

		pNode = pGraph->getNodeWithID(pwe->id2);
		if(old_size2!=D2->size() && pNode->getEdgesSize()>1)
		{
			int id1 = pwe->id2;
			for(tr1::unordered_map<int, void*>::iterator iter1 = pNode->getEdgesIterator();iter1!=pNode->getEdgesEndIterator();iter1++)
			{
				EdgeX* edge = (EdgeX*)(iter1->second);

				int id2 = edge->getOtherNode()->getID();
				if(id2==pwe->id1)
					continue;
				string sig = getSig(id1, id2, pwe->edgeLabel);

				if(arcs.find(sig)!=arcs.end())
					;//delete pwe;
				else
				{
					//create pwe
					pair_with_edge* pwe = new pair_with_edge();
					pwe->id1 = id2;
					pwe->id2 = id1;
					pwe->edgeLabel = edge->getLabel();
					pwe->minDomainSize = domains_values.find(pwe->id1)->second->size();
					if(pwe->minDomainSize<domains_values.find(pwe->id2)->second->size())
					pwe->minDomainSize=domains_values.find(pwe->id2)->second->size();

					arcs.insert(std::pair<string, pair_with_edge*>(sig, pwe));
				}
			}
		}

		arcs.erase(sig);
		delete pwe;
	}

	//clear arcs list
	map<string,pair_with_edge*>::iterator iter = arcs.begin();
	for(;iter!=arcs.end();iter++)
	{
		pair_with_edge* temp_pwe = iter->second;
		delete temp_pwe;
	}
	arcs.clear();
}

bool GraMiCounter::refine(GraphX* graph, tr1::unordered_map<int, tr1::unordered_set<int>*>& domains_values, tr1::unordered_set<int>* D1, tr1::unordered_set<int>* D2, pair_with_edge* pwe, int support)
{
	double edgeLabel = pwe->edgeLabel;
	//set to iterate over the smaller domain
	if(D1->size()>D2->size())
	{
		tr1::unordered_set<int>* temp = D1;
		D1 = D2;
		D2 = temp;
	}
	//go over each node in the domain1
	set<int> new_D2;
	for(tr1::unordered_set<int>::iterator iter = D1->begin();iter!=D1->end();)
	{
		NodeX* node = graph->getNodeWithID((*iter));
		//go over its edges
		bool deleteIt = true;
		for(tr1::unordered_map<int, void*>::iterator iter1 = node->getEdgesIterator();iter1!=node->getEdgesEndIterator();++iter1)
		{
			EdgeX* edge = (EdgeX*)(iter1->second);
			if(edge->getLabel()!=edgeLabel)
				continue;
			tr1::unordered_set<int>::iterator iter_f = D2->find(edge->getOtherNode()->getID());
			if(iter_f!=D2->end())
			{
				deleteIt = false;
				new_D2.insert(*iter_f);
			}
		}

		if(deleteIt)
		{
			iter = D1->erase(iter);
			if(D1->size()<support)
				return true;
		}
		else
			++iter;
	}

	//add qualified nodes to D2
	//option 1 (a bit slower)
	/*D2->clear();
	D2->insert(new_D2.begin(), new_D2.end());
	if(D2->size()<support)
		return true;*/
	//option 2
	for(tr1::unordered_set<int>::iterator iter = D2->begin();iter!=D2->end();)
	{
		if(new_D2.find(*iter)==new_D2.end())
		{
			iter = D2->erase(iter);
			if(D2->size()<support)
				return true;
		}
		else
		{
			iter++;
		}
	}

	return false;
}

/**
 * given the graph, check whether it is Acyclic or not (we assume the graph is connected)
 * @param me
 * @return
 */
bool GraMiCounter::isItAcyclic(GraphX& query)
{
	tr1::unordered_set<int> visited;
	vector<int> toBeVisited;

	int currentNodeID;
	NodeX* currentNode;
	toBeVisited.push_back(0);
	while(visited.size()<query.getNumOfNodes())
	{
		if(toBeVisited.size()==0)
			break;

		currentNodeID = toBeVisited.at(0);
		currentNode = query.getNodeWithID(currentNodeID);

		toBeVisited.erase(toBeVisited.begin());
		visited.insert(currentNodeID);
		//get all neighbor nodes (incoming and outgoing)
		int alreadyVisitedNeighbors = 0;//this should not be more than 1

		//all edges!
		for(tr1::unordered_map<int, void*>::iterator iter = currentNode->getEdgesIterator();iter!=currentNode->getEdgesEndIterator();iter++)
		{
			EdgeX* edge = (EdgeX*)iter->second;
			int otherNodeID = edge->getOtherNode()->getID();

			if(visited.find(otherNodeID)!=visited.end())
			{
				alreadyVisitedNeighbors++;
				if(alreadyVisitedNeighbors>1)
				{
					//cout<<"It is CYCLIC!";
					return false;
				}
			}
			else
			{
				toBeVisited.push_back(otherNodeID);
			}
		}
	}

	return true;
}

/**
 * Populate the ordering 'order' from domains with smaller sizes to those with larger sizes.
 */
void GraMiCounter::setDomainsOrder(tr1::unordered_map<int, tr1::unordered_set<int>*>& domains_values, int* order, Pattern* pattern)
{
	int invalidCol;
	if(Settings::usePredictedInvColumn)
		invalidCol = pattern->getInvalidCol();
	else
		invalidCol = -1;

	tr1::unordered_set<int> added;
	int i = 0;
	if(invalidCol>-1)
	{
		order[0] = invalidCol;
		added.insert(invalidCol);
		i = 1;
	}

	for(;i<domains_values.size();i++)
	{
		//find the minimum
		int min = std::numeric_limits<int>::max();
		int ID4min = -1;

		for(tr1::unordered_map<int, tr1::unordered_set<int>*>::iterator iter = domains_values.begin();iter!=domains_values.end();iter++)
		{
			//get the pattern node
			int currentID = iter->first;
			tr1::unordered_set<int>* currentSet = iter->second;

			long currentDomainValue = currentSet->size()/pattern->getGraph()->getNodeWithID(currentID)->getEdgesSize();

			if(currentDomainValue<min && added.find(currentID)==added.end())
			{
				ID4min = currentID;
				min = currentDomainValue;
			}
		}

		order[added.size()] = ID4min;
		added.insert(ID4min);
	}

	if(Settings::debugMSG)
	{
		for(int ii=0;ii<domains_values.size();ii++)
			cout<<order[ii]<<", ";
		cout<<endl;
	}
}
