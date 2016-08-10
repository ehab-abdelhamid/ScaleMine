/**
 * Copyright 2016 Ehab Abdelhamid, Ibrahim Abdelaziz

CanonicalLabel.cpp

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
#include<iostream>
#include<algorithm>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <tr1/unordered_set>
#include <sstream>

#include"utils.h"
#include"EdgeX.h"
#include"Settings.h"

long CanonicalLabel::elapsed1 = 0;
long CanonicalLabel::elapsed2 = 0;
char CanonicalLabel::sigBuffer[10000];

int mysprintf(char* str, int x)
{
	int i;
	if(x<10)
	{
		*str = x + '0';
		return 1;
	}
	else if(x<100)
	{
		*(str+1) = (x%10) + '0';
		*(str) = (x/10) + '0';
		return 2;
	}
	else if(x<1000) i=2;
	else if(x<10000) i=3;
	else
	{
		cout<<"Error x is too large: "<<x<<endl<<flush;
		exit(0);
	}
	int l = i+1;

	while(x > 0)
	{
		*(str+i) = (x % 10) + '0';
		x = x / 10;
		i = i - 1;
	}

	return l;
}

void printCombinations(vector<vector<NodeWithInfo* >* >* all)
{
	for(vector<vector<NodeWithInfo* >* >::iterator enum1 = all->begin();enum1!=all->end();++enum1)
	{
		vector<NodeWithInfo* >* combination = *enum1;
		for(vector<NodeWithInfo* >::iterator enum11 = combination->begin();enum11!=combination->end();++enum11)
		{
			NodeWithInfo* nInfo = *enum11;
			cout<<nInfo->node->getID()<<", ";
		}
		cout<<endl;
	}
}

/**
 * just for testing correctness
 */
bool CL_tester()
{
	tr1::unordered_map<string, void* > temp;

	cout<<"Canonical Label test started ..."<<endl;

	string g1 = "# t 1\nv 0 1\nv 1 1\nv 2 2\nv 3 2\nv 4 2\ne 1 2 2\ne 0 3 2\ne 3 4 1\ne 4 0 2";
	string g2 = "# t 1\nv 0 1\nv 2 2\nv 3 2\nv 4 2\ne 0 3 2\ne 3 4 1\ne 4 0 2";
	GraphX graph11;
	graph11.loadFromString(g1, temp);
	string CL11 = CanonicalLabel::generate(&graph11);
	GraphX graph22;
	graph22.loadFromString(g2, temp);
	string CL22 = CanonicalLabel::generate(&graph22);
	if(CL11.compare(CL22)==0)
		return false;

	//if(true) return true;

	//the following two should produce the same CL
	GraphX graph1_1;
	string str1 = "# t 1\nv 0 2\nv 1 1\nv 2 1\nv 3 2\ne 0 2 2\ne 0 3 2\ne 1 3 1\ne 1 2 1";
	graph1_1.loadFromString(str1, temp);
	string CL1_1 = CanonicalLabel::generate(&graph1_1);

	GraphX graph1_2;
	string str2 = "# t 1\nv 3 2\nv 1 1\nv 2 1\nv 0 2\ne 3 2 2\ne 3 0 2\ne 1 0 1\ne 1 2 1";
	graph1_2.loadFromString(str2, temp);
	string CL1_2 = CanonicalLabel::generate(&graph1_2);

	if(CL1_1.compare(CL1_2)!=0)
	{
		cout<<"Different!"<<endl<<CL1_1<<endl<<CL1_2<<endl;
		return false;
	}

	//the following two should produce the same CL
	GraphX graph1_3;
	string str3 = "# t 1\nv 0 1\nv 1 1\nv 2 1\nv 3 1\nv 4 1\nv 5 1\nv 6 1\nv 7 1\ne 0 1 1\ne 0 7 1\ne 0 6 1\ne 1 5 1\ne 1 2 1\ne 2 4 1\ne 2 3 1";
	graph1_3.loadFromString(str3, temp);
	string CL1_3 = CanonicalLabel::generate(&graph1_3);

	GraphX graph1_4;
	string str4 = "# t 1\nv 1 1\nv 0 1\nv 7 1\nv 4 1\nv 3 1\nv 5 1\nv 6 1\nv 2 1\ne 1 0 1\ne 1 2 1\ne 1 6 1\ne 0 5 1\ne 0 7 1\ne 7 3 1\ne 7 4 1";
	graph1_4.loadFromString(str4, temp);
	string CL1_4 = CanonicalLabel::generate(&graph1_4);

	if(CL1_3.compare(CL1_4)!=0)
	{
		cout<<"Different!"<<endl<<CL1_3<<endl<<CL1_4<<endl;
		return false;
	}

	return true;
}

string CL_Partition::getPrefix()
{
	string prefix;
	for(vector<NodeWithInfo*>::iterator iter = nodes.begin();iter!=nodes.end();iter++)
	{
		prefix.append(doubleToString((*iter)->node->getLabel()));
		prefix.append("\0");
	}

	return prefix;
}

void CL_Partition::CL_Part_Destructor()
{
	//delete combinations
	for(vector<vector<NodeWithInfo* >* >::iterator iter = combinations.begin();iter!=combinations.end();++iter)
	{
		vector<NodeWithInfo*>* tempV = (*iter);

		delete tempV;
	}

	st_counter--;
}

bool CanonicalLabel::sortPartitions(vector<CL_Partition* >& parts)
{
	bool b = true;
	sort(parts.begin(), parts.begin()+parts.size(), descending);//std::greater<CL_Partition*>());
	int i=0;
	for(vector<CL_Partition*>::iterator enum1 = parts.begin();enum1!=parts.end();enum1++)
	{
		CL_Partition* part = *enum1;
		int oldID = part->getID();
		if(oldID!=i)
			b = false;
		part->setID(i);
		i++;
	}
	return b;
}

/**
 * the main function for generating canonical label for the given graph
 * This function is based on the work published in:
 * "Finding Frequent Patterns in a Large Sparse Graph", DMKD 2005
 * If the graph is huge and takes much time, the returned canonical label is not unique and is identified by prefix of X
 */
string CanonicalLabel::generate(GraphX* graph)
{
	if(graph->getNumOfNodes()==0)
		return "";

	map<int , NodeWithInfo* > allNodes;

	//partition by label and degree
	map<string, CL_Partition* > parts;
	for( tr1::unordered_map<int,NodeX*>::const_iterator enum1=graph->getNodesIterator(); enum1!=graph->getNodesEndIterator(); ++enum1)
	{
		NodeX* node = (*enum1).second;
		NodeWithInfo* nInfo = new NodeWithInfo(node);
		string label = doubleToString(nInfo->node->getLabel());
		int degree = nInfo->node->getEdgesSize();
		string key = degree+"_"+label;

		map<string, CL_Partition* >::iterator tIter = parts.find(key);
		CL_Partition* part;
		if(tIter==parts.end())
		{
			part = new CL_Partition(parts.size());
			parts[key] = part;
		}
		else
		{
			part = (*tIter).second;
		}

		part->addNode(nInfo);
		nInfo->partID = part->getID();
		allNodes[nInfo->node->getID()] = nInfo;
		//cout<<"nInfo ID: "<<nInfo->node->getID()<<", nInfo->partID: "<<nInfo->partID<<endl;
		part->setSortingValues();
	}

	vector<CL_Partition* > partsV;
	MapToVec(parts, partsV);
	parts.clear();
	sortPartitions(partsV);
	//generate Neighbors list, then sort nodes inside each partition
	for(std::vector<CL_Partition*>::iterator enum1 = partsV.begin(); enum1 != partsV.end(); ++enum1)
	{
		CL_Partition* part = *enum1;
		//for each node generate NL
		for(std::vector<NodeWithInfo*>::const_iterator enum2 = part->getNodesEnum(); enum2 != part->getNodesEnd(); ++enum2)
		{
			NodeWithInfo* nInfo = *enum2;
			nInfo->nl = generateNeighborsList(nInfo, allNodes);
			//cout<<"NInfo Node ID: "<<nInfo->toString()<<endl;
		}
		//sort nodes
		part->sortNodes();
	}
	sortPartitions(partsV);

	//iterative partitioning
	while(true)
	{
		//generate new partitions
		bool b = false;
		unsigned int i = 0;
		for(vector<CL_Partition*>::iterator enum1 = partsV.begin(); enum1 != partsV.end(); ++enum1)
		{
			CL_Partition* part = *enum1;
			map<string, CL_Partition*>* newParts = part->getNewParts();
			if(newParts!= NULL)
			{
				//invalidate the 'enum1' iterator
				enum1 = partsV.end();

				vector<CL_Partition* > newPartsArr;
				MapToVec(*newParts, newPartsArr);
				newParts->clear();
				delete newParts;
				for(vector<CL_Partition* >::iterator j = newPartsArr.begin(); j!=newPartsArr.end(); ++j)
				{
					(*j)->setSortingValues();
				}

				sort(newPartsArr.begin(), newPartsArr.begin()+newPartsArr.size(), ascending);

				(*(partsV.begin()+i))->CL_Part_Destructor();
				delete (*(partsV.begin()+i));
				partsV.erase(partsV.begin()+i);

				int j_count = 0;
				for(vector<CL_Partition* >::iterator j = newPartsArr.begin(); j!=newPartsArr.end(); ++j,j_count++)
				{
					(*j)->setID(i+j_count);
					(*j)->setSortingValues();
					partsV.insert(partsV.begin()+i+j_count, *j);
				}

				for(unsigned int j=i+newPartsArr.size();j<partsV.size();j++)
				{
					(*(partsV.begin()+j))->setID(j);
				}
				newPartsArr.clear();

				//generate Neighbors list, then sort nodes inside each partition
				for(vector<CL_Partition*>::iterator enum11 = partsV.begin(); enum11 != partsV.end(); ++enum11)
				{
					CL_Partition* part2 = *enum11;

					int toto = part2->getID();

					//for each node generate NL
					for(vector<NodeWithInfo*>::const_iterator enum2 = part2->getNodesEnum();enum2!=part2->getNodesEnd();++enum2)
					{
						NodeWithInfo* nInfo = *enum2;
						nInfo->nl = generateNeighborsList(nInfo, allNodes);
					}
					//sort nodes
					part2->sortNodes();
				}

				b = true;

				break;
			}
			i++;
		}

		if(i==partsV.size()) break;
	}

	//reset counters and generate signature prefix
	string prefix;
	for(vector<CL_Partition* >::iterator i = partsV.begin();i!=partsV.end();++i)
	{
		if(enablePrint)
		{
			cout<<"Part:"+intToString((*i)->getID())<<endl;
			cout<<(*i)->getCombinations()->size()<<endl;
		}
		(*i)->counter = 0;
		prefix.append((*i)->getPrefix());
		prefix.append(",");
	}

	//generate a fast lookup array for edge labels
	//get the maximum label id
	int maxLabelID = 0;
	for(tr1::unordered_map<int, NodeX*>::const_iterator iter = graph->getNodesIterator();iter!=graph->getNodesEndIterator();++iter)
	{
		int currID = iter->second->getID();
		if(maxLabelID<currID)
			maxLabelID = currID;
	}

	double** ellt = new double*[maxLabelID+1];
	for(int i=0;i<(maxLabelID+1);i++)
		ellt[i] = new double[maxLabelID+1];
	//fill in the edge labels
	for(tr1::unordered_map<int, NodeX*>::const_iterator i=graph->getNodesIterator();i!=graph->getNodesEndIterator();i++)
	{
		int id1 = i->second->getID();
		for(tr1::unordered_map<int, NodeX*>::const_iterator ii=graph->getNodesIterator();ii!=graph->getNodesEndIterator();ii++)
		{
			int id2 = ii->second->getID();
			ellt[id1][id2] = graph->getEdgeLabel(id1, id2);
		}
	}

	for(int i=0;i<partsV.size();i++)
	{
		if(partsV.at(i)->getNumNodes()>10)
		{
			if(Settings::debugMSG)
			{
				cout<<"pre_pre_ number of combinations is very large = ("<<partsV.at(i)->getNumNodes()<<"!), "<<endl;
			}
			stringstream nosig;
			nosig<<"X"<<graph->getNumOfNodes()<<"_"<<graph->getNumOfEdges();
			return nosig.str();
		}
	}

	long numberOfCombinations_ = 1;
	for(int i=1;i<partsV.size();i++)
	{
		numberOfCombinations_ = numberOfCombinations_*partsV.at(i)->getCombinations()->size();
	}

	if(numberOfCombinations_>2000000 || partsV.at(0)->getCombinations()->size()>200000)
	{
		if(Settings::debugMSG)
		{
			cout<<"pre_ number of combinations exceeded: "<<numberOfCombinations_<<", "<<partsV.at(0)->getCombinations()->size()<<endl;
		}

		stringstream nosig;
		nosig<<"X"<<graph->getNumOfNodes()<<"_"<<graph->getNumOfEdges();
		return nosig.str();
	}

	//an optimization for the first partition
	//only store permutations where the first line is the max

	//get the max signature line
	char* maxSigLine = 0;
	CL_Partition* firstPart = partsV.at(0);
	for(vector<vector<NodeWithInfo*>*>::iterator iter = firstPart->getCombinations()->begin();iter!=firstPart->getCombinations()->end();iter++)
	{
		vector<NodeWithInfo*>* singleComb = (*iter);
		char* temp1 = generateCanLabel(singleComb, graph, ellt, true);

		if(maxSigLine==0 || strcmp(temp1, maxSigLine)>0)
		{
			if(maxSigLine!=0)
				delete[] maxSigLine;

			int length = strlen(temp1);
			maxSigLine = new char[length+1];
			strcpy(maxSigLine, temp1);
		}
	}

	//delete any permutation from the first partition that is less than the max
	for(vector<vector<NodeWithInfo*>*>::iterator iter = firstPart->getCombinations()->begin();iter!=firstPart->getCombinations()->end();)
	{
		vector<NodeWithInfo*>* singleComb = (*iter);

		char* temp1 = generateCanLabel(singleComb, graph, ellt, true);

		if(strcmp(temp1, maxSigLine)<0)
		{
			iter = firstPart->getCombinations()->erase(iter);
		}
		else
			iter++;
	}

	long numberOfCombinations = 1;
	for(int i=0;i<partsV.size();i++)
	{
		numberOfCombinations = numberOfCombinations*partsV.at(i)->getCombinations()->size();
	}

	if(numberOfCombinations>20000000)//numberOfCombinations>2)//it was 20000000
	{
		if(Settings::debugMSG)
		{
			cout<<"number of combinations exceeded: "<<numberOfCombinations_<<", "<<partsV.at(0)->getCombinations()->size()<<endl;
		}

		stringstream nosig;
		nosig<<"X"<<graph->getNumOfNodes()<<"_"<<graph->getNumOfEdges();
		return nosig.str();
	}

	//generate permutations and save the max
	char* maxCL = new char[1];
	maxCL[0] = 0;

	vector<NodeWithInfo* >* bestCombination = 0;
	bool b = true;
	int count = 0;
	while(true)
	{
		count++;
		if(count%1000000==0)
		{
			cout<<count<<"/"<<numberOfCombinations<<endl;
		}
		//add a combination
		vector<NodeWithInfo* >* singleCombination = new vector<NodeWithInfo* >();
		for(vector<CL_Partition* >::iterator i = partsV.begin();i!=partsV.end();++i)
		{
			vector<NodeWithInfo*>* temp = *(((*i)->getCombinations())->begin()+(*i)->counter);
			singleCombination->insert(singleCombination->end(), temp->begin(), temp->end());
		}

		char* temp = generateCanLabel(singleCombination, graph, ellt);

		if(enablePrint) cout<<temp<<endl;
		if(bestCombination==0 || strcmp(temp, maxCL)>0)
		{
			delete[] maxCL;
			delete bestCombination;

			int sigLength = strlen(temp);
			maxCL = new char[sigLength+1];
			strcpy (maxCL, temp);

			bestCombination = singleCombination;
		}
		else
		{
			//delete temp;
			delete singleCombination;
		}

		//increase counters
		int i_counter = partsV.size()-1;
		for(vector<CL_Partition* >::reverse_iterator i = partsV.rbegin();i!=partsV.rend();++i,--i_counter)
		{
			if((*i)->counter<(*i)->getCombinations()->size()-1)
			{
				(*i)->counter++;
				break;
			}
			else if(i_counter==0)
				b = false;
			else
				(*i)->counter = 0;
		}

		if(!b)
			break;
	}

	for(int i=0;i<(maxLabelID+1);i++)
		delete[] ellt[i];
	delete[] ellt;

	//cout<<"*"<<maxCL<<endl;

	if(enablePrint) cout<<"Start destruction ...."<<endl;
	for(vector<CL_Partition* >::iterator i = partsV.begin();i!=partsV.end();++i)
	{
		//added to delete NodeWithInfo for this partition
		if((*i)->getCombinations()->size()>0)
		{
			for(vector<NodeWithInfo*>::iterator iter1 = (*(*i)->getCombinations()->begin())->begin();iter1!=(*(*i)->getCombinations()->begin())->end();++iter1)
			{
				delete (*iter1);
			}
		}

		(*i)->CL_Part_Destructor();
		delete (*i);
	}
	partsV.clear();
	if(enablePrint) cout<<"All destruction work finished"<<endl;

	string temp(prefix);
	temp.append(maxCL);
	delete maxCL;
	return temp;
}

string CanonicalLabel::generateNeighborsList(NodeWithInfo* nInfo, map<int , NodeWithInfo* >& allNodes)
{
	string sb = "";

	NodeX* node = nInfo->node;

	vector<PartID_label* > orderedNL;
	for( tr1::unordered_map<int, void*>::const_iterator enum1=node->getEdgesIterator(); enum1!=node->getEdgesEndIterator(); ++enum1)
	{
		EdgeX* edge = (EdgeX*)enum1->second;

		//get the other node partition
		NodeX* oNode = edge->getOtherNode();
		int oNodePartID = (allNodes.find(oNode->getID()))->second->partID;

		//prepare the details
		string oNodeLabel = doubleToString(oNode->getLabel());
		PartID_label* pidl = new PartID_label();
		pidl->label = oNodeLabel;
		pidl->label = pidl->label+doubleToString(edge->getLabel());//added by ehab on 13 Aug 2015
		pidl->partID = oNodePartID;

		//put the details in order
		vector<PartID_label*>::iterator i = orderedNL.begin();
		for(;i!=orderedNL.end();i++)
		{
			PartID_label* currentpidl = *i;
			if(pidl->partID>currentpidl->partID || (pidl->partID==currentpidl->partID && pidl->label.compare(currentpidl->label)>0))
				;
			else
				break;
		}
		orderedNL.insert(i, pidl);
	}

	for(vector<PartID_label*>::iterator i=orderedNL.begin();i!=orderedNL.end();i++)
	{
		sb.append((*i)->toString());
		delete (*i);
	}

	return sb;
}

char* CanonicalLabel::generateCanLabel(vector<NodeWithInfo* >* nodes, GraphX* graph, double** ellt, bool oneLineOnly)
{
	char* sb = CanonicalLabel::sigBuffer;//new char[(nodes->size()-1)*10];
	char* startSB = sb;

	//add the upper triangle thing
	vector<NodeWithInfo* >::iterator nodesEnd = nodes->end();
	vector<NodeWithInfo* >::iterator j;
	vector<NodeWithInfo* >::iterator enum1 = nodes->begin();
	while(true)
	//for(vector<NodeWithInfo* >::iterator enum1 = nodes->begin();enum1!=nodesEnd;++enum1)
	{
		NodeWithInfo* nInfo1 = *enum1;
		double* _label = ellt[nInfo1->node->getID()];

		j = enum1+1;
		if(j==nodesEnd)
			break;
		for(;j!=nodesEnd;++j)
		{
			double label = _label[(*j)->node->getID()];
			if(label==0.0001)
			{
				*sb='0';
				sb++;
			}
			else
			{
				sb+=mysprintf(sb, (int)label);
			}
		}
		*sb=',';
		sb++;

		if(oneLineOnly)
			break;
		++enum1;
	}

	*sb = '\0';
	return startSB;
}

bool ascending(CL_Partition* a, CL_Partition* b)
{
	int r = a->nodeLabel.compare(b->nodeLabel);
	if(r!=0)
	{
		if(r<0) return true;
		else return false;
	}

	if(a->degree!=b->degree)
	{
		if(a->degree-b->degree<0) return true;
		else return false;
	}

	r = a->nl.compare(b->nl);
	if(r!=0)
	{
		if(r<0) return true;
		else return false;
	}

	if(a->partID<b->partID)
		return true;
	else
		return false;
}

bool descending(CL_Partition* a, CL_Partition* b)
{
	return !ascending(a, b);
}

/**
 * this function should be called once all partitions get properly partitioned
 */
void CL_Partition::setSortingValues()
{
	NodeWithInfo* nInfo = *(nodes.begin());
	degree = nInfo->node->getEdgesSize();
	nodeLabel = nInfo->node->getLabel();
	nl = nInfo->nl;
}

CL_Partition::CL_Partition(int partID)
{
	this->partID = partID;
	CL_Partition::st_counter++;
}

void CL_Partition::addNode(NodeWithInfo* nInfo)
{
	nInfo->partID = this->partID;
	nodes.push_back(nInfo);
}

vector<NodeWithInfo*>::const_iterator CL_Partition::getNodesEnum() const
{
	return nodes.begin();
}

int CL_Partition::getNumNodes()
{
	return nodes.size();
}

void CL_Partition::clearNodes()
{
	nodes.clear();
}

void CL_Partition::setID(int partID)
{
	this->partID = partID;
	for(vector<NodeWithInfo*>::iterator enum1 = nodes.begin();enum1!=nodes.end();++enum1)
	{
		NodeWithInfo* nInfo = *enum1;
		nInfo->partID = this->partID;
	}
}

int CL_Partition::getID()
{
	return partID;
}

/**
 * get the partitioning of this partition
 * returns null if no need to partition
 * [Translated to C++]
 */
map<string, CL_Partition*>* CL_Partition::getNewParts()
{
	if(nodes.size()==1)
		return NULL;

	map<string, CL_Partition* >* parts = new map<string, CL_Partition* >();
	for(vector<NodeWithInfo*>::iterator enum1 = nodes.begin(); enum1!=nodes.end(); ++enum1)
	{
		NodeWithInfo* nInfo = *enum1;
		string key = nInfo->nl;

		map<string, CL_Partition* >::iterator temp = parts->find(key);
		CL_Partition* part;
		if(temp==parts->end())
		{
			part = new CL_Partition(parts->size());
			(*parts)[key]= part;
		}
		else
			part = temp->second;
		part->addNode(nInfo);
	}

	if(parts->size()==1)
	{
		parts->begin()->second->CL_Part_Destructor();
		delete parts->begin()->second;
		parts->clear();
		delete parts;
		return NULL;
	}
	else
		return parts;
}

string CL_Partition::toString()
{
	string sb = "";
	sb.append("Partition ID: "+intToString(partID)+"\n");

	for(vector<NodeWithInfo*>::iterator enum1 = nodes.begin(); enum1!=nodes.end(); ++enum1)
	{
		NodeWithInfo* nInfo = *enum1;
		sb.append(nInfo->toString()+"\n");
	}

	return sb;
}

vector<vector<NodeWithInfo* >* >* CL_Partition::getCombinations()
{
	if(combinations.size()>0)
		return &combinations;

	vector<NodeWithInfo* >* notused = new vector<NodeWithInfo* >(nodes);

	//check if all nodes within this partition are only connected to each other
	//collect all nodes in a map
	tr1::unordered_set<int> allNodes;
	for(vector<NodeWithInfo*>::iterator iter = nodes.begin();iter!=nodes.end();iter++)
	{
		allNodes.insert((*iter)->node->getID());
	}
	bool connectedWithin = true;
	for(vector<NodeWithInfo*>::iterator iter = nodes.begin();iter!=nodes.end();iter++)
	{
		NodeX* node = (*iter)->node;
		for(tr1::unordered_map<int, void*>::iterator iter1 = node->getEdgesIterator();iter1!=node->getEdgesEndIterator();iter1++)
		{
			int otherNodeID = ((EdgeX*)iter1->second)->getOtherNode()->getID();
			if(allNodes.find(otherNodeID)==allNodes.end())
			{
				connectedWithin = false;
				break;
			}
		}
		if(!connectedWithin)
			break;
	}

	combinations_fn(notused, connectedWithin);

	delete notused;

	return &combinations;
}

void CL_Partition::combinations_fn(vector<NodeWithInfo* >* notused, bool connectedWithin)
{
	do{
		vector<NodeWithInfo* >* v1 = new vector<NodeWithInfo* >(*notused);
		combinations.insert(combinations.end(), v1);
		//this needs a proof
		if(connectedWithin)
			break;
	} while(std::next_permutation(notused->begin(),notused->end(), NodeWithInfo_descending));
}

void CL_Partition::sortNodes()
{
	sort(nodes.begin(), nodes.begin()+nodes.size(), NodeWithInfo_descending);//, std::greater<NodeWithInfo>());
}

bool NodeWithInfo_ascending (NodeWithInfo* a, NodeWithInfo* b)
{
	int r = a->nl.compare(b->nl);
	if(r<0)
		return true;
	else if(r==0)//this part is not needed for the algorithm, but it is needed for the STL::next_permutation function
	{
		if(a->node->getID()<b->node->getID())
			return true;
		else
			return false;
	}
	else
		return false;
}
bool NodeWithInfo_descending (NodeWithInfo* a, NodeWithInfo* b)
{
	return !NodeWithInfo_ascending(a, b);
}

NodeWithInfo::NodeWithInfo(NodeX* node)
{
	this->node = node;
	NodeWithInfo::st_counter++;
}

string NodeWithInfo::toString()
{
	return "NodeID: "+intToString(node->getID())+"\nPartID: "+intToString(partID)+"\nNeighbors List:"+nl;
}

string PartID_label::toString()
{
	return "(p"+intToString(partID)+","+label+"),";
}
