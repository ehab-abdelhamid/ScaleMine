/**
 * Copyright 2016 Ehab Abdelhamid, Ibrahim Abdelaziz

Miner.cpp

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
 * the base mining class for the master
 */

#include <iostream>
#include <sstream>
#include<sys/time.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include "Miner.h"
#include "utils.h"

long long getmsofday_2()
{
   struct timeval tv;
   struct timezone tz;
   gettimeofday(&tv, &tz);
   return (long long)tv.tv_sec*1000 + tv.tv_usec/1000;
}

int Miner::numIsFreqCalls = 0;
int Miner::numFreqs = 0;
int Miner::numInfreq = 0;

Miner::Miner()
{
	candidates.insert(candidates.begin(), new CLMap());//graphs with zero edges!
	candidates.insert(candidates.begin()+1, new CLMap());//graphs with one edge! no ned also, because we have special treatment for graphs with one edge

	frequentPatterns.insert(frequentPatterns.begin(), new CLMap());//frequent patterns with zero edges!

	infrequentPatterns.insert(infrequentPatterns.begin(), new CLMap());//infrequents with zero edges!
	infrequentPatterns.insert(infrequentPatterns.begin()+1, new CLMap());//infrequents with one edge!
}

/**
 * start the mining process, given the filename (file base name for the partitions)
 * , and the required support threshold
 */
void Miner::initMining(string fileName, int graphType, int support, int numWorkers)
{
	long long start1 = getmsofday_2();

	//first time: send task to all workers
	cout<<"Sending graph filename to all workers: "<<endl;
	for(int i=1;i<numWorkers;i++)
	{
		//Graph Loading Task
		//produceTask(rank, size, str_message);
		int frequency = support;
		char fileName_[500];
		sprintf(fileName_, "f%s,%d\n", fileName.c_str(), frequency);
		int cnt=strlen(fileName_)+1;
		MPI_Send(fileName_, cnt, MPI_CHAR, i, 0, MPI_COMM_WORLD);
		//MPI_Send(&frequency, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		cout<<"Sent to worker: "<<i<<endl;
	}
	cout<<"Finished sending to all workers."<<endl;

	this->support = support;

	//load the graph
	long long start = getmsofday_2();
	loadGraph(fileName, graphType, support, frequentEdges);
	long long end = getmsofday_2();
	graphLoadingElapsed = end - start;

	masterProcessingTime += (getmsofday_2() - start1);
}

void Miner::startMining(int numWorkers)
{
	long long start1 = getmsofday_2();

	//create available workers list
	availableWorker = new bool[numWorkers];
	for(int i=1;i<numWorkers;i++)
	{
		availableWorker[i] = true;
	}

	//create variables for graph loading time, processing time for each worker
	graphloadingTime = new unsigned long[numWorkers];
	processingTime = new unsigned long[numWorkers];
	for(int i=1;i<numWorkers;i++)
	{
		graphloadingTime[i] = 0;
		processingTime[i] = 0;
	}

	CLMap* clmap = new CLMap();
	clmap->addAll(&frequentEdges);
	frequentPatterns.insert(frequentPatterns.begin()+1, clmap);//frequent edges!

	cout<<"Frequent edges:"<<endl;
	frequentEdges.print();

	//extend frequent edges into 2 edges graphs
	//cout<<"Extending frequent edges"<<endl;
	extendFreqEdges();

	//Phase 1: loop over candidates, check frequentness, extend if frequent and add to candidate

	masterProcessingTime += (getmsofday_2() - start1);

	unsigned int smallestCandidate = 2;
	while(true)
	{
		start1 = getmsofday_2();
		//check if we have more candidates with the same current size
		while(candidates.at(smallestCandidate)->getSize()==0)
		{
			smallestCandidate++;
			//if we finished all candidates, break
			if(smallestCandidate>=candidates.size())
				break;
		}

		masterProcessingTime += (getmsofday_2() - start1);

		cout<<"candidates.size = "<<candidates.size()<<", smallestCandidate = "<<smallestCandidate<<endl;
		if(smallestCandidate==candidates.size())
			break;
		//start parallelizing
		CLMap* currentCandidates = candidates.at(smallestCandidate);

		cout<<"*******************"<<endl;

		//send the first set of subgraphs to workers
		//at this step, i am sure no wroker is running
		cout<<currentCandidates->getSize()<<endl;
		for(int i=1;i<numWorkers;i++)
		{
			popACandidate(currentCandidates, currentlyChecking, i);
		}

		//go over all candidates in this list
		while(true)
		{
			MPI_Status status;
			char str_message[2000];
			MPI_Recv(str_message, 2000, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			//ack message recieved
			if(str_message[0]=='l')
			{
				start1 = getmsofday_2();

				printf("Master received ack from %d\n", status.MPI_SOURCE);
				//get elapsed tie for loading the graph at the worker
				unsigned long time;
				sscanf(str_message, "l,%Lu", &time);
				graphloadingTime[status.MPI_SOURCE] = time;

				masterProcessingTime += (getmsofday_2() - start1);

				availableWorker[status.MPI_SOURCE] = true;
			}
			else if(str_message[0]=='t')
			{
				start1 = getmsofday_2();
				Miner::numFreqs++;

				//get the pattern that was sent
				int candidID;
				unsigned long time;
				sscanf(str_message, "t%d,%Lu", &candidID, &time);
				Pattern* candidate = currentlyChecking.find(candidID)->second;
				candidate->setFrequency(support);
				currentlyChecking.erase(candidID);

				char temp[500];
				sprintf(temp, "%d_%d", status.MPI_SOURCE, candidID);
				isFrerqtimes.insert(std::pair<string, unsigned long>(temp, time));

				printf("Master received isFreq[TRUE] from %d\n", status.MPI_SOURCE);

				cout<<"Found to be FREQUENT, with frequency = "<<candidate->getFrequency()<<endl;

				if((frequentPatterns.size()-1)<candidate->getSize())
				{
					cout<<"HHHH-1"<<endl;
					cout<<"candidate->getSize() = "<<candidate->getSize()<<endl;
					cout<<"frequentPatterns.size() = "<<frequentPatterns.size()<<endl;
					frequentPatterns.insert(frequentPatterns.begin()+candidate->getSize(), new CLMap());
				}
				cout<<"HHHH-2"<<endl;
				frequentPatterns[candidate->getSize()]->addPattern(candidate);

				//extend using the current frequent patterns with the same size
				CLMap_Iterator clmapIter = frequentPatterns[candidate->getSize()]->getFirstElement();
				while(clmapIter.pattern!=0)
				{
					Pattern* temp = clmapIter.pattern;
					frequentPatterns[candidate->getSize()]->advanceIterator(clmapIter);

					set<std::pair<PrimaryGraph*, PrimaryGraph* > > joiningPGs = candidate->getJoiningPG(temp);
					for(set<std::pair<PrimaryGraph*, PrimaryGraph* > >::iterator pgsIter = joiningPGs.begin();pgsIter!=joiningPGs.end();++pgsIter)
					{
						PrimaryGraph* candPG = pgsIter->first;
						PrimaryGraph* otherPG = pgsIter->second;

						//find isomorphisms
						vector<map<int, int>* > result;
						tr1::unordered_map<int, tr1::unordered_set<int>*> domains_values;
						otherPG->getGraph()->isIsomorphic(candPG->getGraph(), result, domains_values);

						//for each result
						for(vector<map<int, int>* >::iterator iter2 = result.begin();iter2!=result.end();++iter2)
						{
							map<int, int>* currResult = *iter2;

							//add the new edge
							int scrNodeID = currResult->find(otherPG->getSrcNodeID())->second;
							EdgeX* edge = otherPG->getEdge();
							GraphX* newGraph = new GraphX(candidate->getGraph());
							int newID = newGraph->getNumOfNodes();
							while(newGraph->getNodeWithID(newID)!=NULL)
							{
								newID++;
							}
							newGraph->AddNode(newID, edge->getOtherNode()->getLabel());
							newGraph->addEdge(scrNodeID, newID, edge->getLabel());
							//create the new candidate object
							Pattern* newCandidate = new Pattern(newGraph, false);

							if(processed.getPattern(newCandidate)==0)
								;
							else
							{
								cout<<newCandidate->getGraph()->getID()<<" not added (1)!!!!!!!!!"<<endl;
								continue;
							}

							//remove the same subgraph from expected frequent candidates if exists
							if(expected_frequentPatterns)
							{
								bool b = false;
								for(int i=0;i<expected_frequentPatterns->size();i++)
								{
									if(expected_frequentPatterns->at(i)->getPattern(newCandidate)!=0)
									{
										cout<<newCandidate->getGraph()->getCanonicalLabel()<<" not added (2)!!!!!!!!!"<<endl;
										b = true;
										break;
									}
								}
								if(b)
									continue;
							}
							//remove the same subgraph from expected infrequent candidates if exists
							if(expected_infrequentPatterns)
							{
								bool b = false;
								for(int i=0;i<expected_infrequentPatterns->size();i++)
								{
									if(expected_infrequentPatterns->at(i)->getPattern(newCandidate)!=0)
									{
										cout<<newCandidate->getGraph()->getCanonicalLabel()<<" not added (3)!!!!!!!!!"<<endl;
										b = true;
										break;
									}
								}
								if(b)
									continue;
							}

							if(candidates.size()<=newCandidate->getSize())
							{
								for(int i=candidates.size();i<=newCandidate->getSize();i++)
								{
									candidates.insert(candidates.begin()+i, new CLMap());
								}
							}
							newCandidate->generatePrimaryGraphs();
							bool b = candidates.at(newCandidate->getSize())->addPattern(newCandidate);
							if(!b) delete newCandidate;
							if(candPG->getEdge()->getOtherNode()->getLabel()==otherPG->getEdge()->getOtherNode()->getLabel() &&
									candPG->getSrcNodeID()!=scrNodeID)
							{
								newGraph = new GraphX(candidate->getGraph());
								newGraph->addEdge(candPG->getEdge()->getOtherNode()->getID(), scrNodeID, otherPG->getEdge()->getLabel());
								newCandidate = new Pattern(newGraph, false);
								newCandidate->generatePrimaryGraphs();
								//CHECK THIS LINE IS TAKING SO MUCH TIME WHEN THE PATTERN IS BIG!!!! ---->
								bool b = candidates.at(newCandidate->getSize())->addPattern(newCandidate);
								if(!b) delete newCandidate;
							}
							delete *iter2;
						}
					}
				}

				availableWorker[status.MPI_SOURCE] = true;

				masterProcessingTime += (getmsofday_2() - start1);
			}
			else if(str_message[0]=='f')
			{
				start1 = getmsofday_2();

				Miner::numInfreq++;

				//get the pattern that was sent
				int candidID;
				unsigned long time;
				sscanf(str_message, "f%d,%Lu", &candidID, &time);
				Pattern* candidate = currentlyChecking.find(candidID)->second;
				candidate->setFrequency(0);
				currentlyChecking.erase(candidID);

				char temp[500];
				sprintf(temp, "%d_%d", status.MPI_SOURCE, candidID);
				isFrerqtimes.insert(std::pair<string, unsigned long>(temp, time));

				while(infrequentPatterns.size()<=candidate->getSize())
				{
					infrequentPatterns.push_back(new CLMap());
				}
				infrequentPatterns[candidate->getSize()]->addPattern(candidate);

				availableWorker[status.MPI_SOURCE] = true;

				masterProcessingTime += (getmsofday_2() - start1);
			}

			for(int i=1;i<numWorkers;i++)
			{
				popACandidate(currentCandidates, currentlyChecking, i);//status.MPI_SOURCE);
			}

			if(currentCandidates->getSize()==0 && currentlyChecking.size()==0)
			{
				break;
			}
		}
		//end parallelizing

		cout<<"---------------------------------------------"<<endl;
	}

	printResult();

	cout<<"Number of calls to IsFreq = "<<Miner::numIsFreqCalls<<endl;
	cout<<"Number of inserted candids = "<<this->numInsertedCandids<<endl;
	//printing time per isFreq
	long long allTime = 0;
	cout<<"time for each isFreq call"<<endl;
	for(map<string, long long>::iterator iter1 = isFrerqtimes.begin();iter1!=isFrerqtimes.end();++iter1)
	{
		string str = iter1->first;

		int workerID;
		int candidID;
		sscanf(str.c_str(), "%d_%d", &workerID, &candidID);

		long long time = iter1->second;
		allTime+=time;
		cout<<str<<":"<<time<<endl;

		processingTime[workerID] += time;
	}

	//printing processing time per worker
	for(int i=1;i<numWorkers;i++)
	{
		cout<<"Worker ID: "<<i<<", graph loading took time: "<<graphloadingTime[i]<<endl;
		cout<<"Worker ID: "<<i<<", processing took time: "<<processingTime[i]<<endl;
	}

	cout<<"Overall processing time by workers: "<<allTime<<endl;
	cout<<"Master processing time: "<<masterProcessingTime<<endl;

	cout<<"Graph loading [master] took "<<(graphLoadingElapsed/1000)<<" sec and "<<(graphLoadingElapsed%1000)<<" ms"<<endl;

	delete availableWorker;
	delete graphloadingTime;
	delete processingTime;
}

void Miner::loadGraph(string baseName, int graphType, int support, CLMap& freqEdges)
{
	this->support = support;

	tr1::unordered_map<string, void* > allMPs;

	int nodesCounter = 0;

	//load graph data
	std::stringstream sstmGF;
	sstmGF << baseName;
	string partFileName = sstmGF.str();
	graph = new GraphX(1, graphType);//create a graph with id=its partition id
	graph->setFrequency(support);
	cout<<"Master: created a graph object"<<endl;
	if(!graph->loadFromFile(partFileName, allMPs))
	{
		cout<<"Master: delete a graph object"<<endl;
		delete graph;
	}
	nodesCounter+=graph->getNumOfNodes();

	cout<<"Master: loop starts"<<endl;
	for(tr1::unordered_map<string, void* >::iterator iter = allMPs.begin(); iter!=allMPs.end();++iter)
	{
		if(((Pattern*)((*iter).second))->getFrequency()>=support)
			freqEdges.addPattern((Pattern*)((*iter).second));
		else
			delete (Pattern*)((*iter).second);
	}
	cout<<"Master: loop finishes"<<endl;

	cout<<" data loaded successfully!"<<endl;
}

void Miner::popACandidate(CLMap* currentCandidates, map<int, Pattern*>& currentlyChecking, int destination)
{
	if(!availableWorker[destination])
		return;

	cout<<"Worker "<<destination<<" is waiting"<<endl;
	if(expected_frequentPatterns)
	{
		cout<<"#candids = "<<currentCandidates->getSize()<<", #expectedFreq = "<<getNumElems((vector<map<string, void*>* >*)expected_frequentPatterns)<<", #expectedInfreq = "<<getNumElems((vector<map<string, void*>* >*)expected_infrequentPatterns)<<endl;
	}

	//get a candidate, if there is
	if(currentCandidates->getSize()>0)
	{
		CLMap_Iterator iter = currentCandidates->getFirstElement();
		string key = iter.key;
		Pattern* candidate = iter.pattern;
		currentCandidates->removePattern(candidate);
		sendACandidate(key, candidate, currentlyChecking, destination);

		//remove it from expected_frequents
		if(expected_frequentPatterns)
		{
			for(int i=0;i<expected_frequentPatterns->size();i++)
			{
				expected_frequentPatterns->at(i)->removePattern(candidate);
			}
		}
		//remove it from expected_infrequents
		if(expected_infrequentPatterns)
		{
			for(int i=0;i<expected_infrequentPatterns->size();i++)
			{
				expected_infrequentPatterns->at(i)->removePattern(candidate);
			}
		}
	}
	else
	{
		Pattern* expCandid = 0;
		string key;
		if(expected_frequentPatterns && expected_frequentPatterns->size()>0)
		{
			CLMap* expCandidates = *(expected_frequentPatterns->begin());
			while(expCandidates->getSize()==0)
			{
				expected_frequentPatterns->erase(expected_frequentPatterns->begin());
				if(expected_frequentPatterns->size()==0)
				{
					expCandidates = 0;
					break;
				}
				expCandidates = *(expected_frequentPatterns->begin());
			}

			if(expCandidates)
			{
				CLMap_Iterator tempIter = expCandidates->getFirstElement();
				key = tempIter.key;
				expCandid = tempIter.pattern;
				expCandidates->removePattern(expCandid);

				if(expCandidates->getSize()==0)
					expected_frequentPatterns->erase(expected_frequentPatterns->begin());
			}
		}

		if(expected_infrequentPatterns && expected_infrequentPatterns->size()>0 && expCandid==0)
		{
			CLMap* expCandidates = *(expected_infrequentPatterns->begin());
			while(expCandidates->getSize()==0)
			{
				expected_infrequentPatterns->erase(expected_infrequentPatterns->begin());
				if(expected_infrequentPatterns->size()==0)
				{
					expCandidates = 0;
					break;
				}
				expCandidates = *(expected_infrequentPatterns->begin());
			}

			if(expCandidates)
			{
				CLMap_Iterator tempIter = expCandidates->getFirstElement();
				key = tempIter.key;
				expCandid = tempIter.pattern;
				expCandidates->removePattern(expCandid);

				if(expCandidates->getSize()==0)
					expected_infrequentPatterns->erase(expected_infrequentPatterns->begin());
			}
		}

		if(expCandid)
		{
			expCandid->makeIDNegative();
			sendACandidate(key, expCandid, currentlyChecking, destination);
		}
	}
}

void Miner::sendACandidate(string key, Pattern* candidate, map<int, Pattern*>& currentlyChecking, int destination)
{
	cout<<"Worker "<<destination<<" will be busy!"<<endl;
	long long start = getmsofday_2();
	cout<<"Check frequency for:"<<endl;
	cout<<"CL:"<<(*(candidate->getGraph())).getCanonicalLabel()<<endl;
	cout<<*(candidate->getGraph())<<endl<<flush;

	//send isFreq request for the current candidate to worker
	currentlyChecking.insert(std::pair<int, Pattern*>(candidate->getID(), candidate));

	ostringstream tempOS;
	tempOS<<*(candidate->getGraph());
	char graphStr[2000];
	sprintf(graphStr, "s,%d,%s\t", candidate->getID(),tempOS.str().c_str());
	int cnt=strlen(graphStr)+1;
	MPI_Status status;
	MPI_Send(graphStr, cnt, MPI_CHAR, destination, 0, MPI_COMM_WORLD);
	printf("Sending subgraph[%d] to worker: %d for frequency check\n", candidate->getID(),destination);
	printf("%s\n", graphStr);

	processed.addPattern(candidate);
	availableWorker[destination] = false;

	Miner::numIsFreqCalls++;
}

void Miner::extendFreqEdges()
{
	CLMap* twoEdgesCandidate = new CLMap();
	candidates.insert(candidates.begin()+2, twoEdgesCandidate);

	long long totalElapsed = 0;

	CLMap_Iterator iter1 = frequentEdges.getFirstElement();
	while(iter1.pattern!=0)
	//for(map<string, Pattern*>::iterator iter1 = frequentEdges.begin();iter1!=frequentEdges.end();++iter1)
	{
		Pattern* edge1 = iter1.pattern;//->second;

		double l1_0 = edge1->getGraph()->getNodeWithID(0)->getLabel();
		double l1_1 = edge1->getGraph()->getNodeWithID(1)->getLabel();

		CLMap_Iterator iter2 = iter1.getCopy();
		while(iter2.pattern!=0)
		//for(map<string, Pattern*>::iterator iter2 = iter1;iter2!=frequentEdges.end();++iter2)
		{
			Pattern* edge2 = iter2.pattern;//->second;
			frequentEdges.advanceIterator(iter2);

			double l2_0 = edge2->getGraph()->getNodeWithID(0)->getLabel();
			double l2_1 = edge2->getGraph()->getNodeWithID(1)->getLabel();
			double edge2Label = edge2->getGraph()->getEdgeLabel(0,1);

			//connectType refers to the way they can connect, -q means they can not connect, 0 means 0 and 0, 1 means 0 and 1, 2 means 1 and 0, and 3 means 1 and 1
			int lastConnectType = -1;
			while(true)
			{
				int connectType = -1;

				//check how can they (edge1 and edge2) connect.
				if(l1_0==l2_0 && lastConnectType<0)
					connectType = 0;
				else if(l1_0==l2_1 && lastConnectType<1)
					connectType = 1;
				else if(l1_1==l2_0 && lastConnectType<2)
					connectType = 2;
				else if(l1_1==l2_1 && lastConnectType<3)
					connectType = 3;

				//could not find an extension by edge1 and edge2
				if(connectType==-1)
					break;

				lastConnectType = connectType;

				long long start = getmsofday_2();

				Pattern* candidate = new Pattern(edge1);
				candidate->invalidateFrequency();

				long long end = getmsofday_2();
				totalElapsed += (end-start);

				switch(connectType)
				{
				case 0://00
					candidate->extend(0, 2, l2_1, edge2Label);
					break;
				case 1://01
					candidate->extend(0, 2, l2_0, edge2Label);
					break;
				case 2://10
					candidate->extend(1, 2, l2_1, edge2Label);
					break;
				case 3://11
					candidate->extend(1, 2, l2_0, edge2Label);
					break;
				}

				if(twoEdgesCandidate->getPattern(candidate)!=0)
					delete candidate;
				else
				{
					candidate->generatePrimaryGraphs();

					bool b = twoEdgesCandidate->addPattern(candidate);
					if(!b) delete candidate;
				}
			}

		}

		frequentEdges.advanceIterator(iter1);
	}
}

void Miner::setFrequentEdges(CLMap& freqEdges)
{
	this->frequentEdges.addAll(&freqEdges);
}

void Miner::printCandidates()
{
	print(candidates);
}

void Miner::printResult()
{
	print(frequentPatterns);
}

void Miner::removePattern(Pattern* pattern, vector<CLMap* >& data)
{
	CLMap* tempList = data.at(pattern->getSize());
	tempList->removePattern(pattern);
}

void Miner::print(vector<CLMap* >& data)
{
	//count the total number of frequent patterns
	int size = 0;
	for(unsigned int i=0;i<data.size();i++)
	{
		size+=data[i]->getSize();
	}
	int count = 1;
	cout<<"[Miner] There are "<<size<<" frequent patterns, and they are:"<<endl;
	for(unsigned int i=0;i<data.size();i++)
	{
		cout<<"With "<<(i)<<" edges:"<<endl;
		data[i]->print();
	}
}

Miner::~Miner()
{
	if(graph!=0)
		delete graph;

	vect_map_destruct(frequentPatterns);
	vect_map_destruct(infrequentPatterns);
	vect_map_destruct(candidates);
}
