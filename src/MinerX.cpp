/**
 * Copyright 2016 Ehab Abdelhamid, Ibrahim Abdelaziz

MinerX.cpp

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

#include <iostream>
#include <sstream>
#include<sys/time.h>
#include <vector>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "MinerX2.h"
#include "utils.h"
#include "Settings.h"

int MinerX::numIsFreqCalls = 0;
int MinerX::numFreqs = 0;
int MinerX::numInfreq = 0;

MinerX::MinerX(int numWorkers, int nThreads)
{
	this->numWorkers = numWorkers * nThreads;

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
void MinerX::initMining(string fileName, int graphType, int support, int numWorkers)
{
	long long start1 = getmsofday();

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
		cout<<"Sent to worker: "<<i<<endl;
	}
	cout<<"Finished sending to all workers."<<endl;

	this->support = support;

	//set the frequency counter module
	//load the graph
	long long start = getmsofday();
	loadGraph(fileName, graphType, support, frequentEdges);
	long long end = getmsofday();
	graphLoadingElapsed = end - start;

	masterProcessingTime += (getmsofday() - start1);
}

void MinerX::startMining()
{
	long long start1 = getmsofday();

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

	cout<<"#Frequent edges = "<<frequentEdges.getSize()<<endl;
	cout<<"Frequent edges:"<<endl;
	frequentEdges.print();

	//extend frequent edges into 2 edges graphs
	cout<<"Extending frequent edges"<<endl;

	long long start = getmsofday();

	extendFreqEdges();

	long long end = getmsofday();

	masterProcessingTime += (getmsofday() - start1);
	//Start of New loop
	while(true)
	{
		//do we have busy workers?
		bool doWeHaveABusyWorker = false;
		for(int i=1;i<numWorkers;i++)
		{
			if(!availableWorker[i])
			{
				doWeHaveABusyWorker = true;
				break;
			}
		}

		cout<<doWeHaveABusyWorker<<endl;

		//get a message from a worker
		if(doWeHaveABusyWorker)
		{
			MPI_Status status;
			char str_message[2000];
			MPI_Recv(str_message, 2000, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			start1 = getmsofday();

			printf("The following message is recieved at master: %s", str_message);

			//ack message recieved
			if(str_message[0]=='l')
			{
				start1 = getmsofday();

				printf("Master received ack from %d\n", status.MPI_SOURCE);
				//get elapsed tie for loading the graph at the worker
				unsigned long time;
				sscanf(str_message, "l,%Lu", &time);
				graphloadingTime[status.MPI_SOURCE] = time;

				masterProcessingTime += (getmsofday() - start1);

				availableWorker[status.MPI_SOURCE] = true;
			}
			else if(str_message[0]=='t')
			{
				start1 = getmsofday();
				MinerX::numFreqs++;

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

				patternFoundFrequent(candidate);

				availableWorker[status.MPI_SOURCE] = true;

				masterProcessingTime += (getmsofday() - start1);
			}
			else if(str_message[0]=='f')
			{
				start1 = getmsofday();

				MinerX::numInfreq++;

				//get the pattern that was sent
				int candidID;
				unsigned long time;
				sscanf(str_message, "f%d,%Lu", &candidID, &time);
				Pattern* candidate = currentlyChecking.find(candidID)->second;
				candidate->setFrequency(0);
				cout<<"Predicted Time for candidate["<<candidate->getID()<<"] = '"<<candidate->getPredictedTime()<<"'"<<endl;
				currentlyChecking.erase(candidID);

				char temp[500];
				sprintf(temp, "%d_%d", status.MPI_SOURCE, candidID);
				isFrerqtimes.insert(std::pair<string, unsigned long>(temp, time));

				patternFoundInFrequent(candidate);

				availableWorker[status.MPI_SOURCE] = true;

				masterProcessingTime += (getmsofday() - start1);
			}//the case when a task is partitioned
			else if(str_message[0]=='m')
			{
				start1 = getmsofday();

				//get the pattern that was sent
				int candidID;
				unsigned long time;
				//message format: m(candidateid),(elapsed),mniCol[0],mniCol[1],...,mniCol[n],0
				char tempStr[100];
				sscanf(str_message, "m%d,%Lu,%s", &candidID, &time, tempStr);

				int c = 0;
				while(true)
				{
					char temp[500];
					sprintf(temp, "%d_%d_%d", status.MPI_SOURCE, candidID, c);
					if(isFrerqtimes.find(temp)==isFrerqtimes.end())
					{
						isFrerqtimes.insert(std::pair<string, unsigned long>(temp, time));
						break;
					}
					c++;
				}

				Pattern* candidate = currentlyChecking.find(candidID)->second;

				//populate the mniTable
				int* mniTable = new int[candidate->getGraph()->getNumOfNodes()];
				char* pch = strtok (tempStr,",");
				mniTable[0] = atoi(pch);
				for(int i=1;i<candidate->getGraph()->getNumOfNodes();i++)
				{
					pch = strtok (NULL, ",");
					mniTable[i] = atoi(pch);
				}

				int freq = candidate->subTaskDone(mniTable, 0, this->support);

				delete[] mniTable;

				cout<<"Recieved result for subtask for candidate["<<candidate->getID()<<", freq = "<<freq<<endl;
				printf("Master received subtask from %d\n", status.MPI_SOURCE);

				//>-1 means that all subtasks for this pattern are done, you can rely on the resulted frequency
				if(freq>-1)
				{
					currentlyChecking.erase(candidID);
					if(Settings::debugMSG)
						cout<<"Candidate: "<<candidID<<" is erased from currentlyChecking."<<endl;

					if(freq>=support)
					{
						this->patternFoundFrequent(candidate);
						if(Settings::debugMSG)
							cout<<"Pattern added to frequents"<<endl;
					}
					else
					{
						patternFoundInFrequent(candidate);
						if(Settings::debugMSG)
							cout<<"Pattern added to infrequents"<<endl;
					}
				}

				availableWorker[status.MPI_SOURCE] = true;

				masterProcessingTime += (getmsofday() - start1);
			}
		}

		for(int i=1;i<numWorkers;i++)
		{
			if(!availableWorker[i])
			{
				cout<<"worker: "<<i<<" is not available. numWorkers = "<<numWorkers<<endl;
				continue;
			}
			else
			{
				cout<<"worker: "<<i<<" is expecting a task."<<endl;
			}
			//get a free candidate
			Pattern* candidate = 0;
			string patternKey;
			//get a candidate from the list of exact candidates
			for(int i=0;i<candidates.size();i++)
			{
				CLMap* tempListCandids = candidates.at(i);
				if(tempListCandids->getSize()>0)
				{
					candidate = tempListCandids->getFirstElement().pattern;
					patternKey = tempListCandids->getFirstElement().key;

					if(processed.getPattern(candidate)!=0)
					{
						candidate = 0;
						tempListCandids->removePattern(candidate);
						continue;
					}

					//remove it from expected_frequents
					if(expected_frequentPatterns)
					{
						for(int i=0;i<expected_frequentPatterns->size();i++)
						{
							Pattern* tempPattern = expected_frequentPatterns->at(i)->getPattern(candidate);
							if(tempPattern!=0)
							{
								candidate->borrowTimeInfor(tempPattern, this->numWorkers);
								expected_frequentPatterns->at(i)->removePattern(tempPattern);
							}
						}
					}
					//remove it from expected_infrequents
					if(expected_infrequentPatterns)
					{
						for(int i=0;i<expected_infrequentPatterns->size();i++)
						{
							Pattern* tempPattern = expected_infrequentPatterns->at(i)->getPattern(candidate);
							if(tempPattern!=0)
							{
								candidate->borrowTimeInfor(tempPattern, this->numWorkers);
								expected_infrequentPatterns->at(i)->removePattern(tempPattern);
							}
						}
					}

					candidate->setSubtaskTaken();

					cout<<"candidate->subtasking: "<<candidate->getSubtaskingValue()<<"/"<<candidate->getSubtaskingValueFixed()<<" for candidate ID: "<<candidate->getID()<<endl;

					//remove it from the list of exact candidates (if it has no more subtasks)
					if(candidate->subtasksFinished())
						tempListCandids->removePattern(candidate);

					break;
				}
			}

			//if we could not find an exact candidate, go to the approximated ones
			if(!candidate)
			{
				//try to get one from the expected frequent subgraphs
				vector<CLMap*>* expectedLists[2];
				expectedLists[0] = expected_frequentPatterns;
				expectedLists[1] = expected_infrequentPatterns;

				for(int i=0;i<2;i++)
				{
					if(expectedLists[i] && expectedLists[i]->size()>0)
					{
						while(expectedLists[i]->size())
						{
							CLMap* expCandidates = *(expectedLists[i]->begin());
							if(expCandidates->getSize()==0)
							{
								expectedLists[i]->erase(expectedLists[i]->begin());
								continue;
							}
							//get the first element from the current list
							CLMap_Iterator tempIter = expCandidates->getFirstElement();
							candidate = tempIter.pattern;//expCandidates->begin()->second;
							patternKey = tempIter.key;//expCandidates->begin()->first;

							candidate->setSubtaskTaken();
							cout<<"candidate->subtasking: "<<candidate->getSubtaskingValue()<<endl;

							if(candidate->subtasksFinished())
								expCandidates->removePattern(candidate);

							if(processed.getPattern(candidate)!=0)
							{
								expCandidates->removePattern(candidate);
								candidate = 0;
								continue;
							}
							else
							{
								break;
							}
						}
					}

					if(candidate)
					{
						candidate->makeIDNegative();
						break;
					}
				}
			}

			//send the selected candidate to a worker
			if(candidate)
			{
				sendACandidate(patternKey, candidate, currentlyChecking, i);
			}
			else
			{
				cout<<"No task for this worker "<<i<<" is available"<<endl;
				break;
			}
		}

		cout<<"currentlyChecking.size() = "<<currentlyChecking.size()<<", getNumElems(candidates) = "<<getNumElems((vector<map<string, void*>* >*)&candidates)<<endl;
		cout<<"Currently checking: ";
		for(map<int, Pattern* >::iterator iter = currentlyChecking.begin(); iter!=currentlyChecking.end();iter++)
		{
			Pattern* tempPattern = iter->second;
			cout<<tempPattern->getID()<<", ";
		}
		cout<<endl;

		//we are done, we cannot get more candidates
		if(currentlyChecking.size()==0 && getNumElems((vector<map<string, void*>* >*)&candidates)==0)
		{
			masterProcessingTime += (getmsofday() - start1);
			break;
		}

		masterProcessingTime += (getmsofday() - start1);
	}
	//End of new loop

	printResult();

	cout<<"Number of calls to IsFreq = "<<MinerX::numIsFreqCalls<<endl;
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

void MinerX::loadGraph(string baseName, int graphType, int support, CLMap& freqEdges)
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
	if(!graph->loadFromFile(partFileName, allMPs))
	{
		delete graph;
	}
	nodesCounter+=graph->getNumOfNodes();

	for(tr1::unordered_map<string, void* >::iterator iter = allMPs.begin(); iter!=allMPs.end();++iter)
	{
		if(((Pattern*)((*iter).second))->getFrequency()>=support)
			freqEdges.addPattern((Pattern*)((*iter).second));
		else
			delete (Pattern*)((*iter).second);
	}

	cout<<"Master: graph loaded successfully!"<<endl;
}

void MinerX::sendACandidate(string key, Pattern* candidate, map<int, Pattern*>& currentlyChecking, int destination)
{
	if(Settings::debugMSG)
		cout<<"Worker "<<destination<<" will be busy!"<<endl;
	long long start = getmsofday();
	if(Settings::debugMSG)
	{
		cout<<"Check frequency for:"<<endl;
		cout<<"CL:"<<(*(candidate->getGraph())).getCanonicalLabel()<<endl;
		cout<<*(candidate->getGraph())<<endl<<flush;
	}

	//send isFreq request for the current candidate to worker
	currentlyChecking.insert(std::pair<int, Pattern*>(candidate->getID(), candidate));

	//the search partition to process, -1 indicates no partitioning for this pattern
	int partToProcess;
	int numParts;
	if(candidate->getSubtaskingValueFixed()>1)
	{
		partToProcess = candidate->getSubtaskingValue();
		numParts = candidate->getSubtaskingValueFixed();
	}
	else
	{
		partToProcess = -1;
		numParts = -1;
	}

	ostringstream tempOS;
	tempOS<<*(candidate->getGraph());
	char graphStr[2000];
	sprintf(graphStr, "s,%d,-1,%d,%d,%s\t", candidate->getID(), partToProcess, numParts, tempOS.str().c_str());
	int cnt=strlen(graphStr)+1;
	MPI_Status status;
	MPI_Send(graphStr, cnt, MPI_CHAR, destination, 0, MPI_COMM_WORLD);

	if(Settings::debugMSG)
	{
		printf("Sending subgraph[%d] to worker: %d for frequency check\n", candidate->getID(),destination);
		printf("Message is: [%s]\n", graphStr);
	}

	//add it to processed only when all subtasks are done
	if(candidate->getSubtaskingValue()==0)
		processed.addPattern(candidate);

	availableWorker[destination] = false;

	MinerX::numIsFreqCalls++;
}

void MinerX::extendFreqEdges()
{
	while(candidates.size()<=2)
	{
		candidates.push_back(new CLMap());
	}
	CLMap* twoEdgesCandidate = candidates.at(2);

	CLMap_Iterator iter1 = frequentEdges.getFirstElement();
	while(iter1.pattern!=0)
	{
		Pattern* edge1 = iter1.pattern;

		double l1_0 = edge1->getGraph()->getNodeWithID(0)->getLabel();
		double l1_1 = edge1->getGraph()->getNodeWithID(1)->getLabel();

		CLMap_Iterator iter2 = iter1.getCopy();

		while(iter2.pattern!=0)
		{
			Pattern* edge2 = iter2.pattern;
			frequentEdges.advanceIterator(iter2);

			double l2_0 = edge2->getGraph()->getNodeWithID(0)->getLabel();
			double l2_1 = edge2->getGraph()->getNodeWithID(1)->getLabel();

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

				Pattern* candidate = new Pattern(edge1);
				candidate->invalidateFrequency();

				switch(connectType)
				{
				case 0://00
					candidate->extend(0, 2, edge2->getGraph()->getNodeWithID(1)->getLabel(), edge2->getGraph()->getEdgeLabel(0,1));
					break;
				case 1://01
					candidate->extend(0, 2, edge2->getGraph()->getNodeWithID(0)->getLabel(), edge2->getGraph()->getEdgeLabel(0,1));
					break;
				case 2://10
					candidate->extend(1, 2, edge2->getGraph()->getNodeWithID(1)->getLabel(), edge2->getGraph()->getEdgeLabel(0,1));
					break;
				case 3://11
					candidate->extend(1, 2, edge2->getGraph()->getNodeWithID(0)->getLabel(), edge2->getGraph()->getEdgeLabel(0,1));
					break;
				}

				if(twoEdgesCandidate->getPattern(candidate)!=0)
				{
					delete candidate;
				}
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

void MinerX::setFrequentEdges(CLMap& freqEdges)
{
	this->frequentEdges.addAll(&freqEdges);
}

void MinerX::setSupport(int numWorkers, int support)
{
	this->support = support;

	if(Settings::debugMSG)
		cout<<"Sending new support value to workers ..."<<endl;
	for(int i=1;i<numWorkers;i++)
	{
		int frequency = support;
		char supportStr[500];
		sprintf(supportStr, "t%d\n", frequency);
		int cnt=strlen(supportStr)+1;
		MPI_Send(supportStr, cnt, MPI_CHAR, i, 0, MPI_COMM_WORLD);

		if(Settings::debugMSG)
			cout<<"New support sent to worker: "<<i<<endl;
	}
		cout<<"Finished sending to all workers."<<endl;
}

void MinerX::setExpectedPatterns(vector<CLMap* >* exFreqPatt, vector<CLMap* >* exInfreqPatt)
{
	expected_frequentPatterns = exFreqPatt;
	expected_frequentPatterns->at(1)->clear();
	expected_infrequentPatterns = exInfreqPatt;

	//get the median time, think about other options to select the time
	vector<CLMap* >* allPatterns[2];// = new vector*[2];
	allPatterns[0] = expected_frequentPatterns;
	allPatterns[1] = expected_infrequentPatterns;

	if(Settings::divideBigTasks)
	{
		unsigned long int value;

		//to use median
		/*vector<int> times;
		for(int l=0;l<2;l++)
		{
			for(vector<map<string, Pattern*>* >::iterator iter = allPatterns[l]->begin();iter!=allPatterns[l]->end();iter++)
			{
				map<string, Pattern*>* tempMap = (*iter);
				for(map<string, Pattern*>::iterator iter1 = tempMap->begin(); iter1!=tempMap->end();iter1++)
				{
					Pattern* pattern = iter1->second;
					unsigned long int tempTime = pattern->getPredictedTime();

					//put in order
					vector<int>::iterator iter2 = times.begin();
					for(;iter2!=times.end();iter2++)
					{
						if(tempTime<(*iter2))
							break;
					}
					times.insert(iter2, tempTime);
				}
			}
		}

		value = times.at(times.size()/2);
		if(value==0) value = 1;
		*/

		//use average
		value = 0;
		int count = 0;
		for(int l=0;l<2;l++)
		{
			for(vector<CLMap* >::iterator iter = allPatterns[l]->begin();iter!=allPatterns[l]->end();iter++)
			{
				CLMap* tempMap = (*iter);
				CLMap_Iterator iter1 = tempMap->getFirstElement();
				while(iter1.pattern!=0)
				//for(map<string, Pattern*>::iterator iter1 = tempMap->begin(); iter1!=tempMap->end();iter1++)
				{
					Pattern* pattern = iter1.pattern;//->second;
					tempMap->advanceIterator(iter1);
					unsigned long int tempTime = pattern->getPredictedTime();

					value+=tempTime;
					count++;
				}
			}
		}

		value = ((double)value)/count;

		if(Settings::debugMSG)
			cout<<"The selected average value is: "<<value<<endl;

		//use the median to assign how many tasks per pattern is needed
		for(int l=0;l<2;l++)
		{
			for(vector<CLMap* >::iterator iter = allPatterns[l]->begin();iter!=allPatterns[l]->end();iter++)
			{
				CLMap* tempMap = (*iter);
				CLMap_Iterator iter1 = tempMap->getFirstElement();
				while(iter1.pattern!=0)
				//for(map<string, Pattern*>::iterator iter1 = tempMap->begin(); iter1!=tempMap->end();iter1++)
				{
					Pattern* pattern = iter1.pattern;//->second;
					tempMap->advanceIterator(iter1);
					unsigned long int tempTime = pattern->getPredictedTime();

					if(tempTime>value*2)
					{
						int numTasks = ceil(tempTime/value);
						if(numTasks>(numWorkers-1))
							numTasks = (numWorkers-1);
						if(Settings::debugMSG)
							cout<<"numTasks:"<<numTasks<<endl;
						pattern->setSubtasking(numTasks, this->numWorkers);
					}
				}
			}
		}
	}
}

void MinerX::printCandidates()
{
	print(candidates);
}

void MinerX::printResult()
{
	print(frequentPatterns);
}

void MinerX::removePattern(Pattern* pattern, vector<CLMap* >& data)
{
	CLMap* tempList = data.at(pattern->getSize());
	tempList->removePattern(pattern);
}

int MinerX::getNumOfExactPatterns(map<int, Pattern*>& mapP)
{
	int count = 0;
	for(map<int, Pattern*>::iterator iter = mapP.begin(); iter!=mapP.end();iter++)
	{
		if(iter->first>=0)
			count++;
	}
	return count;
}

void MinerX::print(vector<CLMap* >& data)
{
	//count the total number of frequent patterns
	int size = 0;
	for(unsigned int i=0;i<data.size();i++)
	{
		size+=data[i]->getSize();
	}
	int count = 1;
	cout<<"Miner found "<<size<<" frequent patterns, and they are:"<<endl;
	for(unsigned int i=0;i<data.size();i++)
	{
		cout<<"With "<<(i)<<" edges:"<<endl;
		data[i]->print(count);
		count+=data[i]->getSize();
	}
}

void MinerX::patternFoundFrequent(Pattern* candidate)
{
	try
	{

	long long start1 = getmsofday();

	if(Settings::debugMSG)
		cout<<"Found to be FREQUENT, with frequency = "<<candidate->getFrequency()<<endl;

	//add to the frequent patterns list
	while(frequentPatterns.size()<=candidate->getSize())
	{
		frequentPatterns.push_back(new CLMap());
	}

	frequentPatterns[candidate->getSize()]->addPattern(candidate);

	//extend using the current frequent patterns with the same size
	CLMap_Iterator iter1 = frequentPatterns[candidate->getSize()]->getFirstElement();
	while(iter1.pattern!=0)
	{
		Pattern* temp = iter1.pattern;//->second;
		frequentPatterns[candidate->getSize()]->advanceIterator(iter1);

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

				if(Settings::debugMSG)
					cout<<"NEW CANDID generated"<<endl;

				if(checkAgainstOtherPatterns(newCandidate))
				{
					delete newGraph;
					delete newCandidate;
				}
				else
				{
					while(candidates.size()<=newCandidate->getSize())
					{
						candidates.push_back(new CLMap());
					}

					newCandidate->generatePrimaryGraphs();
					if(candidates.at(newCandidate->getSize())->getPattern(newCandidate)==0)
					{
						if(Settings::debugMSG)
							cout<<"Pattern: "<<newCandidate->getID()<<" is added to candidates."<<endl;
						bool b = candidates.at(newCandidate->getSize())->addPattern(newCandidate);
						if(!b)
						{
							delete newGraph;
							delete newCandidate;
						}
					}
				}

				//if the two edges have the same label, we need to assume the same label belong to one node
				if(candPG->getEdge()->getOtherNode()->getLabel()==otherPG->getEdge()->getOtherNode()->getLabel() &&
						candPG->getSrcNodeID()!=scrNodeID)
				{
					newGraph = new GraphX(candidate->getGraph());
					newGraph->addEdge(candPG->getEdge()->getOtherNode()->getID(), scrNodeID, otherPG->getEdge()->getLabel());
					newCandidate = new Pattern(newGraph, false);

					if(checkAgainstOtherPatterns(newCandidate))
					{
						delete *iter2;
						delete newGraph;
						delete newCandidate;
						continue;
					}

					newCandidate->generatePrimaryGraphs();

					while(candidates.size()<=newCandidate->getSize())
					{
						candidates.push_back(new CLMap());
					}

					//CHECK THIS LINE IS TAKING SO MUCH TIME WHEN THE PATTERN IS BIG!!!! ---->
					if(candidates.at(newCandidate->getSize())->getPattern(newCandidate)==0)
					{
						bool b = candidates.at(newCandidate->getSize())->addPattern(newCandidate);
						if(!b)
						{
							delete newGraph;
							delete newCandidate;
						}

						if(Settings::debugMSG)
							cout<<"Pattern: "<<newCandidate->getID()<<" is added to candidates."<<endl;
					}
				}
				delete *iter2;
			}
		}
	}

	MinerX2::patternFreqTime+=(getmsofday()-start1);

	}
	catch(std::bad_alloc& exc)
	{
	  cout<<"Error: Bad allocation found in (patternFoundToBeFrequent)!"<<endl<<flush;
	  exit(0);
	}
}

/**
 * check th enew candidate against other lists
 * return true to avoid adding it
 */
bool MinerX::checkAgainstOtherPatterns(Pattern* newCandidate)
{
	//do not add it if it is already done
	if(processed.getPattern(newCandidate)==0)
		;
	else
	{
		if(Settings::debugMSG)
			cout<<newCandidate->getGraph()->getCanonicalLabel()<<" not added (1)!!!!!!!!!"<<endl;
		return true;
	}

	//do not add it if it is already currently checking
	//think about the case where it is not part of the currentlyChecking list, but it still has remaining subtasks
	bool foundInCurrChecking = false;
	for(map<int, Pattern* >::iterator iter = currentlyChecking.begin(); iter!=currentlyChecking.end();iter++)
	{
		Pattern* tempPattern = iter->second;
		if(tempPattern->getGraph()->isTheSameWith(newCandidate->getGraph()))
		{
			foundInCurrChecking = true;
			break;
		}
	}
	if(foundInCurrChecking)
		return true;

	//erase the same candidate from expected frequents
	if(expected_frequentPatterns)
	{
		for(int i=0;i<expected_frequentPatterns->size();i++)
		{
			Pattern* pattern = expected_frequentPatterns->at(i)->getPattern(newCandidate);
			if(pattern!=0)
			{
				newCandidate->borrowTimeInfor(pattern, this->numWorkers);
				expected_frequentPatterns->at(i)->removePattern(newCandidate);
				break;
			}
		}
	}
	//erase the same candidate from expected infrequents
	if(expected_infrequentPatterns)
	{
		for(int i=0;i<expected_infrequentPatterns->size();i++)
		{
			Pattern* pattern = expected_infrequentPatterns->at(i)->getPattern(newCandidate);
			if(pattern!=0)
			{
				newCandidate->borrowTimeInfor(pattern, this->numWorkers);
				expected_infrequentPatterns->at(i)->removePattern(newCandidate);
				break;
			}
		}
	}

	return false;
}

void MinerX::patternFoundInFrequent(Pattern* candidate)
{
	while(infrequentPatterns.size()<=candidate->getSize())
	{
		infrequentPatterns.push_back(new CLMap());
	}

	infrequentPatterns[candidate->getSize()]->addPattern(candidate);

	//remove it and its supergraphs from the set of expected_infrequentPatterns
	for(int i=candidate->getSize();i<expected_infrequentPatterns->size();i++)
	{
		CLMap* currentList = expected_infrequentPatterns->at(i);
		CLMap_Iterator iter = currentList->getFirstElement();
		while(iter.pattern!=0)
		{
			Pattern* currentcandid = iter.pattern;//->second;

			if(currentcandid->getGraph()->getNumOfNodes()<candidate->getGraph()->getNumOfNodes())
			{
				currentList->advanceIterator(iter);
				continue;
			}

			vector<map<int, int>* > result;
			tr1::unordered_map<int, tr1::unordered_set<int>*> domains_values;

			if(Settings::debugMSG)
			{
				cout<<"candidate:"<<endl;
				cout<<(*(candidate->getGraph()))<<endl<<flush;

				cout<<"currentcandid:"<<endl;
				cout<<(*(currentcandid->getGraph()))<<endl<<flush;
			}

			candidate->getGraph()->isIsomorphic(currentcandid->getGraph(), result, domains_values);

			if(result.size()>0)
			{
				Pattern* tempPattern = iter.pattern;
				currentList->advanceIterator(iter);
				currentList->removePattern(tempPattern);
				this->currentlyChecking.erase(currentcandid->getID());

				for(vector<map<int, int>* >::iterator iter1 = result.begin();iter1!=result.end();iter1++)
				{
					delete (*iter1);
				}
				result.clear();
			}
			else
			{
				currentList->advanceIterator(iter);
			}
		}
	}
}

bool MinerX::foundInCurrentlyChecking(Pattern* newCandidate)
{
	for(map<int, Pattern* >::iterator iter = currentlyChecking.begin(); iter!=currentlyChecking.end();iter++)
	{
		Pattern* tempPattern = iter->second;
		if(tempPattern->getGraph()->isTheSameWith(newCandidate->getGraph()))
		{
			return true;
		}
	}

	return false;
}

MinerX::~MinerX()
{
	/*for(map<string, Pattern*>::iterator iter = frequentEdges.begin();iter!=frequentEdges.end();iter++)
	{
		delete iter->second;
	}*/

	cout<<"DM1_"<<flush;
	vect_map_destruct(frequentPatterns);
	cout<<"DM2_"<<flush;
	vect_map_destruct(infrequentPatterns);
	cout<<"DM3_"<<flush;
	vect_map_destruct(candidates);
	cout<<"DM4_"<<flush;
}
