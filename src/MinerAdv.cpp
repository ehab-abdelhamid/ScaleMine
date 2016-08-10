/**
 * Copyright 2016 Ehab Abdelhamid, Ibrahim Abdelaziz

MinerAdv.cpp

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
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include "MinerAdv.h"
#include "utils.h"
#include "Settings.h"

/**
 * start the mining process, given the filename (file base name for the partitions)
 * , and the required support threshold
 */
void MinerAdv::startMining(string fileName, int graphType, int support, int numMachine, int nThreads)
{
	this->nThreads = nThreads;
	int numWorkers = (numMachine)*nThreads;
	//create available workers list
	bool* availableWorker = new bool[numWorkers];
	for(int i=1;i<numWorkers;i++)
	{
		availableWorker[i] = true;
	}

	int numberOfCreatedCandids = 0;

	numOfVisitedNodes = 0;
	numIterations = 0;

	//first time: send task to all workers

	cout<<"MinerAdv: Sending graph filename to all workers: "<<endl;

	for(int i=1;i<numMachine;i++)
	{
		//Graph Loading Task
		//produceTask(rank, size, str_message);
		int frequency = support;
		char fileName_[500];
		sprintf(fileName_, "f%s,%d\n", fileName.c_str(), frequency);
		int cnt=strlen(fileName_)+1;
		MPI_Send(fileName_, cnt, MPI_CHAR, i, 0, MPI_COMM_WORLD);

		if(Settings::debugMSG)
			cout<<"Sent to worker: "<<i<<endl;

		for(int j=0;j<nThreads;j++)
			availableWorker[i*nThreads+j] = false;
	}
	if(Settings::debugMSG)
		cout<<"Finished sending to all workers."<<endl;

	this->support = support;

	//load the graph
	long long start = getmsofday();
	loadGraph(fileName, graphType, support, frequentEdges);

	//waiting for all workers to load the graph
	int tempNumWorkers = numMachine-1;
	while(tempNumWorkers>0)
	{
		MPI_Status status;
		char str_message[2000];
		MPI_Recv(str_message, 2000, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		int srcThread = status.MPI_SOURCE*nThreads+status.MPI_TAG;
		availableWorker[srcThread] = true;

		if(Settings::debugMSG)
		{
			cout<<"Master recieved masg: "<<str_message<<endl;
		}

		if(str_message[0]=='l')
		{
			if(Settings::debugMSG)
				printf("Master received ack from %d(%d:%d)\n", status.MPI_SOURCE*nThreads+status.MPI_TAG, status.MPI_SOURCE, status.MPI_TAG);

			for(int j=0;j<nThreads;j++)
				availableWorker[status.MPI_SOURCE*nThreads+j] = true;
		}

		tempNumWorkers--;
		cout<<"Remaining workers to complete loading: "<<tempNumWorkers<<endl;
	}

	if(frequentEdges.getSize()==0)
	{
		cout<<"No frequent patterns found! Exiting"<<endl;
		exit(0);
	}

	long long end = getmsofday();
	long long elapsed = end - start;

	cout<<"Loading took "<<(elapsed/1000)<<" sec and "<<(elapsed%1000)<<" ms"<<endl;


	Settings::graphLoadingTime = elapsed;

	CLMap* clmap = new CLMap();
	clmap->addAll(&frequentEdges);

	frequentPatterns.insert(frequentPatterns.begin()+1, clmap);//frequent edges!

	if(Settings::debugMSG)
	{
		cout<<"#frequent edges  = "<<frequentEdges.getSize()<<endl;
		cout<<"Frequent edges:"<<endl;
		frequentEdges.print();
	}

	//extend frequent edges into 2 edges graphs

	start = getmsofday();

	extendFreqEdges();

	end = getmsofday();

	cout<<"Start mining approximate ..."<<endl;

	bool firstTime = true;

	while(true)
	{
		//check if we have more candidates with the same current size
		if(currentlyChecking.size()==0 && getNumElems((vector<map<string, void*>* >*)&candidates)==0)
		{
			break;
		}

		//send available candidates to available workers
		for(int i=nThreads;i<numWorkers;i++)
		{
			if(availableWorker[i])
			{
				if(popACandidateApprox(candidates, currentlyChecking, i))
				{
					availableWorker[i] = false;
				}
				else
				{
					break;
				}
			}
		}

		MPI_Status status;
		char str_message[2000];

		MPI_Recv(str_message, 2000, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		int srcThread = status.MPI_SOURCE*nThreads+status.MPI_TAG;
		availableWorker[srcThread] = true;
		firstTime = false;

		if(Settings::debugMSG)
		{
			cout<<"Master recieved masg: "<<str_message<<endl;
		}

		//ack message recieved
		if(str_message[0]=='l')
		{
			if(Settings::debugMSG)
				printf("Master received ack from %d(%d:%d)\n", status.MPI_SOURCE*nThreads+status.MPI_TAG, status.MPI_SOURCE, status.MPI_TAG);

			for(int j=0;j<nThreads;j++)
				availableWorker[status.MPI_SOURCE*nThreads+j] = true;
		}
		else if(str_message[0]=='t')
		{
			//get the pattern that was sent
			int candidID;
			unsigned long time;
			sscanf(str_message, "t%d,%Lu,%Lu", &candidID, &time);
			Pattern* candidate = currentlyChecking.find(candidID)->second;
			candidate->setFrequency(support);
			currentlyChecking.erase(candidID);
			if(Settings::debugMSG)
				cout<<"Removing candidate#"<<candidID<<" from currentlyChecking."<<endl;

			cout<<"Predicted Time = "<<candidate->getPredictedTime()<<endl;

			char temp[500];
			sprintf(temp, "%d_%d", srcThread, candidID);
			isFrerqtimes.insert(std::pair<string, unsigned long>(temp, time));

			if(Settings::debugMSG)
			{
				printf("Master received isFreq[TRUE] from %d\n", srcThread);
				cout<<"Found to be FREQUENT, with frequency = "<<candidate->getFrequency()<<endl;
			}
			//add to the frequent patterns list
			if((frequentPatterns.size()-1)<candidate->getSize())
			{
				frequentPatterns.insert(frequentPatterns.begin()+candidate->getSize(), new CLMap());
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
						if(candidates.size()<=newCandidate->getSize())
						{
							candidates.insert(candidates.begin()+newCandidate->getSize(), new CLMap());
						}
						newCandidate->generatePrimaryGraphs();
						bool b = candidates.at(newCandidate->getSize())->addPattern(newCandidate);
						if(!b) delete newCandidate;
						//if the two edges have the same label, we need to assume the same label belong to one node
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
		}
		else if(str_message[0]=='f')
		{
			//get the pattern that was sent
			int candidID;
			unsigned long time;
			sscanf(str_message, "f%d,%Lu,%Lu", &candidID, &time);
			Pattern* candidate = currentlyChecking.find(candidID)->second;

			if(Settings::debugMSG)
				cout<<"Predicted Time = "<<candidate->getPredictedTime()<<endl;
			candidate->setFrequency(0);
			currentlyChecking.erase(candidID);

			if(Settings::debugMSG)
				cout<<"Removing candidate#"<<candidID<<" from currentlyChecking."<<endl;

			char temp[500];
			sprintf(temp, "%d_%d", srcThread, candidID);
			isFrerqtimes.insert(std::pair<string, unsigned long>(temp, time));

			while(infrequentPatterns.size()<=candidate->getSize())
			{
				infrequentPatterns.push_back(new CLMap());
			}
			infrequentPatterns[candidate->getSize()]->addPattern(candidate);
		}
		//we got results from approximate freq computation
		else if(str_message[0]=='a')
		{
			if(str_message[1]=='t')
			{
				//get the pattern that was sent
				int candidID;
				int frequency;
				unsigned long time;
				unsigned long predictedTime;
				int exact;
				unsigned long numVisitedNodes;
				unsigned long numIters;
				sscanf(str_message, "at%d,%d,%Lu,%Lu,%d,%Lu,%Lu", &candidID, &frequency, &time, &predictedTime, &exact, &numVisitedNodes, &numIters);

				numOfVisitedNodes += numVisitedNodes;
				numIterations += numIters;

				map<int, Pattern* >::iterator f_iter = currentlyChecking.find(candidID);
				if(f_iter==currentlyChecking.end())
				{
					cout<<"[T] Pattern#"<<candidID<<" is not found in currentlyChecking"<<endl;
					cout<<flush;
					exit(1);
				}
				Pattern* candidate = f_iter->second;
				candidate->setFrequency(frequency);
				candidate->setPredictedTime(predictedTime);
				candidate->setPredictedFreq(frequency);
				candidate->setMaxIters(numIters);
				if(exact==1)
				{
					if(Settings::debugMSG)
						cout<<"Set result exact for this pattern, freq = "<<frequency<<endl;
					candidate->setResultExact();
				}

				if(Settings::debugMSG)
					cout<<"Predicted Time: "<<predictedTime<<" for candidateID: "<<candidID<<endl;
				currentlyChecking.erase(candidID);

				if(Settings::debugMSG)
					cout<<"Removing candidate#"<<candidID<<" from currentlyChecking."<<endl;

				char temp[500];
				sprintf(temp, "%d_%d", srcThread, candidID);
				isFrerqtimes.insert(std::pair<string, unsigned long>(temp, time));

				if(Settings::debugMSG)
				{
					printf("Master received isFreqApprox[TRUE] from %d.\n", srcThread);
					cout<<"Found to be FREQUENT, with frequency = "<<candidate->getFrequency()<<endl;
				}
				//add to the frequent patterns list
				if((frequentPatterns.size()-1)<candidate->getSize())
				{
					frequentPatterns.insert(frequentPatterns.begin()+candidate->getSize(), new CLMap());
				}
				frequentPatterns[candidate->getSize()]->addPattern(candidate);

				if(MinerAdv::maxNumCandids==-1 || numberOfCreatedCandids<MinerAdv::maxNumCandids)
				{
					if(Settings::debugMSG)
					{
						cout<<"Extending the frequent pattern ..."<<endl;
					}

					//extend using the current frequent patterns with the same size
					CLMap_Iterator iter1 = frequentPatterns[candidate->getSize()]->getFirstElement();
					while(iter1.pattern!=0)
					//for(map<string, Pattern*>::iterator iter1 = frequentPatterns[candidate->getSize()]->begin();iter1!=frequentPatterns[candidate->getSize()]->end();++iter1)
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
								if(candidates.size()<=newCandidate->getSize())
								{
									candidates.insert(candidates.begin()+newCandidate->getSize(), new CLMap());
								}
								newCandidate->generatePrimaryGraphs();
								if(candidates.at(newCandidate->getSize())->getPattern(newCandidate)==0)
								{
									bool b = candidates.at(newCandidate->getSize())->addPattern(newCandidate);
									if(!b) delete newCandidate;
									else
										numberOfCreatedCandids++;
								}
								else
								{
									delete newCandidate;
								}
								//if the two edges have the same label, we need to assume the same label belong to one node
								if(candPG->getEdge()->getOtherNode()->getLabel()==otherPG->getEdge()->getOtherNode()->getLabel() &&
										candPG->getSrcNodeID()!=scrNodeID)
								{
									newGraph = new GraphX(candidate->getGraph());
									newGraph->addEdge(candPG->getEdge()->getOtherNode()->getID(), scrNodeID, otherPG->getEdge()->getLabel());
									newCandidate = new Pattern(newGraph, false);
									newCandidate->generatePrimaryGraphs();
									//CHECK THIS LINE IS TAKING SO MUCH TIME WHEN THE PATTERN IS BIG!!!! ---->
									if(Settings::debugMSG)
									{
										cout<<"candidates.size = "<<candidates.size()<<", newCandidate->getSize() = "<<newCandidate->getSize()<<endl;
									}

									if(candidates.at(newCandidate->getSize())->getPattern(newCandidate)==0)
									{
										bool b = candidates.at(newCandidate->getSize())->addPattern(newCandidate);
										if(!b) delete newCandidate;
										else
											numberOfCreatedCandids++;
									}
									else
									{
										delete newCandidate;
									}
								}
								delete *iter2;
							}
						}
					}
					if(Settings::debugMSG)
					{
						cout<<"Finished Extending the frequent pattern."<<endl<<flush;
					}
				}
				else
				{
					if(Settings::debugMSG)
					{
						cout<<"Not extending the patterns as we already generated many candidates:"<<numberOfCreatedCandids<<endl;
					}
				}
			}
			else if(str_message[1]=='f')
			{
				//get the pattern that was sent
				int candidID;
				int frequency;
				unsigned long time;
				unsigned long predictedTime;
				int invalidCol;
				int predictedValids;
				int exact;
				unsigned long numVisitedNodes;
				unsigned long numIters;
				sscanf(str_message, "af%d,%d,%Lu,%Lu,%d,%d,%d,%Lu,%Lu", &candidID, &frequency, &time, &predictedTime, &invalidCol, &predictedValids, &exact, &numVisitedNodes, &numIters);

				numOfVisitedNodes += numVisitedNodes;
				numIterations += numIters;

				map<int, Pattern* >::iterator f_iter = currentlyChecking.find(candidID);
				if(f_iter==currentlyChecking.end())
				{
					cout<<"[F] Pattern#"<<candidID<<" is not found in currentlyChecking"<<endl;
					cout<<flush;
					exit(1);
				}

				Pattern* candidate = f_iter->second;
				candidate->setFrequency(frequency);
				candidate->setPredictedTime(predictedTime);
				candidate->setInvalidCol(invalidCol, predictedValids);
				candidate->setMaxIters(numIters);
				if(exact==1)
				{
					if(Settings::debugMSG)
						cout<<"Set result exact for this pattern"<<endl;
					candidate->setResultExact();
				}
				currentlyChecking.erase(candidID);
				if(Settings::debugMSG)
					cout<<"Removing candidate#"<<candidID<<" from currentlyChecking."<<endl;

				if(Settings::debugMSG)
					cout<<"Predicted Time: "<<predictedTime<<" for candidateID: "<<candidID<<", invalid col: "<<invalidCol<<" with #valid nodes = "<<predictedValids<<endl;

				char temp[500];
				sprintf(temp, "%d_%d", srcThread, candidID);
				isFrerqtimes.insert(std::pair<string, unsigned long>(temp, time));

				while(infrequentPatterns.size()<=candidate->getSize())
				{
					infrequentPatterns.push_back(new CLMap());
				}
				infrequentPatterns[candidate->getSize()]->addPattern(candidate);
			}
		}

		if(Settings::debugMSG)
		{
			cout<<"currentlyChecking.size() = "<<currentlyChecking.size()<<", getNumElems(candidates) = "<<getNumElems((vector<map<string, void*>* >*)&candidates)<<endl;
			cout<<"Currently checking: ";
			for(map<int, Pattern* >::iterator iter = currentlyChecking.begin(); iter!=currentlyChecking.end();iter++)
			{
				Pattern* tempPattern = iter->second;
				cout<<tempPattern->getID()<<", ";
			}
			cout<<endl;
		}
		else
		{
//			int r = rand()%100;
//			if(r==0)
//			{
//				cout<<"#CANDIDATES = "<<getNumElems((vector<map<string, void*>* >*)&candidates)<<". #currentlyChecking = "<<currentlyChecking.size()<<endl;
//			}
		}
		//end parallelizing

		if(Settings::debugMSG)
			cout<<"---------------------------------------------"<<endl;
	}

	cout<<"Mining approximate DONE."<<endl;

	delete[] availableWorker;
}

bool MinerAdv::popACandidateApprox(vector<CLMap*>& candidates, map<int, Pattern*>& currentlyChecking, int destination)
{
	CLMap* currentCandidates = 0;
	for(vector<CLMap*>::iterator iter = candidates.begin();iter!=candidates.end();iter++)
	{
		if((*iter)->getSize()!=0)
		{
			currentCandidates = (*iter);
			break;
		}
	}

	if(currentCandidates==0)
	{
		return false;
	}

	if(Settings::debugMSG)
		cout<<"Popping a candidate ..."<<endl;
	//get a candidate, if there is
	if(currentCandidates->getSize()>0)
	{
		CLMap_Iterator iter = currentCandidates->getFirstElement();
		string key = iter.key;
		Pattern* candidate = iter.pattern;
		currentCandidates->removePattern(candidate);

		sendACandidateApprox(key, candidate, currentlyChecking, destination);
	}

	return true;
}

void MinerAdv::sendACandidateApprox(string key, Pattern* candidate, map<int, Pattern*>& currentlyChecking, int destination)
{
	long long start = getmsofday();
	if(Settings::debugMSG)
	{
		cout<<"Check frequency for:"<<endl;
		cout<<"CL:"<<(*(candidate->getGraph())).getCanonicalLabel()<<endl;
		cout<<*(candidate->getGraph())<<endl<<flush;
	}

	//send isFreq request for the current candidate to worker
	if(currentlyChecking.find(candidate->getID())!=currentlyChecking.end())
	{
		cout<<"candidate->getID() = "<<candidate->getID()<<endl<<flush;
		exit(0);
	}

	currentlyChecking.insert(std::pair<int, Pattern*>(candidate->getID(), candidate));

	ostringstream tempOS;
	tempOS<<*(candidate->getGraph());

	char graphStr[2000];

	sprintf(graphStr, "a,%d,%s\t", candidate->getID(),tempOS.str().c_str());
	int cnt=strlen(graphStr)+1;
	MPI_Status status;

	int rank = destination/nThreads;
	int tagID = destination%nThreads;

	//MPI_Send(graphStr, cnt, MPI_CHAR, destination, 0, MPI_COMM_WORLD);
	MPI_Send(graphStr, cnt, MPI_CHAR, rank, tagID, MPI_COMM_WORLD);

	if(Settings::debugMSG)
	{
		printf("Sending subgraph[%d] to worker: %d for approximate frequency check (by sampling)\n", candidate->getID(),destination);
		printf("%s\n", graphStr);
	}

	Miner::numIsFreqCalls++;
}

void MinerAdv::printTotalExpectedTime()
{
	vector<CLMap* >* allPatterns[2];
	allPatterns[0] = &frequentPatterns;
	allPatterns[1] = &infrequentPatterns;

	unsigned long int totatTime = 0;

	//add to candidates
	for(int l=0;l<2;l++)
	{

		for(vector<CLMap* >::iterator iter = allPatterns[l]->begin();iter!=allPatterns[l]->end();iter++)
		{
			CLMap* tempMap = (*iter);
			CLMap_Iterator iter1 = tempMap->getFirstElement();
			while(iter1.pattern!=0)
			{
				Pattern* pattern = iter1.pattern;//->second;
				tempMap->advanceIterator(iter1);
				totatTime+=pattern->getPredictedTime();
			}
		}
	}
}
