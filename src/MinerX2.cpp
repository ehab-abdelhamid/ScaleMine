/**
 * Copyright 2016 Ehab Abdelhamid, Ibrahim Abdelaziz

MinerX2.cpp

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

#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include "utils.h"
#include "MinerX2.h"
#include "Settings.h"
#include "GraMiCounter.h"
#include "QueueSimulator.h"

void MinerX2::setExpectedPatterns(vector<CLMap* >* exFreqPatt, vector<CLMap* >* exInfreqPatt)
{
	expected_frequentPatterns = exFreqPatt;
	expected_infrequentPatterns = exInfreqPatt;

	vector<CLMap* >* allPatterns[2];
	allPatterns[0] = expected_frequentPatterns;
	allPatterns[1] = expected_infrequentPatterns;

	int numPatterns = 0;

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

				while(candidates.size()<=pattern->getSize())
				{
					candidates.push_back(new CLMap());
				}

				CLMap* tempMap = candidates.at(pattern->getSize());
				Pattern* newPattern = new Pattern(pattern);
				newPattern->generatePrimaryGraphs();

				//added for phase 2 to work
				newPattern->borrowTimeInfor(pattern, this->numWorkers);
				newPattern->setPredictedFreq(pattern->getPredictedFreq());

				bool b= tempMap->addPattern(newPattern);
				if(!b) delete newPattern;
				else
				{
					numPatterns++;
				}

			}
		}
	}

	if(Settings::divideBigTasks)
	{
		//use average
		timeThreshold = 0;
		switch(Settings::divideByFunc)
		{
		case 0:
		case 2:
			{
				int count = 0;
				for(int l=0;l<2;l++)
				{
					for(vector<CLMap* >::iterator iter = allPatterns[l]->begin();iter!=allPatterns[l]->end();iter++)
					{
						CLMap* tempMap = (*iter);
						CLMap_Iterator iter1 = tempMap->getFirstElement();
						while(iter1.pattern!=0)
						{
							Pattern* pattern = iter1.pattern;
							tempMap->advanceIterator(iter1);

							unsigned long int tempTime = pattern->getPredictedTime();

							if(tempTime>0)
							{
								timeThreshold+=tempTime;
								count++;
							}
						}
					}
				}

				timeThreshold = ((double)timeThreshold)/count;
				//use the threshold to assign how many tasks per pattern is needed
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
							unsigned long int tempTime = pattern->getPredictedTime();

							if(Settings::debugMSG)
							{
								cout<<"Pattern#"<<pattern->getID()<<endl;
								cout<<"tempTime = "<<tempTime<<", timeThreshold = "<<timeThreshold<<endl;
							}

							if(tempTime>timeThreshold*2)
							{
								int numTasks = ceil(tempTime/timeThreshold);

								//limit the number of partitions for expected frequent subgraphs (to allow good improvement by the optimization of side effects)
								//No need for this in the case of infrequent subgraphs where only one column is required to be processed
								if(pattern->getInvalidCol()==-1)
								{
									if(numTasks>(numWorkers-1))
										numTasks = (numWorkers-1);
								}
								else
								{
									if(numTasks>((numWorkers-1)*2))
										numTasks = ((numWorkers-1)*2);
								}

								if(Settings::debugMSG)
									cout<<"numTasks:"<<numTasks<<endl;
								pattern->setSubtasking(numTasks, this->numWorkers);
							}
						}
					}
				}
			}
			break;
		case 1:
			{
				//use function 1: minimize: threshold+overall_execution_time*#subtasks*overhead_constant
				timeThreshold = 0;
				unsigned long int cumulativeTime = 0;
				int tempCounter = 0;
				//start with an initial value = average

				vector<Pattern* > selectedForOverheadCompute;

				bool onePass = false;

				int counter = 10000;//to prevent forever looping

				while(selectedForOverheadCompute.size()<5)
				{
					counter--;
					if(counter<0)
						break;

					int r = rand()%numPatterns;

					for(vector<CLMap* >::iterator iter = candidates.begin();iter!=candidates.end();iter++)
					{
						CLMap* tempMap = (*iter);
						CLMap_Iterator iter1 = tempMap->getFirstElement();
						while(iter1.pattern!=0)
						{
							Pattern* pattern = iter1.pattern;//->second;
							tempMap->advanceIterator(iter1);
							unsigned long int tempTime = pattern->getPredictedTime();

							if(!onePass)
							{
								if(Settings::debugMSG)
								{
									cout<<"Adding predicted time = "<<tempTime<<endl;
								}
								cumulativeTime+=tempTime;
								tempCounter++;
							}

							r--;
							if(r==0)
							{
								//only select patterns with non-exact predictions, of size more than one edge, and with expected time < 60 seconds
								if(pattern->getGraph()->getNumOfEdges()>1 && !pattern->isResultExact() && pattern->getPredictedTime()<60000)
								{
									selectedForOverheadCompute.push_back(pattern);
								}
							}
						}
					}

					onePass = true;
				}

				if(Settings::debugMSG)
				{
					cout<<"selectedForOverheadCompute.size() = "<<selectedForOverheadCompute.size()<<endl<<flush;
					cout<<"Cumulative expected time is: "<<cumulativeTime<<endl;
				}

				timeThreshold = ((double)cumulativeTime)/tempCounter;//selectedForOverheadCompute.size();

				if(Settings::debugMSG)
				{
					cout<<"Initial value of timeThreshold is: "<<timeThreshold<<endl;
				}

				//computing overhead value
				if(Settings::debugMSG)
					cout<<"Phase II starting: computing overhead value"<<endl;

				unsigned long phast2TimeStart = getmsofday();

				double tempOverhead = 0;
				for(vector<Pattern* >::iterator iter = selectedForOverheadCompute.begin();iter!=selectedForOverheadCompute.end();iter++)
				{
					Pattern* pattern = (*iter);

					//get frequency using a single subtask
					unsigned long startTime = getmsofday();

					if(Settings::debugMSG)
					{
						cout<<"Master starts checking frequency (Single task)"<<endl;
						cout<<"Pattern ID: "<<pattern->getID()<<endl;
						cout<<(*pattern->getGraph())<<endl;
						cout<<"Predicted frequency = "<<pattern->getPredictedFreq()<<". Predicted invalid column: "<<pattern->getInvalidCol()<<". Predicted valids = "<<pattern->getPredictedValids()<<". Predicted time: "<<pattern->getPredictedTime()<<endl;
						cout<<"IS result exact? "<<pattern->isResultExact()<<endl;
					}

					int* mniTable;

					int freq = GraMiCounter::isFrequent_adv(graph, pattern, this->support, -1, -1, mniTable);

					if(Settings::debugMSG)
					{
						cout<<"Master Finished checking frequency (Single task). Frequency is: "<<freq<<endl<<flush;
					}

					pattern->setFrequency(freq);
					pattern->setResultExact();

					unsigned long elapsed_singleTask = getmsofday() - startTime;

					//--------------------------------

					//get frequency using 4 subtasks
					startTime = getmsofday();

					int subTaskNum = 4;
					if(Settings::debugMSG)
						printf("Master Starts checking frequency (multiple subtasks)\n");

					int* mniTable_;

					mniTable_ = new int[pattern->getGraph()->getNumOfNodes()];

					for(int i=0;i<subTaskNum;i++)
					{
						GraMiCounter::isFrequent_adv(graph, pattern, this->support, i, subTaskNum, mniTable_);
					}

					delete[] mniTable_;

					if(Settings::debugMSG)
					{
						printf("Master Finished checking frequency (multi tasks).\n");
					}

					unsigned long elapsed_MultiTask = getmsofday() - startTime;

					if(elapsed_MultiTask>elapsed_singleTask)
					{
						tempOverhead += (elapsed_MultiTask-elapsed_singleTask)/((double)elapsed_singleTask)/(subTaskNum-1);
						cout<<"elapsed_MultiTask = "<<elapsed_MultiTask<<", elapsed_singleTask = "<<elapsed_singleTask<<", diff = "<<(elapsed_MultiTask-elapsed_singleTask)<<", tempOverhead = "<<tempOverhead<<endl;
					}
					else
					{
						cout<<"elapsed_MultiTask<=elapsed_singleTask ["<<elapsed_MultiTask<<", "<<elapsed_singleTask<<"]. tempOverhead = "<<tempOverhead<<endl;
					}
				}

				if(selectedForOverheadCompute.size()>0)
					overhead = tempOverhead/selectedForOverheadCompute.size();
				else
				{
					if(Settings::debugMSG)
						cout<<"Overhead is not modified"<<endl;
				}

				Settings::phase2Time+=(getmsofday()-phast2TimeStart);

				if(Settings::debugMSG)
					cout<<"Phase II Finished: overhead = "<<overhead<<endl;

				if(overhead>0)
				{
					//try lower and upper values, compute the function, select the minimum, until no improvement is found (hill climbing)
					int changeUpper = timeThreshold/2;
					int changeLower = timeThreshold/2;
					while(true)
					{
						long lower = timeThreshold-changeLower;
						long upper = timeThreshold+changeUpper;

						long lowerScore = function1(lower, cumulativeTime, overhead, allPatterns);
						long upperScore = function1(upper, cumulativeTime, overhead, allPatterns);
						long currentScore = function1(timeThreshold, cumulativeTime, overhead, allPatterns);

						bool scoreChanged = false;

						if(lowerScore<currentScore)
						{
							timeThreshold = lower;
							currentScore = lowerScore;
							scoreChanged = true;
							changeUpper = changeUpper / 2;
						}

						if(upperScore<currentScore)
						{
							timeThreshold = upper;
							currentScore = upperScore;
							scoreChanged = true;

							changeLower = changeLower / 2;
						}

						if(!scoreChanged)
						{
							changeUpper = changeUpper / 2;
							changeLower = changeLower / 2;

							if(changeUpper<1 && changeLower<1)
								break;
						}
					}
				}

				if(Settings::debugMSG)
					cout<<"The selected timeThreshold is: "<<timeThreshold<<endl;

				//use the threshold to assign how many tasks per pattern is needed
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

							unsigned long int tempTime = pattern->getPredictedTime();

							if(Settings::debugMSG)
							{
								cout<<"tempTime = "<<tempTime<<", timeThreshold = "<<timeThreshold<<endl;
							}

							if(tempTime>timeThreshold*2)
							{
								int numTasks = ceil(tempTime/timeThreshold);

								if(Settings::debugMSG)
									cout<<"numTasks:"<<numTasks<<endl;
								pattern->setSubtasking(numTasks, this->numWorkers);
							}
						}
					}
				}
			}
			break;
		case 3:
			{
				divideFFun3(allPatterns);
			}
			break;
		case 4:
			{
				divideFFun4(allPatterns);
			}
			break;
		}
	}
}

void MinerX2::initWithoutApprox(string fileName, int graphType, int support, int numMachine, int nThreads)
{
	this->nThreads = nThreads;
	numWorkers = (numMachine)*nThreads;
	//create available workers list
	bool* availableWorker = new bool[numWorkers];
	for(int i=1;i<numWorkers;i++)
	{
		availableWorker[i] = true;
	}

	//first time: send task to all workers

	cout<<"MinerX2: Sending graph filename to all workers: "<<endl;

	int loadCount = 0;

	for(int i=1;i<numMachine;i++)
	{
		int frequency = support;
		char fileName_[500];
		sprintf(fileName_, "f%s,%d\n", fileName.c_str(), frequency);
		int cnt=strlen(fileName_)+1;
		MPI_Send(fileName_, cnt, MPI_CHAR, i, 0, MPI_COMM_WORLD);

		if(Settings::debugMSG)
			cout<<"Sent to worker: "<<i<<endl;

		for(int j=0;j<nThreads;j++)
			availableWorker[i*nThreads+j] = false;

		loadCount++;
	}
	if(Settings::debugMSG)
		cout<<"Finished sending to all workers."<<endl;

	this->support = support;

	//load the graph
	long long start = getmsofday();
	loadGraph(fileName, graphType, support, frequentEdges);
	long long end = getmsofday();
	long long elapsed = end - start;

	cout<<"Master Loading took "<<(elapsed/1000)<<" sec and "<<(elapsed%1000)<<" ms"<<endl;

	cout<<"Waiting for all workers to finish loading ..."<<endl;
	long long start1 = getmsofday();
	while(loadCount>0)
	{
		MPI_Status status;
		char str_message[2000];
		MPI_Recv(str_message, 2000, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		int srcThread = status.MPI_SOURCE*nThreads+status.MPI_TAG;


		if(Settings::debugMSG)
			printf("The following message is recieved at master: %s", str_message);

		//ack message recieved
		if(str_message[0]=='l')
		{
			start1 = getmsofday();

			if(Settings::debugMSG)
				printf("Master received ack from %d\n", srcThread);

			loadCount--;
		}
		else
		{
			cout<<"jsgadudgfu767"<<endl;
			exit(0);
		}
	}
	masterProcessingTime += (getmsofday() - start1);

	if(frequentEdges.getSize()==0)
	{
		cout<<"No frequent patterns found! Exiting"<<endl;
		exit(0);
	}

	Settings::graphLoadingTime = elapsed;

	if(Settings::debugMSG)
	{
		cout<<"#frequent edges  = "<<frequentEdges.getSize()<<endl;
		cout<<"Frequent edges:"<<endl;
		frequentEdges.print();
	}
}

int MinerX2::getNumAvailWorkers()
{
	int numWorkers = numMachines*nThreads;
	int count = 0;
	for(int i=nThreads;i<numWorkers;i++)
	{
		if(availableWorker[i])
			count++;
	}
	return count;
}

void MinerX2::startMining()
{
	numWorkers = numMachines*nThreads;
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

	//extend frequent edges into 2 edges graphs
	long long start = getmsofday();
	extendFreqEdges();
	long long end = getmsofday();

	masterProcessingTime += (getmsofday() - start1);
	//Start of New loop
	while(true)
	{
		//do we have busy workers?
		bool doWeHaveABusyWorker = false;
		if(Settings::debugMSG)
			cout<<"Busy workers: ";

		for(int i=nThreads;i<numWorkers;i++)
		{
			if(!availableWorker[i])
			{
				doWeHaveABusyWorker = true;
				if(Settings::debugMSG)
					cout<<i<<", ";
				else
					break;
			}
		}
		if(Settings::debugMSG)
			cout<<endl;

		if(Settings::debugMSG)
			cout<<"doWeHaveABusyWorker = "<<doWeHaveABusyWorker<<endl;

		//get a message from a worker
		if(doWeHaveABusyWorker)
		{
			MPI_Status status;
			char str_message[2000];
			MPI_Recv(str_message, 2000, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			if(Settings::debugMSG)
				cout<<"[master] a msg is recieved"<<endl<<flush;

			int srcThread = status.MPI_SOURCE*nThreads+status.MPI_TAG;
			start1 = getmsofday();

			if(Settings::debugMSG)
				cout<<"The following message is recieved at master: "<<str_message<<endl<<flush;

			//ack message recieved
			if(str_message[0]=='l')
			{
				start1 = getmsofday();

				if(Settings::debugMSG)
					printf("Master received ack from %d\n", srcThread);
				//get elapsed tie for loading the graph at the worker
				unsigned long time;
				sscanf(str_message, "l,%Lu", &time);
				graphloadingTime[srcThread] = time;

				masterProcessingTime += (getmsofday() - start1);

				for(int j=0;j<nThreads;j++)
					availableWorker[status.MPI_SOURCE*nThreads+j] = true;
				//availableWorker[srcThread] = true;
			}
			else if(str_message[0]=='t')
			{
				start1 = getmsofday();

				//get the pattern that was sent
				int candidID;
				unsigned long time;
				sscanf(str_message, "t%d,%Lu", &candidID, &time);
				map<int, Pattern* >::iterator iter = currentlyChecking.find(candidID);

				if(iter==currentlyChecking.end())
				{
					cout<<"Pattern#"<<candidID<<" is not found"<<endl;
				}
				else
				{
					Pattern* candidate = iter->second;

					candidate->setFrequency(support);

					if(Settings::debugMSG)
						cout<<"Predicted Time for candidate["<<candidate->getID()<<"] = '"<<candidate->getPredictedTime()<<"'"<<endl;
					currentlyChecking.erase(candidID);

					if(candidate->getInvalidCol()>-1)
					{
						//adding candidate again
						processed.removePattern(candidate);

						while(candidates.size()<=candidate->getSize())
						{
							candidates.push_back(new CLMap());
						}

						candidate->generatePrimaryGraphs();

						if(candidates.at(candidate->getSize())->getPattern(candidate)==0)
						{
							candidates.at(candidate->getSize())->addPattern(candidate);
						}

						if(Settings::debugMSG)
							cout<<"Pattern added to candidates AGAIN [2]!"<<endl;

						//reset mni
						if(Settings::debugMSG)
							cout<<"candidate->getPredictedTime() = "<<candidate->getPredictedTime()<<", this->timeThreshold = "<<this->timeThreshold<<endl;

						int numTasks = candidate->getPredictedTime()/this->timeThreshold;//numWorkers-1;

						candidate->setSubtasking(numTasks, this->numWorkers);
						candidate->setInvalidCol(-1, -1);
						candidate->postponeExpensiveNodes = false;

						//resetting invalid column for the predicted pattern
						Pattern* predictedPattern = candidate->predictedPattern;
						if(predictedPattern==0)
							predictedPattern = this->getPredictedCandidate(candidate);

						if(predictedPattern==NULL)
						{
							cout<<"Error7535"<<endl;
							exit(0);
						}
						else
						{
							predictedPattern->setInvalidCol(-1, -1);
						}
					}
					else
					{
						MinerX::numFreqs++;

						char temp[500];
						sprintf(temp, "%d_%d", srcThread, candidID);
						isFrerqtimes.insert(std::pair<string, unsigned long>(temp, time));

						if(Settings::debugMSG)
							printf("Master received isFreq[TRUE] from %d\n", srcThread);

						patternFoundFrequent(candidate);
					}
				}

				availableWorker[srcThread] = true;

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

				map<int, Pattern* >::iterator iter = currentlyChecking.find(candidID);
				if(iter==currentlyChecking.end())
				{
					cout<<"Pattern#"<<candidID<<" is not found"<<endl;
				}
				else
				{
					Pattern* candidate = iter->second;
					candidate->setFrequency(0);
					if(Settings::debugMSG)
						cout<<"Predicted Time for candidate["<<candidate->getID()<<"] = '"<<candidate->getPredictedTime()<<"'"<<endl;
					currentlyChecking.erase(candidID);

					char temp[500];
					sprintf(temp, "%d_%d", srcThread, candidID);
					isFrerqtimes.insert(std::pair<string, unsigned long>(temp, time));

					patternFoundInFrequent(candidate);
				}
				availableWorker[srcThread] = true;

				masterProcessingTime += (getmsofday() - start1);
			}//the case when a task is partitioned
			else if(str_message[0]=='m')
			{
				start1 = getmsofday();

				//get the pattern that was sent
				int candidID;
				unsigned long time;
				//message format: m(candidateid),(elapsed),mniCol[0],mniCol[1],...,mniCol[n],0
				char tempStr[1000];

				sscanf(str_message, "m%d,%Lu,%s", &candidID, &time, tempStr);

				int c = 0;
				while(true)
				{
					char temp[500];
					sprintf(temp, "%d_%d_%d", srcThread, candidID, c);
					if(isFrerqtimes.find(temp)==isFrerqtimes.end())
					{
						isFrerqtimes.insert(std::pair<string, unsigned long>(temp, time));
						break;
					}
					c++;
				}

				map<int, Pattern* >::iterator iter = currentlyChecking.find(candidID);
				if(iter==currentlyChecking.end())
				{
					if(Settings::debugMSG)
						cout<<"Pattern#"<<candidID<<" is not found"<<endl;
				}
				else
				{
					Pattern* candidate = iter->second;
					//populate the mniTable
					int* mniTable = new int[candidate->getGraph()->getNumOfNodes()];

					char * brkth;
					char* pch = strtok_r (tempStr,",", &brkth);

					mniTable[0] = atoi(pch);

					for(int i=1;i<candidate->getGraph()->getNumOfNodes();i++)
					{
						pch = strtok_r (NULL, ",", &brkth);
						mniTable[i] = atoi(pch);
					}

					//populate the postponed_mniTable
					int* ppn_mniTable = new int[candidate->getGraph()->getNumOfNodes()];
					pch = strtok_r (NULL, ",", &brkth);
					ppn_mniTable[0] = atoi(pch);
					for(int i=1;i<candidate->getGraph()->getNumOfNodes();i++)
					{
						pch = strtok_r (NULL, ",", &brkth);
						ppn_mniTable[i] = atoi(pch);
						GraMiCounter::numPostponedNodes+=ppn_mniTable[i];
					}

					int freq = candidate->subTaskDone(mniTable, ppn_mniTable, this->support);

					delete[] mniTable;
					delete[] ppn_mniTable;

					if(Settings::debugMSG)
					{
						cout<<"Recieved result for subtask for candidate["<<candidate->getID()<<", freq = "<<freq<<endl<<flush;
						printf("Master received subtask from %d\n", srcThread);
					}

					//>-1 means that all subtasks for this pattern are done, and you can rely on the resulted frequency
					if(freq>-1)
					{
						currentlyChecking.erase(candidID);

						if(Settings::debugMSG)
							cout<<"Candidate: "<<candidID<<" is erased from currentlyChecking."<<endl;
						if(candidate->getInvalidCol()==-1)
						{
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
						else
						{
							if(freq>=support)
							{
								//adding candidate again
								processed.removePattern(candidate);

								while(candidates.size()<=candidate->getSize())
								{
									candidates.push_back(new CLMap());
								}

								candidate->generatePrimaryGraphs();

								if(candidates.at(candidate->getSize())->getPattern(candidate)==0)
								{
									candidates.at(candidate->getSize())->addPattern(candidate);
								}

								if(Settings::debugMSG)
									cout<<"Pattern added to candidates AGAIN!"<<endl;

								//reset mni
								int numTasks = candidate->getPredictedTime()/this->timeThreshold;//numWorkers-1;
								int nAvailWorkers = getNumAvailWorkers();
								if(nAvailWorkers>numTasks)
									numTasks = nAvailWorkers;

								candidate->setSubtasking(numTasks, this->numWorkers);
								candidate->setInvalidCol(-1, -1);

								//resetting invalid column for the predicted pattern
								Pattern* predictedPattern = candidate->predictedPattern;
								if(predictedPattern==0)
									predictedPattern = this->getPredictedCandidate(candidate);

								if(predictedPattern==NULL)
								{
									cout<<"Error6565"<<endl;
									exit(0);
								}
								else
								{
									predictedPattern->setInvalidCol(-1, -1);
								}
							}
							else
							{
								patternFoundInFrequent(candidate);
								if(Settings::debugMSG)
									cout<<"Pattern (using invalidCol) is added to infrequents"<<endl;
							}
						}
					}
					//the case when postponed Nodes require us to repeat processing
					if(freq==-2)
					{
						//adding candidate again
						if(Settings::debugMSG)
						{
							cout<<"Postponed nodes are used, but it affects the result, need to re-check AGAIN without the postponed nodes option."<<endl;
						}

						candidate->postponeExpensiveNodes = false;
						processed.removePattern(candidate);

						while(candidates.size()<=candidate->getSize())
						{
							candidates.push_back(new CLMap());
						}

						candidate->generatePrimaryGraphs();

						if(candidates.at(candidate->getSize())->getPattern(candidate)==0)
						{
							candidates.at(candidate->getSize())->addPattern(candidate);
						}

						if(Settings::debugMSG)
							cout<<"Pattern added to candidates AGAIN!"<<endl;

						//reset mni
						int numTasks = candidate->getSubtaskingValueFixed();
						int nAvailWorkers = getNumAvailWorkers();
						if(nAvailWorkers>numTasks)
							numTasks = nAvailWorkers;

						candidate->setSubtasking(numTasks, this->numWorkers);
						//candidate->setInvalidCol(-1, -1);
					}
				}

				availableWorker[srcThread] = true;

				masterProcessingTime += (getmsofday() - start1);
			}
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

					map<int, Pattern* >::iterator iter = currentlyChecking.find(candidID);

					if(iter==currentlyChecking.end())
					{
						cout<<"Pattern#"<<candidID<<" is not found"<<endl;
					}
					else
					{
						Pattern* candidate = iter->second;
						candidate->setFrequency(frequency);
						candidate->setPredictedTime(predictedTime);
						candidate->setMaxIters(numIters);
						if(exact==1)
						{
							if(Settings::debugMSG)
								cout<<"Set result exact for this pattern"<<endl;
							candidate->setResultExact();
							candidate->setFrequency(support);
						}
						//candidate->setSubtasking(predictedTime/timeThreshold);

						unsigned long int tempTime = candidate->getPredictedTime();

						if(tempTime>timeThreshold*2)
						{
							int numTasks = ceil(tempTime/timeThreshold);
							if(numTasks>(numWorkers-1))
								numTasks = (numWorkers-1);
							if(Settings::debugMSG)
								cout<<"numTasks:"<<numTasks<<endl;
							candidate->setSubtasking(numTasks, this->numWorkers);
						}
						else
						{
							candidate->setSubtasking(1, this->numWorkers);
						}

						if(Settings::debugMSG)
							cout<<"Predicted Time: "<<predictedTime<<" for candidateID: "<<candidID<<endl;
						currentlyChecking.erase(candidID);

						char temp[500];
						sprintf(temp, "%d_%d", srcThread, candidID);
						isFrerqtimes.insert(std::pair<string, unsigned long>(temp, time));

						if(Settings::debugMSG)
						{
							printf("Master received isFreqApprox[TRUE] from %d.\n", srcThread);
							cout<<"Found to be FREQUENT, with frequency = "<<candidate->getFrequency()<<endl;
						}
						//add to the frequent patterns list
						while((expected_frequentPatterns->size())<=candidate->getSize())
						{
							expected_frequentPatterns->push_back(new CLMap());
						}
						expected_frequentPatterns->at(candidate->getSize())->addPattern(candidate);

						//*********************NEW ********************
						//set the original candidate available for processing
						CLMap* tempListCandids = candidates.at(candidate->getSize());
						CLMap_Iterator iter = tempListCandids->getFirstElement();

						while(iter.pattern!=0)
						{
							Pattern* tempCandidate = iter.pattern;//->second;
							string patternKey = iter.key;//->first;

							if(tempCandidate->getGraph()->getCanonicalLabel().compare(candidate->getGraph()->getCanonicalLabel())==0)
							{
								tempCandidate->inTheProcessOfPrediction = false;
								break;
							}

							tempListCandids->advanceIterator(iter);
						}
						//***************NEW END*****************


						if(Settings::divideByFunc==3)
						{
							if(Settings::debugMSG)
							{
								cout<<"Reevaluating division threshold"<<endl;
							}
							vector<CLMap* >* allPatterns[2];
							allPatterns[0] = expected_frequentPatterns;
							allPatterns[1] = expected_infrequentPatterns;
							divideFFun3(allPatterns);
						}
					}

					availableWorker[srcThread] = true;
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

					map<int, Pattern* >::iterator iter = currentlyChecking.find(candidID);

					if(iter==currentlyChecking.end())
					{
						cout<<"Pattern#"<<candidID<<" is not found"<<endl;
					}
					else
					{
						Pattern* candidate = iter->second;
						candidate->setFrequency(frequency);
						candidate->setPredictedTime(predictedTime);
						candidate->setInvalidCol(invalidCol, predictedValids);
						candidate->setMaxIters(numIters);
						if(exact==1)
						{
							if(Settings::debugMSG)
								cout<<"Set result exact for this pattern"<<endl;
							candidate->setResultExact();
							candidate->setFrequency(0);
						}

						if(Settings::debugMSG)
							cout<<"Predicted Time: "<<predictedTime<<" for candidateID: "<<candidID<<", invalid col: "<<invalidCol<<" with #valid nodes = "<<predictedValids<<endl;

						unsigned long tempTime = candidate->getPredictedTime();

						if(tempTime>timeThreshold*2)
						{
							int numTasks = ceil(tempTime/timeThreshold);

							if(invalidCol==-1)
							{
								if(numTasks>(numWorkers-1))
									numTasks = (numWorkers-1);
							}
							else
							{
								if(numTasks>((numWorkers-1)*2))
									numTasks = ((numWorkers-1)*2);
							}

							if(Settings::debugMSG)
								cout<<"numTasks:"<<numTasks<<endl;
							candidate->setSubtasking(numTasks, this->numWorkers);
						}
						else
						{
							candidate->setSubtasking(1, this->numWorkers);
						}

						currentlyChecking.erase(candidID);

						char temp[500];
						sprintf(temp, "%d_%d", srcThread, candidID);
						isFrerqtimes.insert(std::pair<string, unsigned long>(temp, time));

						while(expected_infrequentPatterns->size()<=candidate->getSize())
						{
							expected_infrequentPatterns->push_back(new CLMap());
						}
						expected_infrequentPatterns->at(candidate->getSize())->addPattern(candidate);

						//*********************NEW ********************
						//set the original candidate available for processing
						CLMap* tempListCandids = candidates.at(candidate->getSize());
						CLMap_Iterator iter = tempListCandids->getFirstElement();

						while(iter.pattern!=0)
						{
							Pattern* tempCandidate = iter.pattern;//->second;
							string patternKey = iter.key;//->first;

							if(tempCandidate->getGraph()->getCanonicalLabel().compare(candidate->getGraph()->getCanonicalLabel())==0)
							{
								tempCandidate->inTheProcessOfPrediction = false;
								break;
							}

							tempListCandids->advanceIterator(iter);
						}
						//***************NEW END*****************

						if(Settings::divideByFunc==3)
						{
							if(Settings::debugMSG)
							{
								cout<<"Reevaluating division threshold"<<endl;
							}
							vector<CLMap* >* allPatterns[2];
							allPatterns[0] = expected_frequentPatterns;
							allPatterns[1] = expected_infrequentPatterns;
							divideFFun3(allPatterns);
						}
					}

					availableWorker[srcThread] = true;
				}
			}
		}

		long start1 = getmsofday();

		for(int i=nThreads;i<numWorkers;i++)
		{
			if(!availableWorker[i])
			{
				if(Settings::debugMSG)
					cout<<"worker: "<<i<<" is not available. numWorkers = "<<numWorkers<<endl;
				continue;
			}
			else
			{
				if(Settings::debugMSG)
					cout<<"worker: "<<i<<" is expecting a task."<<endl;
			}

			string patternKey;
			int evalType; //0->exact, 1->approximate
			Pattern* candidate = selectCandidate(patternKey, evalType);

			//send the selected candidate to a worker
			if(candidate)
			{
				sendACandidate(patternKey, candidate, currentlyChecking, i, evalType);
			}
			else
			{
				if(Settings::debugMSG)
					cout<<"No task for this worker "<<i<<" is available"<<endl;
				break;
			}
		}

		MinerX2::sendCandidsTime+=(getmsofday()-start1);

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
//				cout<<"#CANDIDATES = "<<getNumElems((vector<map<string, void*>* >*)&candidates)<<", #currentlyChecking = "<<currentlyChecking.size()<<endl;
//			}
		}

		//we are done, we cannot get more candidates
		if(currentlyChecking.size()==0 && getNumElems((vector<map<string, void*>* >*)&candidates)==0)
		{
			masterProcessingTime += (getmsofday() - start1);
			break;
		}

		masterProcessingTime += (getmsofday() - start1);
	}
	//End of new loop

	cout<<"Mining Finished, results:"<<endl;
	printResult();

	delete availableWorker;
	delete graphloadingTime;
	delete processingTime;
}

void MinerX2::sendACandidate(string key, Pattern* candidate, map<int, Pattern*>& currentlyChecking, int destination, bool evalType)
{
	if(Settings::debugMSG)
		cout<<"Worker "<<destination<<" will be busy!"<<endl;

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
	char c;
	int cID = candidate->getID();
	if(evalType==0)
	{
		c = 's';
		sprintf(graphStr, "%c,%d,%d,%d,%d,%d,%d,%Lu,%s\t", c, cID, candidate->getInvalidCol(), candidate->getPredictedValids(), partToProcess, numParts, candidate->postponeExpensiveNodes, candidate->getMaxIters(), tempOS.str().c_str());
	}
	else
	{
		c = 'a';
		sprintf(graphStr, "%c,%d,%s\t", c, cID, tempOS.str().c_str());
		//send approximate requests with negative IDs (just not over-write exact patterns) ... [no need now!]
//		if(cID>0)
//			cID = cID*-1;
	}
	int cnt=strlen(graphStr)+1;
	MPI_Status status;
	int rank = destination/nThreads;
	int tagID = destination%nThreads;
	MPI_Send(graphStr, cnt, MPI_CHAR, rank, tagID, MPI_COMM_WORLD);

	if(Settings::debugMSG)
	{
		printf("Sending subgraph[%d] [%s] to worker: %d for frequency check\n", candidate->getID(), candidate->getGraph()->getCanonicalLabel().c_str(),destination);
		printf("Message is: [%s]\n", graphStr);
	}

	//add it to processed only when all subtasks are done
	if(evalType==0 && candidate->getSubtaskingValue()==0)//I added the evalType==0 as a new condition (2/12/2015)
		processed.addPattern(candidate);

	availableWorker[destination] = false;

	MinerX::numIsFreqCalls++;
}

Pattern* MinerX2::getPredictedCandidate(Pattern* candidate)
{
	vector<CLMap*>* expectedLists[2];
	expectedLists[0] = expected_frequentPatterns;
	expectedLists[1] = expected_infrequentPatterns;

	for(int i=0;i<2;i++)
	{
		if(expectedLists[i] && expectedLists[i]->size()>candidate->getSize() && expectedLists[i]->size()>0)
		{
			//for(vector<CLMap*>::iterator iterV = expectedLists[i]->begin();iterV!=expectedLists[i]->end();iterV++)
			{
				CLMap* tempMap = expectedLists[i]->at(candidate->getSize());// (*iterV);
				CLMap_Iterator iterM = tempMap->getFirstElement();

				//for(map<string, Pattern*>::iterator iterM = tempMap->begin();iterM!=tempMap->end();iterM++)
				while(iterM.pattern!=0)
				{
					Pattern* tempPattern = iterM.pattern;//->second;
					tempMap->advanceIterator(iterM);

					if(tempPattern->getGraph()->isTheSameWith(candidate->getGraph()))
						return tempPattern;
				}
			}
		}
	}

	return NULL;
}

void MinerX2::patternFoundInFrequent(Pattern* candidate)
{
	long long start1 = getmsofday();

	while(infrequentPatterns.size()<=candidate->getSize())
	{
		infrequentPatterns.push_back(new CLMap());
	}

	infrequentPatterns[candidate->getSize()]->addPattern(candidate);

	//remove its supergraphs from the set of candidates (if there exist)
	for(int i=candidate->getSize()+1;i<candidates.size();i++)
	{
		CLMap* currentList = candidates.at(i);
		CLMap_Iterator iter = currentList->getFirstElement();

		//currentList->print();

		bool needToInspectHigherLevels = false;

		while(iter.pattern!=0)
		{
			Pattern* currentcandid = iter.pattern;

			if(currentcandid->getGraph()->getNumOfNodes()<candidate->getGraph()->getNumOfNodes())
			{
				currentList->advanceIterator(iter);
				continue;
			}

			vector<map<int, int>* > result;
			tr1::unordered_map<int, tr1::unordered_set<int>*> domains_values;

			//do isomorphism only when node labels are satisfied
			if(currentcandid->getGraph()->CanSatisfyNodeLabels(candidate->getGraph()))
				candidate->getGraph()->isIsomorphic(currentcandid->getGraph(), result, domains_values, -1, -1, false);

			if(result.size()>0)
			{
				currentList->advanceIterator(iter);
				currentList->removePattern(currentcandid);

				this->currentlyChecking.erase(currentcandid->getID());

				if(Settings::debugMSG)
				{
					cout<<"Erasing candidate# "<<currentcandid->getID()<<endl;
					cout<<"iter.pattern = "<<iter.pattern<<endl;
				}

				delete currentcandid;

				needToInspectHigherLevels = true;

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

			if(!needToInspectHigherLevels)
				break;
		}
	}

	MinerX2::patternInfreqTime+=(getmsofday()-start1);
}

bool MinerX2::checkAgainstOtherPatterns(Pattern* newCandidate)
{
	//do not add it if it is already done
	if(processed.getPattern(newCandidate)==0)
		;
	else
	{
		if(Settings::debugMSG)
			cout<<newCandidate->getID()<<" with CL: "<<newCandidate->getGraph()->getCanonicalLabel()<<" not added (1)!!!!!!!!!"<<endl;
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
	{
		if(Settings::debugMSG)
			cout<<newCandidate->getID()<<" with CL: "<<newCandidate->getGraph()->getCanonicalLabel()<<" not added (1)!!!!!!!!!"<<endl;
		return true;
	}

	return false;
}

Pattern* MinerX2::selectCandidate(string& bestPatternKey, int& evalType)
{
	//get a free candidate
	Pattern* candidate = 0;
	string patternKey;
	Pattern* bestCandidate = 0;

	//get a candidate from the list of exact candidates
	for(int i=0;i<candidates.size();i++)
	{
		//adaptive partitioning
		if(Settings::divideByFunc==2)
		{
			divideF2Counter++;
			if(divideF2Counter>10)
			{
				//redivide the remaining tasks
				int numSubtasks = this->getNumberOfExpectedSubTasks();
				cout<<"CXC. numSubtasks = "<<numSubtasks<<", this->numWorkers = "<<this->numWorkers<<endl;
				if(numSubtasks<this->numWorkers)
				{
					if(Settings::debugMSG)
					{
						cout<<"We need more subtasks"<<endl;
					}
				}

				divideF2Counter = 0;
			}
		}


		CLMap* tempListCandids = candidates.at(i);

		if(tempListCandids->getSize()>0)
		{
			CLMap_Iterator iter = tempListCandids->getFirstElement();

			while(iter.pattern!=0)
			{
				candidate = iter.pattern;
				patternKey = iter.key;

				tempListCandids->advanceIterator(iter);

				//special case, if this pattern a one edge candidate, then there is no need to evaluate it as it is already frequent
				if(candidate->getGraph()->getNumOfEdges()==1)
				{
					tempListCandids->removePattern(candidate);
					if(Settings::debugMSG)
						cout<<"Pattern: "<<candidate->getID()<<", has been removed as it is a graph with ONE edge."<<endl;

					candidate = 0;
					continue;
				}

				//if this pattern is under prediction
				if(candidate->inTheProcessOfPrediction)
				{
					candidate=0;
					continue;
				}

				//in case it is already processed before, remove it and continue
				if(processed.getPattern(candidate)!=0)
				{
					tempListCandids->removePattern(candidate);
					if(Settings::debugMSG)
						cout<<"Pattern: "<<candidate->getID()<<", has been removed as it is seen in processed."<<endl;

					candidate = 0;
					continue;
				}

				//search for its prediction information. If not available, conduct a sampling-based execution.
				long long start11 = getmsofday();
				Pattern* predictedPattern = candidate->predictedPattern;
				if(predictedPattern==0)
					predictedPattern = this->getPredictedCandidate(candidate);

				MinerX2::predictedPatternTime+=(getmsofday()-start11);

				if(Settings::predictOnTheGo)
				{
					if(predictedPattern==NULL)
					{
						if(Settings::debugMSG)
							cout<<"Pattern: "<<candidate->getID()<<", has not been predicted."<<endl;

						candidate->inTheProcessOfPrediction = true;
						//case that this pattern has not been predicted
						candidate = new Pattern(candidate);
						candidate->makeIDNegative();

						evalType = 1;

						bestCandidate = candidate;
						bestPatternKey = patternKey;

						if(Settings::debugMSG)
							cout<<"Candidate#"<<candidate->getID()<<" is sent for estimations ...."<<endl;
						break;
					}
					else
					{
						predictedPattern->selected = true;
						//case that this pattern is already predicted
						evalType = 0;
						candidate->borrowTimeInfor(predictedPattern, this->numWorkers);
						candidate->predictedPattern = predictedPattern;

					}
				}
				else
				{
					evalType = 0;
					candidate->setPredictedTime(-1);

					candidate->setSubtasking(1, this->numWorkers);
					candidate->setInvalidCol(-1, -1);

					candidate->setMaxIters(-1);
				}

				bestCandidate = candidate;
				bestPatternKey = patternKey;

				break;
			}
		}

		if(bestCandidate!=0)
		{
			if(bestCandidate->isResultExact())
			{
				if(Settings::debugMSG)
					cout<<"Pattern# "<<bestCandidate->getID()<<" was exactly approximated! we can use its frequency."<<endl;
				if(bestCandidate->getFrequency()>=support)
				{
					if(Settings::debugMSG)
						cout<<"FOUND FREQUENT"<<endl;
					patternFoundFrequent(bestCandidate);
				}
				else
				{
					if(Settings::debugMSG)
						cout<<"FOUND INFREQUENT"<<endl;
					patternFoundInFrequent(bestCandidate);
				}

				if(Settings::debugMSG)
					cout<<"Removing candidate#"<<bestCandidate->getID()<<" from the list of candidates."<<endl;
				tempListCandids->removePattern(bestCandidate);

				bestCandidate=0;
				i--;
				continue;
			}
			else
			{
				bestCandidate->setSubtaskTaken();
				if(Settings::debugMSG)
					cout<<"candidate->subtasking: "<<bestCandidate->getSubtaskingValue()<<"/"<<bestCandidate->getSubtaskingValueFixed()<<" for candidate ID: "<<bestCandidate->getID()<<endl;
			}

			//remove it from the list of exact candidates (if it has no more subtasks)
			if(bestCandidate->subtasksFinished() && evalType==0)
			{
				if(Settings::debugMSG)
					cout<<"Removing candidate#"<<bestCandidate->getID()<<" from the list of candidates."<<endl;
				tempListCandids->removePattern(bestCandidate);
			}

			break;
		}
	}

	return bestCandidate;
}

long MinerX2::function1(long timeThreshold, long cumulativeTime, double overhead, vector<CLMap* >** allPatterns)
{
	long subTasks = 0;

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

				double tempTime = pattern->getPredictedTime();

				if(tempTime<timeThreshold)
					subTasks++;
				else
					subTasks+=ceil(tempTime/timeThreshold);
			}
		}
	}

	long score = timeThreshold+cumulativeTime*subTasks*overhead;

	if(Settings::debugMSG)
	{
		cout<<"cumulativeTime = "<<cumulativeTime<<", overhead = "<<overhead<<", subTasks = "<<subTasks<<endl;
		cout<<"For timeThreshold = "<<timeThreshold<<", score is: "<<score<<endl;
	}

	return score;
}

int MinerX2::getNumberOfExpectedSubTasks()
{
	int counter = 0;

	vector<CLMap*>* expectedLists[2];
	expectedLists[0] = expected_frequentPatterns;
	expectedLists[1] = expected_infrequentPatterns;

	for(int i=0;i<2;i++)
	{
		if(expectedLists[i] && expectedLists[i]->size()>0)
		{
			for(vector<CLMap*>::iterator iterV = expectedLists[i]->begin();iterV!=expectedLists[i]->end();iterV++)
			{
				CLMap* tempMap = (*iterV);
				CLMap_Iterator iterM = tempMap->getFirstElement();

				while(iterM.pattern!=0)
				{
					Pattern* tempPattern = iterM.pattern;//->second;
					tempMap->advanceIterator(iterM);

					if(tempPattern->getSubtaskingValue()==1)
					{
						counter+=tempPattern->getSubtaskingValueFixed();
					}
				}
			}
		}
	}

	return counter;
}

void MinerX2::divideFFun3(vector<CLMap* >* allPatterns[2])
{
	int totalExpectedTime = 0;
	for(int l=0;l<2;l++)
	{
		for(vector<CLMap* >::iterator iter = allPatterns[l]->begin();iter!=allPatterns[l]->end();iter++)
		{
			CLMap* tempMap = (*iter);
			CLMap_Iterator iter1 = tempMap->getFirstElement();
			while(iter1.pattern!=0)
			{
				Pattern* pattern = iter1.pattern;
				tempMap->advanceIterator(iter1);

				if(pattern->selected)
					continue;

				unsigned long int tempTime = pattern->getPredictedTime();

				totalExpectedTime+=tempTime;
			}
		}
	}

	timeThreshold = totalExpectedTime/(2*(this->numWorkers-1));
	if(Settings::debugMSG)
	{
		cout<<"The total expected time is: "<<totalExpectedTime<<endl;
		cout<<"The selected timeThreshold is: "<<timeThreshold<<endl;
	}

	//use the threshold to assign how many tasks per pattern is needed
	for(int l=0;l<2;l++)
	{
		for(vector<CLMap* >::iterator iter = allPatterns[l]->begin();iter!=allPatterns[l]->end();iter++)
		{
			CLMap* tempMap = (*iter);
			CLMap_Iterator iter1 = tempMap->getFirstElement();
			while(iter1.pattern!=0)
			{
				Pattern* pattern = iter1.pattern;
				tempMap->advanceIterator(iter1);

				if(pattern->selected)
					continue;

				unsigned long int tempTime = pattern->getPredictedTime();

				if(Settings::debugMSG)
				{
					cout<<"Pattern#"<<pattern->getID()<<endl;
					cout<<"tempTime = "<<tempTime<<", timeThreshold = "<<timeThreshold<<endl;
				}

				if(tempTime>timeThreshold*2)
				{
					int numTasks = ceil(tempTime/timeThreshold);

					if(numTasks>(numWorkers-1))
						numTasks = (numWorkers-1);

					if(Settings::debugMSG)
						cout<<"numTasks:"<<numTasks<<endl;

					//reset subtasking
					pattern->setSubtasking(1, this->numWorkers);

					pattern->setSubtasking(numTasks, this->numWorkers);
				}
			}
		}
	}
}

void MinerX2::divideFFun4(vector<CLMap* >* allPatterns[2])
{
	//iterate to get the total time
	int totalNumOfPatterns = 0;
	for(vector<CLMap*>::iterator iter = candidates.begin();iter!=candidates.end();iter++)
	{
		totalNumOfPatterns+=(*iter)->getSize();
	}

	long* times = new long[totalNumOfPatterns];
	int pointer = 0;

	if(Settings::debugMSG)
		cout<<"Number of predicted times = "<<totalNumOfPatterns<<endl;

	//iterate over candidates, ordered by size
	//add each one to the times array
	for(int i=0;i<candidates.size();i++)
	{
		CLMap* tempListCandids = new CLMap();
		tempListCandids->addAll(candidates.at(i));

		while(tempListCandids->getSize()>0)
		{
			CLMap_Iterator iter = tempListCandids->getFirstElement();
			Pattern* candidate = iter.pattern;//->second;
			string patternKey = iter.key;//->first;

			Pattern* predictedPattern = candidate->predictedPattern;
			if(predictedPattern==0)
				predictedPattern = this->getPredictedCandidate(candidate);

			if(predictedPattern==NULL)
			{
				if(Settings::debugMSG)
					cout<<"This candidate has no predictions, setting its predicted time to 0";
				candidate->setPredictedTime(0);
			}
			else
			{
				candidate->borrowTimeInfor(predictedPattern, this->numWorkers);
			}

			times[pointer] = candidate->getPredictedTime();
			if(Settings::debugMSG)
				cout<<"Adding to times: "<<candidate->getPredictedTime()<<" at pointer = "<<pointer<<endl;
			pointer++;

			tempListCandids->removePattern(candidate);
		}

		delete tempListCandids;
	}

	//do simulation to get the dividing threshold
	timeThreshold = QueueSimulator::simulate(times, totalNumOfPatterns, (this->numMachines-1)*this->nThreads);
	if(Settings::debugMSG)
	{
		cout<<"[SIMULATION BASED] The selected timeThreshold is: "<<timeThreshold<<endl;
	}

	//use the threshold to assign how many tasks per pattern is needed
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

				if(pattern->selected)
					continue;

				unsigned long int tempTime = pattern->getPredictedTime();

				if(Settings::debugMSG)
				{
					cout<<"Pattern#"<<pattern->getID()<<endl;
					cout<<"tempTime = "<<tempTime<<", timeThreshold = "<<timeThreshold<<endl;
				}

				if(tempTime>timeThreshold*2)
				{
					int numTasks = ceil(tempTime/timeThreshold);

					if(numTasks>((this->numWorkers-1)*this->nThreads))
						numTasks = ((this->numWorkers-1)*this->nThreads);

					if(Settings::debugMSG)
						cout<<"numTasks:"<<numTasks<<endl;

					//reset subtasking
					pattern->setSubtasking(1, this->numWorkers);

					pattern->setSubtasking(numTasks, this->numWorkers);
				}
			}
		}
	}

	delete[] times;
}

void MinerX2::setSupport(int numMachines, int support)
{
	this->support = support;

	if(Settings::debugMSG)
		cout<<"Sending new support value to workers ..."<<endl;
	for(int i=nThreads;i<numMachines*nThreads;i++)
	{
		int frequency = support;
		char supportStr[500];
		sprintf(supportStr, "t%d\n", frequency);
		int cnt=strlen(supportStr)+1;
		int rank = i/nThreads;
		int tagID = i%nThreads;
		MPI_Send(supportStr, cnt, MPI_CHAR, rank, tagID, MPI_COMM_WORLD);

		if(Settings::debugMSG)
			cout<<"New support sent to worker: "<<i<<endl;
	}
}
