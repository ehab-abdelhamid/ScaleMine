/**
 * Copyright 2016 Ehab Abdelhamid, Ibrahim Abdelaziz

Worker.cpp

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
 * represents a compute node
 */
#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include <sstream>
#include<sys/time.h>
#include "Worker.h"
#include "utils.h"
#include "GraMiCounter.h"
#include "Settings.h"
#include <boost/thread.hpp>

#include <unistd.h>

void clearAllMP(tr1::unordered_map<string, void* >& allMP)
{
	for(tr1::unordered_map<string, void* >::iterator iter = allMP.begin();iter!=allMP.end();iter++)
	{
		delete iter->second;
	}
	allMP.clear();
}

void Worker::start(int machineID, int nThreads)
{
	this->machineID = machineID;
	this->nThreads = nThreads;
	char str_message[BUFFERSIZE];
	boost::thread** threads;

	//initialize
	//get filename from master
	if(Settings::debugMSG)
		cout<<"Worker is waiting for messages ..."<<endl;

	MPI_Status status;
	MPI_Recv(str_message, BUFFERSIZE, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

	if(Settings::debugMSG)
		cout<<"Worker recieved load msg:"<<str_message<<endl<<flush;

	if(str_message[0]=='f')
	{
		unsigned long startTime = getmsofday();
		int pos = getPos(str_message, ',');
		char fileName[500];
		strncpy (fileName, str_message+1, pos-1);
		fileName[pos-1] = '\0';
		char frequencyStr[10];
		strncpy (frequencyStr, str_message+pos+1, strlen(str_message)-pos);
		frequency = atoi(frequencyStr);
		if(Settings::debugMSG)
			printf("Worker %d: Recieved File Name: %s, Frequency %d\n", machineID, fileName, frequency);
		graph = new GraphX(1, 0);//create a graph with id=its partition id
		graph->setFrequency(frequency);
		string partFileName(fileName);
		tr1::unordered_map<string, void* > allMPs;
		if(!graph->loadFromFile(partFileName, allMPs))
		{
			printf("Graph cannot be loaded at worker: %d\n",machineID);
			delete graph;
			return;
		}
		unsigned long elapsed = getmsofday() - startTime;
		//send notification to master
		if(Settings::debugMSG)
			printf("Graph loaded at worker: %d\n",machineID);

		char str_message_1[BUFFERSIZE];

		sprintf(str_message_1, "l,%Lu",elapsed);
		int cnt=strlen(str_message_1)+1;

		//initialize threads
		threads = new boost::thread*[nThreads];
		for(int i=0;i<nThreads;i++)
		{
			try
			{
				threads[i] = new boost::thread(&Worker::processing, this, i);
			}
			catch(boost::thread_exception ex)
			{
				cout<<ex.what()<<endl;
				cout<<"Thread Exception"<<endl<<flush;
			}
		}

		MPI_Send(str_message_1, cnt, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
	}
	else
	{
		cout<<"uqskgkjsdfk"<<endl;
		exit(0);
	}

	for(int i=0;i<nThreads;i++)
	{
		threads[i]->join();
	}
}

void *Worker::processing(Worker* worker, int id)
{
	usleep(1000000);
	int threadID = worker->machineID*worker->nThreads+id;
	//receive task
	while(1)
	{
		if(Settings::debugMSG)
			cout<<"Worker is waiting for messages ..."<<endl;

		char* str_message = new char[BUFFERSIZE+1];
		MPI_Status status;
		MPI_Recv(str_message, BUFFERSIZE, MPI_CHAR, 0, id, MPI_COMM_WORLD, &status);

		if(Settings::debugMSG)
			cout<<"Worker recieved msg:"<<str_message<<endl<<flush;

		//I recieve a subgraph for isfrequency (exact)
		if(str_message[0]=='s')
		{
			tr1::unordered_map<string, void* > allMPs;
			//load graph data
			std::stringstream sstmGF;
			int candidateID;
			int invalidCol;
			int predictedValids;
			char subgraphTemp[500];
			int subTaskNum;
			int subTaskMax;
			int postponeExpensiveNodes;
			unsigned long avgItersPerNode;
			sscanf (str_message,"s,%d,%d,%d,%d,%d,%d,%Lu,%[^\t]",&candidateID,&invalidCol,&predictedValids,&subTaskNum,&subTaskMax,&postponeExpensiveNodes,&avgItersPerNode,subgraphTemp);
			if(Settings::debugMSG)
				printf("Worker %d: Recieved IsFreq for subgraph ID: %d. InvalidCol = %d, predicted valids = %d. It is asked to process part: %d/%d. Using postponedNodes? %d\n", threadID, candidateID, invalidCol, predictedValids, subTaskNum, subTaskMax, postponeExpensiveNodes);
			sstmGF << subgraphTemp;
			string subgraphStr = sstmGF.str();
			GraphX* subgraph = new GraphX(1, false);//create a graph with id=its partition id
			subgraph->loadFromString(subgraphStr, allMPs);
			clearAllMP(allMPs);
			if(Settings::debugMSG)
				cout<<"Worker "<<threadID<<" loaded pattern "<<candidateID<<" with #nodes = "<<subgraph->getNumOfNodes()<<" and #edges = "<<subgraph->getNumOfEdges()<<endl;
			Pattern* newCandidate = new Pattern(subgraph, false);
			newCandidate->setInvalidCol(invalidCol, predictedValids);
			bool isFreq;
			unsigned long startTime = getmsofday();

			if(Settings::debugMSG)
				printf("Worker %d: Starts checking frequency\n", threadID);

			int* mniTable = new int[subgraph->getNumOfNodes()];

			tr1::unordered_map<int, tr1::unordered_set<int>* >* postponedNodes = 0;
			if(postponeExpensiveNodes)
				postponedNodes = new tr1::unordered_map<int, tr1::unordered_set<int>* >();

			if(postponeExpensiveNodes && invalidCol==-1)
			{
				cout<<"Error: cannot postpone nodes while invalid column is not set."<<endl<<flush;
				exit(0);
			}

			long maxIters = avgItersPerNode*10;
			if(Settings::debugMSG)
				cout<<"avgItersPerNode = "<<avgItersPerNode<<", maxIters = "<<maxIters<<endl;

			isFreq = worker->isFrequent(newCandidate, subTaskNum, subTaskMax, mniTable, postponedNodes, maxIters);

			//decide whether we need to conduct it again or not
			if(subTaskNum==-1)
			{
				int freq1 = mniTable[0];
				for(int i=1;i<subgraph->getNumOfNodes();i++)
				{
					if(freq1<mniTable[i])
						freq1 = mniTable[i];
				}

				int value = 0;
				if(postponedNodes!=0)
				{
					tr1::unordered_map<int, tr1::unordered_set<int>* >::iterator iter = postponedNodes->find(0);
					if(iter!=postponedNodes->end())
						value = iter->second->size();
				}

				int freq2 = mniTable[0]+value;

				for(int i=1;i<subgraph->getNumOfNodes();i++)
				{
					int value = 0;
					if(postponedNodes!=0)
					{
						tr1::unordered_map<int, tr1::unordered_set<int>* >::iterator iter = postponedNodes->find(i);
						if(iter!=postponedNodes->end())
							value = iter->second->size();
					}
					if(freq2<mniTable[i]+value)
						freq2 = mniTable[i]+value;
				}


				if(freq1!=freq2)
				{
					if(isFreq)
					{
						isFreq = worker->isFrequent(newCandidate, subTaskNum, subTaskMax, mniTable, 0, maxIters);
					}
					else
					{
						if(freq2>=worker->frequency)
							isFreq = worker->isFrequent(newCandidate, subTaskNum, subTaskMax, mniTable, 0, maxIters);
					}
				}
			}

			if(Settings::debugMSG)
			{
				printf("Worker %d: Finished checking frequency. Result is: %d\n", threadID, isFreq);
				printf("Subtask Num:%d", subTaskNum);
			}

			unsigned long elapsed = getmsofday() - startTime;

			if(Settings::debugMSG)
				printf("Time taken is for candidate [%d]: %Lu\n", candidateID, elapsed);

			if(subTaskNum==-1)
			{
				if(isFreq)
				{
					sprintf(str_message, "t%d,%Lu", candidateID, elapsed);
					int cnt=strlen(str_message)+1;
					MPI_Send(str_message, cnt, MPI_CHAR, 0, id, MPI_COMM_WORLD);
				}
				else
				{
					sprintf(str_message, "f%d,%Lu", candidateID, elapsed);
					int cnt=strlen(str_message)+1;
					MPI_Send(str_message, cnt, MPI_CHAR, 0, id, MPI_COMM_WORLD);
				}
			}
			else
			{
				//message format: m(candidateid),(elapsed),mniCol[0],mniCol[1],...,mniCol[n],0
				sprintf(str_message, "m%d,%Lu", candidateID, elapsed);

				//put the mni table
				for(int i=0;i<subgraph->getNumOfNodes();i++)
				{
					sprintf(str_message, "%s,%d", str_message, mniTable[i]);
				}

				//put the postponed mni table
				for(int i=0;i<subgraph->getNumOfNodes();i++)
				{
					int value = 0;
					if(postponeExpensiveNodes!=0)
					{
						tr1::unordered_map<int, tr1::unordered_set<int>* >::iterator iter = postponedNodes->find(i);
						if(iter!=postponedNodes->end())
							value = iter->second->size();
					}

					sprintf(str_message, "%s,%d", str_message, value);
				}

				int cnt=strlen(str_message)+1;
				MPI_Send(str_message, cnt, MPI_CHAR, 0, id, MPI_COMM_WORLD);
			}

			if(Settings::debugMSG)
				cout<<"[Worker] sent a message: "<<str_message<<endl<<flush;

			delete[] mniTable;

			if(postponedNodes!=0)
			{
				for(tr1::unordered_map<int, tr1::unordered_set<int>* >::iterator iterr = postponedNodes->begin();iterr!=postponedNodes->end();iterr++)
				{
					iterr->second->clear();
					delete iterr->second;
				}
				delete postponedNodes;
			}

			delete subgraph;
			delete newCandidate;
		}
		//I recieve a subgraph for isfrequency (approximate)
		else if(str_message[0]=='a')
		{
			tr1::unordered_map<string, void* > allMPs;
			//load graph data
			std::stringstream sstmGF;
			int candidateID;
			char subgraphTemp[500];
			sscanf (str_message,"a,%d,%[^\t]",&candidateID,subgraphTemp);
			if(Settings::debugMSG)
				printf("Worker %d: Recieved IsFreq [approximate] for subgraph ID: %d. subgraph str=%s\n", worker->machineID, candidateID, subgraphTemp);
			sstmGF << subgraphTemp;
			string subgraphStr = sstmGF.str();
			GraphX* subgraph = new GraphX(1, false);//create a graph with id=its partition id
			subgraph->loadFromString(subgraphStr, allMPs);
			clearAllMP(allMPs);
			//printf("I am a worker with subgraph\n'%s'\n", subgraphStr.c_str());
			Pattern* newCandidate = new Pattern(subgraph, false);
			bool isFreq;
			unsigned long startTime = getmsofday();

			if(Settings::debugMSG)
				printf("Worker %d: Starts checking approximate frequency\n", (worker->machineID*worker->nThreads+id));
			int invalidCol;
			int predictedValids;
			bool exact;//this is set to true if the approximate function can return exact result (by using one of the first optimizations)

			unsigned long numVisitedNodes;
			unsigned long numIterations;
			unsigned long predictedTime;
			int freq = worker->getApproxFreq(newCandidate, invalidCol, predictedValids, exact, numVisitedNodes, numIterations, predictedTime);

			if(freq>=worker->frequency)
				isFreq = true;
			else
				isFreq = false;

			if(Settings::debugMSG)
				printf("Worker %d: Finished checking approximate frequency. Frequency: %d. Frequent?: %d, and invalid column is: %d with predicted number of valid nodes = %d. Num of visited nodes = %Lu. Num of iterations = %Lu\n", threadID, freq, isFreq, invalidCol, predictedValids, numVisitedNodes, numIterations);

			unsigned long elapsed = getmsofday() - startTime;

			if(Settings::debugMSG)
				printf("Time taken is: %Lu. Predicted total time is: %lu\n", elapsed);

			if(isFreq)
			{
				sprintf(str_message, "at%d,%d,%Lu,%Lu,%d,%Lu,%Lu", candidateID, freq, elapsed, predictedTime,exact,numVisitedNodes, numIterations);
				int cnt=strlen(str_message)+1;
				MPI_Send(str_message, cnt, MPI_CHAR, 0, id, MPI_COMM_WORLD);
			}
			else
			{
				sprintf(str_message, "af%d,%d,%Lu,%Lu,%d,%d,%d,%Lu,%Lu", candidateID, freq, elapsed, predictedTime, invalidCol, predictedValids,exact, numVisitedNodes, numIterations);
				int cnt=strlen(str_message)+1;
				MPI_Send(str_message, cnt, MPI_CHAR, 0, id, MPI_COMM_WORLD);
			}

			delete subgraph;
			delete newCandidate;
		}
		else if(str_message[0]=='t')
		{
		}

		delete[] str_message;
	}
}

bool Worker::isFrequent(Pattern* candidate, int subTaskNum, int subTaskMax, int* mniTable, tr1::unordered_map<int, tr1::unordered_set<int>* >* postponedNodes, unsigned long maxIters)
{
	int freq = GraMiCounter::isFrequent_adv(graph, candidate, frequency, subTaskNum, subTaskMax, mniTable, postponedNodes, maxIters);

	if(freq>=frequency)
		return true;

	return false;
}

int Worker::getApproxFreq(Pattern* candidate, int& invalidCol, int& validNodesPredicted, bool& exact, unsigned long& numVisiteNodes, unsigned long& numIterations, unsigned long& predictedTime)
{
	vector<unsigned long> listOfNumOfIters;
	int freq = GraMiCounter::isFrequent_approx(graph, candidate, frequency, invalidCol, validNodesPredicted, true, exact, numVisiteNodes, numIterations, listOfNumOfIters, predictedTime);

	//http://www.mathwords.com/o/outlier.htm
	//get the outlier threshold to set it to the numIterations
	if(listOfNumOfIters.size()==0)
	{
		numIterations = Settings::postponeNodesAfterIterations;
		if(Settings::debugMSG)
			cout<<"No node is sampled"<<endl;
	}
	else
	{
		unsigned long min = listOfNumOfIters[0];
		unsigned long max = listOfNumOfIters[listOfNumOfIters.size()-1];
		unsigned long firstQuartile = listOfNumOfIters[listOfNumOfIters.size()*0.25];
		unsigned long median = listOfNumOfIters[listOfNumOfIters.size()*0.5];
		unsigned long thirdQuartile = listOfNumOfIters[listOfNumOfIters.size()*0.75];
		unsigned long interquartileRange = thirdQuartile - firstQuartile;
		unsigned long outlier = thirdQuartile + 1.5*interquartileRange;
		numIterations = outlier*10;
	}

	return freq;
}
