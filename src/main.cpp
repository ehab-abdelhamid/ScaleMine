/**
 * Copyright 2016 Ehab Abdelhamid, Ibrahim Abdelaziz

main.cpp

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

#include<iostream>
#include<mpi.h>
#include <cstdlib>
#include"Settings.h"
#include"Worker.h"
#include"MinerAdv.h"
#include"MinerX2.h"
#include"utils.h"
#include "GraMiCounter.h"
#include "CLMap.h"
#include "QueueSimulator.h"
#include "Pattern.h"

#include<vector>
#include<map>

using namespace std;

int CL_Partition::st_counter;
int NodeWithInfo::st_counter;
int PartID_label::st_counter;
int Pattern::st_counter;
int Set_Iterator::st_SI;
bool Settings::useSearchSpacePrediction;
bool Settings::divideBigTasks;
bool Settings::predictOnTheGo;
long Settings::graphLoadingTime;
long Settings::phase2Time;
bool Settings::debugMSG;
bool Settings::showMemUse;
int Settings::divideByFunc;
int Settings::maxSubgraphSize;
int Settings::maxNumNodes;
int Settings::fixedNumSubtasks;
bool Settings::showNumCandids;
double Settings::minImbalance;
bool Settings::usePredictedInvColumn;
bool Settings::smartBreak;

int Settings::minNumberOfSamples;
int Settings::maxNumberOfSamples;

long MinerAdv::maxNumCandids;
long Settings::postponeNodesAfterIterations;
bool Settings::postponeExpensiveNodes;
int GraMiCounter::numPostponedNodes = 0;

//time measures for the matser work
long MinerX2::genPrimaryGraphTime = 0;
long MinerX2::patternFreqTime = 0;
long MinerX2::patternInfreqTime = 0;
long MinerX2::sendCandidsTime = 0;
long MinerX2::predictedPatternTime = 0;

int clientFirstTask = 0;
#define BUFFERSIZE 500


int main( int argc, char *argv[] )
{
	//set application parameters
	GraMiCounter::numSideEffectNodes = 0;
	MinerAdv::maxNumCandids = -1;
	Settings::fixedNumSubtasks = -1;
	Settings::useSearchSpacePrediction = true;
	Settings::divideBigTasks = true;
	Settings::predictOnTheGo = true;
	Settings::debugMSG = false;
	Settings::showMemUse = false;
	Settings::postponeExpensiveNodes = true;
	Settings::postponeNodesAfterIterations = 10000000;
	Settings::divideByFunc = 4;
	Settings::showNumCandids = false;
	Settings::minImbalance = 15;
	Settings::usePredictedInvColumn = true;
	Settings::smartBreak = true;

	//The miner routine
	Set_Iterator::st_SI = 0;
	long long start = getmsofday();

	//application default parameters
	string fileName = "../Datasets/test.lg";
	int graphType = 0;//undirected graph
	int support = 9460;

	int nThreads = 1;

	//load user-given parameter values
	//get the number of threads
	char * argThreads = getCmdOption(argv, argv + argc, "-threads");
	if(argThreads)
	{
		nThreads = atoi(argThreads);
	}

	//load graph file
	char * argfilename = getCmdOption(argv, argv + argc, "-file");
	if(argfilename)
	{
		fileName = string(argfilename);
	}

	//get user-given support threshold
	char * argSupport = getCmdOption(argv, argv + argc, "-freq");
	if(argSupport)
	{
		support = atoi(argSupport);
	}

	//parameter to set the maximum subgraph size (in terms of the number of edges)
	char * argMaxSize = getCmdOption(argv, argv + argc, "-maxSize");
	if(argMaxSize)
	{
		Settings::maxSubgraphSize = atoi(argMaxSize);
	}
	else
		Settings::maxSubgraphSize = -1;

	//parameter to set the maximum subgraph size (in terms of the number of nodes)
	char * argMaxNodes = getCmdOption(argv, argv + argc, "-maxNodes");
	if(argMaxNodes)
	{
		Settings::maxNumNodes = atoi(argMaxNodes);
	}
	else
		Settings::maxNumNodes = -1;

	//parameter to set the maximum subgraph size (in terms of the number of edges)
	char * argFixedSubtasks = getCmdOption(argv, argv + argc, "-fixedSubtasks");
	if(argFixedSubtasks)
	{
		Settings::fixedNumSubtasks = atoi(argFixedSubtasks);
	}
	else
		Settings::fixedNumSubtasks = -1;

	//parameter to set the maximum number of samples
	char * argMaxSamples = getCmdOption(argv, argv + argc, "-maxSamples");
	if(argMaxSamples)
	{
		Settings::maxNumberOfSamples = atoi(argMaxSamples);
	}
	else
		Settings::maxNumberOfSamples = 300;

	//parameter to set the maximum number of samples
	char * argMinSamples = getCmdOption(argv, argv + argc, "-minSamples");
	if(argMinSamples)
	{
		Settings::minNumberOfSamples = atoi(argMinSamples);
	}
	else
		Settings::minNumberOfSamples = 60;

	//initialize MPI
	int i;
	int rank;
	int numWorkers;
	MPI_Status    status;
	srand (time(NULL));

	int provided = MPI::Init_thread(argc, argv, MPI_THREAD_MULTIPLE);
	if(provided != MPI_THREAD_MULTIPLE){
		cout<<"Provided = "<<provided<<endl;
		cout<<"ERROR: Cannot initialize MPI with the desire level of support";
	}
	rank = MPI::COMM_WORLD.Get_rank();
	MPI_Comm_size(MPI_COMM_WORLD, &numWorkers);

	//start the master process
	int cnt=0;
	if(rank==0)
	{
		long long elapsedApprox = 0;
		GraMiCounter::numByPassedNodes = 0;

		//start approximate mining
		MinerAdv* miner = new MinerAdv();

		int approxSupport = support - (support*Settings::errorMargin);

		cout<<"Inexact search space building starts .................................."<<endl;

		long long startApprox = getmsofday();

		miner->startMining(fileName, graphType, approxSupport, numWorkers, nThreads);

		long long endApprox = getmsofday();
		elapsedApprox = endApprox - startApprox;

		miner->printTotalExpectedTime();

		//do exact mining using the approximated search space
		MinerX* exactMiner;

		cout<<"Exact mining starts .................................."<<endl;
		exactMiner = new MinerX2(numWorkers, nThreads);

		exactMiner->setInputGraph(miner->getInputGraph());
		miner->setInputGraph(0);

		exactMiner->setSupport(numWorkers, support);

		if(miner->numOfVisitedNodes>0)
			exactMiner->setAvgIterPerNode(miner->numIterations/miner->numOfVisitedNodes);
		else
			exactMiner->setAvgIterPerNode(0);

		//borrow values from the approximate miner
		exactMiner->setExpectedPatterns(miner->getFrequentPatterns(), miner->getInfrequentPatterns());

		exactMiner->setFrequentEdges(miner->getFrequentEdges());

		//start mining

		exactMiner->startMining();

		long long end = getmsofday();

		long long elapsed = end - start;

		cout<<"Graph loading [master] took "<<(Settings::graphLoadingTime/1000)<<" sec and "<<(Settings::graphLoadingTime%1000)<<" ms"<<endl;

		cout<<"Mining took "<<(elapsed/1000)<<" sec and "<<(elapsed%1000)<<" ms"<<endl;

		cout<<"Finished!";

		MPI::Finalize();

		delete miner;
		delete exactMiner;
	}
	else	//start workers
	{
		Worker worker;
		worker.start(rank, nThreads);
	}

	return 0;
}
