/**
 * Copyright 2016 Ehab Abdelhamid, Ibrahim Abdelaziz

QueueSimulator.cpp

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
 * This class is for the event-driven tasks simulation
 */

#include <math.h>
#include<iostream>
#include <stdlib.h>
#include <cstring>
#include "QueueSimulator.h"
#include "Settings.h"

using namespace std;

/**
 * Strat the simulation given the list of tasks and the number of workers
 */
long QueueSimulator::simulate(long* times, int numTimes, int nWorkers)
{
	long totalTime = 0;
	long maxTime = 0;
	for(int i=0;i<numTimes;i++)
	{
		totalTime+=times[i];
		if(maxTime<times[i])
			maxTime = times[i];
		if(Settings::debugMSG)
			cout<<"Total Time = "<<totalTime<<","<<times[i]<<endl;
	}

	long startThreshold;
	startThreshold = maxTime;//use the maximum task time as initial value

	long timePerWorker = totalTime/nWorkers;

	long threshold = startThreshold;

	//iterate until the imbalance value is small enough
	long eval;
	while(true)
	{
		eval = evaluate(times, numTimes, nWorkers, threshold, timePerWorker);

		if(Settings::debugMSG)
			cout<<eval<<" , "<<Settings::minImbalance<<endl;

		//the break condition
		if(eval<Settings::minImbalance)
			break;

		threshold = threshold - (threshold*0.1);
	}

	return threshold;
}

/**
 * evaulate a specific setup of how tasks are divided
 */
double QueueSimulator::evaluate(long* times, int numTimes, int nWorkers, long threshold, long avg)
{
	//initialize workers times
	long* workerTime = new long[nWorkers];
	for(int i=0;i<nWorkers;i++)
		workerTime[i] = 0;

	//go over each time, distribute it over the workerTimes
	for(int i=0;i<numTimes;i++)
	{
		int numSubtasks = ceil(((double)times[i])/threshold);
		if(numSubtasks>nWorkers)
			numSubtasks = nWorkers;
		for(int j=0;j<numSubtasks;j++)
		{
			long valueToAdd = times[i]/numSubtasks;

			//select the worker time with lowest value
			int index = 0;
			long min = workerTime[0];
			for(int k = 1;k<nWorkers;k++)
			{
				if(min>workerTime[k])
				{
					index = k;
					min = workerTime[k];
				}
			}

			workerTime[index]+=valueToAdd;
		}
	}

	//get the difference between the maximum time and the ideal average time per worker
	int index = 0;
	long max = workerTime[0];
	for(int i=0;i<nWorkers;i++)
	{
		if(Settings::debugMSG)
			cout<<"Worker"<<i<<": "<<workerTime[i]<<",";
		if(max<workerTime[i])
		{
			index = i;
			max = workerTime[i];
		}
	}
	if(Settings::debugMSG)
		cout<<endl;

	delete[] workerTime;

	double r = ((((double)max)/avg)-1)*100;
	if(Settings::debugMSG)
		cout<<"Returning: "<<r<<endl;

	//return max-avg;
	return r;
}
