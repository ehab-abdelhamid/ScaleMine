/**
 * Copyright 2016 Ehab Abdelhamid, Ibrahim Abdelaziz

Settings.h

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

#ifndef SETTINGS_H_
#define SETTINGS_H_


class Settings
{
public:
	static bool useSearchSpacePrediction;
	static int maxNumberOfSamples;
	static int minNumberOfSamples;
	static const int samplesBucket = 20;
	//below are used for the approximate results mode .. the goal is to get more accurate results
	static const int minNumberOfSamples_accur = 120;
	static int samplesBucket_accur;

	static const bool fullCount = true;
	static const double errorMargin = 0;
	static bool divideBigTasks;
	static bool predictOnTheGo;//true means I can estimate new patterns that were never estimated before
	static const double validsPerTaskMargin = 0.8;
	static long graphLoadingTime;
	static long phase2Time;
	static bool debugMSG;
	static bool showMemUse;
	static int divideByFunc; //0- use the avergae
							 //1- use a value to minimize the following: threshold+overall_execution_time*#subtasks*overhead_constant
	static bool postponeExpensiveNodes;
	static long postponeNodesAfterIterations;
	static int maxSubgraphSize;
	static int maxNumNodes;
	static bool usePredictedInvColumn;
	static bool smartBreak;

	static int fixedNumSubtasks;
	static bool stopAfterApproximation;
	static bool showNumCandids;//mainly used for showing how the number of candidates increase as the time goes
	static double minImbalance;
};

#endif /* SETTINGS_H_ */
