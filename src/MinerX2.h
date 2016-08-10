/**
 * Copyright 2016 Ehab Abdelhamid, Ibrahim Abdelaziz

MinerX2.h

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

#ifndef MINERX2_H_
#define MINERX2_H_

#include "MinerX.h"

class MinerX2: public MinerX
{
private:
	int divideF2Counter = 0;
	int nThreads;
	int numMachines;

protected:
	Pattern* getPredictedCandidate(Pattern* );
	void sendACandidate(string key, Pattern* candidate, map<int, Pattern*>& currentlyChecking, int destination, bool evalType);
	double timeThreshold;
	bool checkAgainstOtherPatterns(Pattern* candidate);
	Pattern* selectCandidate(string& patternKey, int& evalType);
	long function1(long, long, double, vector<CLMap* >** );
	double overhead = 0.000484;//for twitter = 0.00011, for patents = 0.000484
	int getNumberOfExpectedSubTasks();
	void divideFFun3(vector<CLMap* >** );
	void divideFFun4(vector<CLMap* >**);
	int getNumAvailWorkers();

public:
	MinerX2(int numMachines, int nThreads): MinerX(numMachines, nThreads) {this->nThreads = nThreads;this->numMachines = numMachines;}
	void startMining();
	void setExpectedPatterns(vector<CLMap* >* exFreqPatt, vector<CLMap* >* exInfreqPatt);
	void setSupport(int numWorkers, int support);
	void initWithoutApprox(string fileName, int graphType, int support, int numMachine, int nThreads);
	void patternFoundInFrequent(Pattern* );

	static long genPrimaryGraphTime;
	static long patternFreqTime;
	static long patternInfreqTime;
	static long sendCandidsTime;
	static long predictedPatternTime;
};


#endif /* MINERX2_H_ */
