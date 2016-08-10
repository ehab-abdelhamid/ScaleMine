/**
 * Copyright 2016 Ehab Abdelhamid, Ibrahim Abdelaziz

QueueSimulator.h

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

#ifndef QUEUESIMULATOR_H_
#define QUEUESIMULATOR_H_


class QueueSimulator
{
public:
	static long simulate(long* times, int numTimes, int nWorkers);
	static double evaluate(long* times, int numTimes, int nWorkers, long threshold, long avgTimePerWorker);
};


#endif /* QUEUESIMULATOR_H_ */
