/**
 * Copyright 2016 Ehab Abdelhamid, Ibrahim Abdelaziz

utils.h

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

#ifndef UTILS_H_
#define UTILS_H_

#include<map>
#include<vector>
#include<string>
#include<sys/time.h>
#include "Pattern.h"
#include"CanonicalLabel.h"
#include "CLMap.h"

using namespace std;

void MapToVec(map<string, CL_Partition* >& , vector<CL_Partition* >& );
string intToString(int a);
string doubleToString(double a);
int getPos(char* str, char c);
long long getmsofday();
char* getCmdOption(char ** begin, char ** end, const std::string & option);
bool cmdOptionExists(char** begin, char** end, const std::string& option);
int getNumElems(vector<map<string, void*>* >* );
void vect_map_destruct(vector<map<string, Pattern*>* > vm);
void vect_map_destruct(vector<CLMap* > vm);
long getTotalSystemMemory();
void process_mem_usage(double& vm_usage, double& resident_set);

#endif /* UTILS_H_ */
