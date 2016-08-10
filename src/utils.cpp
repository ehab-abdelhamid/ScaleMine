/**
 * Copyright 2016 Ehab Abdelhamid, Ibrahim Abdelaziz

utils.cpp

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

#include <fstream>
#include<sstream>
#include <string.h>
#include <algorithm>
#include <unistd.h>
#include"utils.h"

/**
 * A function to copy map content into a vector
 */
void MapToVec(map<string, CL_Partition* >& m, vector<CL_Partition* >& v) {
    for( map<string, CL_Partition* >::const_iterator it = m.begin(); it != m.end(); ++it ) {
    	v.push_back( it->second );
    }
}

string intToString(int a)
{
	stringstream sstmGF;
	sstmGF << a;
	return sstmGF.str();
}

string doubleToString(double a)
{
	stringstream sstmGF;
	sstmGF << a;
	return sstmGF.str();
}

int getPos(char* str, char c)
{
	for(int i=0;i<strlen(str);i++)
	{
		if(str[i]==c)
			return i;
	}
	return -1;
}

//get time in milli seconds
long long getmsofday()
{
   struct timeval tv;
   struct timezone tz;
   gettimeofday(&tv, &tz);
   return (long long)tv.tv_sec*1000 + tv.tv_usec/1000;
}

char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}

int getNumElems(vector<map<string, void*>* >* a)
{
	int counter = 0;
	for(int i=0;i<a->size();i++)
	{
		counter+=a->at(i)->size();
	}
	return counter;
}

void vect_map_destruct(vector<map<string, Pattern*>* > vm)
{
	for(vector<map<string, Pattern*>* >::iterator iter1 = vm.begin();iter1!=vm.end();++iter1)
	{
		for(map<string, Pattern*>::iterator iter2 = (*iter1)->begin(); iter2!=(*iter1)->end();iter2++)
		{
			delete (iter2->second);

		}

		(*iter1)->clear();
		delete (*iter1);

	}
}

void vect_map_destruct(vector<CLMap* > vm)
{
	for(vector<CLMap* >::iterator iter1 = vm.begin();iter1!=vm.end();++iter1)
	{
		(*iter1)->deleteObjects();
		(*iter1)->clear();
		delete (*iter1);
	}
}

/**
 * get current memory usage
 */
void process_mem_usage(double& vm_usage, double& resident_set)
{
   using std::ios_base;
   using std::ifstream;
   using std::string;

   vm_usage     = 0.0;
   resident_set = 0.0;

   // 'file' stat seems to give the most reliable results
   //
   ifstream stat_stream("/proc/self/stat",ios_base::in);

   // dummy vars for leading entries in stat that we don't care about
   //
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;

   // the two fields we want
   //
   unsigned long vsize;
   long rss;

   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
               >> utime >> stime >> cutime >> cstime >> priority >> nice
               >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

   stat_stream.close();

   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
   vm_usage     = vsize / 1024.0;
   resident_set = rss * page_size_kb;
}

long getTotalSystemMemory()
{
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size;
}
