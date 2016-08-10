/*
 * CLMap.h
 *
 *  Created on: Dec 16, 2015
 *      Author: abdelhe
 */

/**
 * a helper structure for subgraph isomorphism function
 */
#ifndef CLMAP_H_
#define CLMAP_H_

#include <map>
#include<list>
#include <string.h>
#include "Pattern.h"

using namespace std;

class CLMap_Iterator
{
public:
	string key;
	Pattern* pattern;
	map<string, Pattern*>::iterator mapPIter;
	map<string, list<Pattern*>* >::iterator mapVIter;
	list<Pattern* >::iterator vectIter;
	CLMap_Iterator getCopy();
};

class CLMap
{
private:
	map<string, Pattern*> patterns_sig;
	map<string, list<Pattern* >* > patterns_nosig;
	Pattern* exists(Pattern*, list<Pattern*>* );
	bool remove(Pattern*, list<Pattern*>* );
	bool theSame(Pattern*, Pattern* );
	int size;

public:
	CLMap();
	bool addPattern(Pattern* pattern);
	void removePattern(Pattern* pattern);
	void addAll(CLMap* clmap);
	Pattern* getPattern(Pattern* pattern);
	int getSize() { return size; }
	void print(int counter = 1);
	CLMap_Iterator getFirstElement();
	void advanceIterator(CLMap_Iterator& );
	void clear();
	void deleteObjects();
};



#endif /* CLMAP_H_ */
