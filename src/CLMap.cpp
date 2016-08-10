/*
 * CLMap.cpp
 *
 *  Created on: Dec 16, 2015
 *      Author: abdelhe
 */

/**
 * a helper structure for subgraph isomorphism function
 */
#include<vector>
#include <list>
#include <iostream>
#include <stdlib.h>

#include "CLMap.h"
#include "Settings.h"

CLMap_Iterator CLMap_Iterator::getCopy()
{
	CLMap_Iterator iter;

	iter.pattern = pattern;
	iter.key = key;
	iter.mapPIter = mapPIter;
	iter.mapVIter = mapVIter;
	iter.vectIter = vectIter;

	return iter;
}

CLMap::CLMap()
{
	size = 0;
}

/**
 * add the pattern to the list of patterns
 * When the signature creation takes forever, use the list and compare using subgraph isomorphism
 */
bool CLMap::addPattern(Pattern* pattern)
{
	//ignore this pattern if its size exceeds the given limit
	if((Settings::maxSubgraphSize>-1 && pattern->getGraph()->getNumOfEdges()>Settings::maxSubgraphSize) ||
			(Settings::maxNumNodes>-1 && pattern->getGraph()->getNumOfNodes()>Settings::maxNumNodes))
		return false;

	string sig = pattern->getGraph()->getCanonicalLabel();
	if(sig.at(0)=='X')
	{
		list<Pattern* >* v;
		map<string,list<Pattern* >* >::iterator it = patterns_nosig.find(sig);
		if(it!=patterns_nosig.end())
		{
			v = it->second;
			if(exists(pattern, v))
				return false;
		}
		else
		{
			v = new list<Pattern* >();
			patterns_nosig.insert(std::pair<string, list<Pattern* >* >(sig, v));
		}

		v->push_back(pattern);
	}
	else
	{
		if(patterns_sig.find(sig)!=patterns_sig.end())
			return false;

		patterns_sig.insert(std::pair<string,Pattern* >(sig,pattern));
	}

	size++;
	return true;
}

void CLMap::removePattern(Pattern* pattern)
{
	string sig = pattern->getGraph()->getCanonicalLabel();
	if(sig.at(0)=='X')
	{
		list<Pattern* >* v;
		map<string,list<Pattern* >* >::iterator it = patterns_nosig.find(sig);
		if(it!=patterns_nosig.end())
		{
			v = it->second;
			if(remove(pattern, v))
			{
				size--;
				if(v->size()==0)
				{
					patterns_nosig.erase(sig);
				}
			}
		}
	}
	else
	{
		int a = patterns_sig.erase(sig);
		size-=a;
	}
}

/**
 * search for the given pattern, if not exist then the returned value = 0
 */
Pattern* CLMap::getPattern(Pattern* pattern)
{
	string sig = pattern->getGraph()->getCanonicalLabel();

	if(sig.at(0)=='X')
	{
		std::map<string,list<Pattern* >* >::iterator it = patterns_nosig.find(sig);
		if(it==patterns_nosig.end())
		{
			return 0;
		}

		return exists(pattern, it->second);
	}
	else
	{
		std::map<string,Pattern* >::iterator it = patterns_sig.find(sig);
		if(it==patterns_sig.end())
		{
			return 0;
		}
		else
		{
			return it->second;
		}
	}
}

/**
 * search for the pattern in the given list, if found return a pointer to its similar pattern
 */
Pattern* CLMap::exists(Pattern* pattern, list<Pattern*>* v)
{
	if(Settings::debugMSG)
	{
		cout<<"CLMap::exists starts..."<<endl;
	}

	for(list<Pattern*>::iterator iter = v->begin();iter!=v->end();iter++)
	{
		Pattern* temp = *iter;
		if(theSame(pattern, temp))
		{
			if(Settings::debugMSG)
			{
				cout<<"CLMap::exists finished, pattern: "<<temp->getID()<<" found."<<endl;
			}

			return temp;
		}
	}

	if(Settings::debugMSG)
	{
		cout<<"CLMap::exists finished, no pattern found."<<endl;
	}
	return 0;
}

/*
 * remove a pattern if found, otherwise no effect and return false
 * return true only if the pattern is removed
 */
bool CLMap::remove(Pattern* pattern, list<Pattern*>* v)
{
	for(list<Pattern*>::iterator iter = v->begin();iter!=v->end();iter++)
	{
		Pattern* temp = *iter;
		if(theSame(pattern, temp))
		{
			v->erase(iter);
			return true;
		}
	}
	return false;
}

//check if they are the same by performing subgraph isomorphism
bool CLMap::theSame(Pattern* p1, Pattern* p2)
{
	return p1->getGraph()->isTheSameWith(p2->getGraph());
}

void CLMap::addAll(CLMap* clmap)
{
	for(map<string, Pattern*>::iterator iter = clmap->patterns_sig.begin();iter!=clmap->patterns_sig.end();iter++)
	{
		patterns_sig.insert(std::pair<string,Pattern* >(iter->first, iter->second));
		size++;
	}

	for(map<string, list<Pattern*>* >::iterator iter = clmap->patterns_nosig.begin();iter!=clmap->patterns_nosig.end();iter++)
	{
		list<Pattern* >* v = new list<Pattern* >();
		v->insert(v->begin(), iter->second->begin(), iter->second->end());
		patterns_nosig.insert(std::pair<string,list<Pattern*>* >(iter->first, v));
		size+=v->size();
	}
}

CLMap_Iterator CLMap::getFirstElement()
{
	CLMap_Iterator iter;

	if(patterns_sig.size()>0)
	{
		iter.key = patterns_sig.begin()->first;
		iter.pattern = patterns_sig.begin()->second;
		iter.mapPIter = patterns_sig.begin();
		return iter;
	}
	else if(patterns_nosig.size()>0 && patterns_nosig.begin()->second->size()>0)
	{
		iter.mapPIter = patterns_sig.end();

		iter.key = patterns_nosig.begin()->first;
		iter.mapVIter = patterns_nosig.begin();
		iter.pattern = *(patterns_nosig.begin()->second->begin());
		iter.vectIter = patterns_nosig.begin()->second->begin();
		return iter;
	}

	iter.pattern = 0;
	if(size>0)
	{
		cout<<"Error56333: size="<<size<<" though I cannot find a first pattern."<<endl;
		exit(0);
	}

	return iter;
}

void CLMap::advanceIterator(CLMap_Iterator& currentIter)
{
	bool b = false;

	//iterating over the map of sigs
	if(currentIter.mapPIter != patterns_sig.end())
	{
		currentIter.mapPIter++;
		if(currentIter.mapPIter != patterns_sig.end())
		{
			currentIter.key = currentIter.mapPIter->first;
			currentIter.pattern = currentIter.mapPIter->second;
			return;
		}

		currentIter.mapVIter = patterns_nosig.begin();
		b = true;
	}

	if(currentIter.mapVIter==patterns_nosig.end())
	{
		currentIter.pattern = 0;
		return;
	}

	if(b)
	{
		currentIter.vectIter = currentIter.mapVIter->second->begin();
	}
	else
	{
		currentIter.vectIter++;
	}

	if(currentIter.vectIter==currentIter.mapVIter->second->end())
	{
		currentIter.mapVIter++;
		if(currentIter.mapVIter==patterns_nosig.end())
		{
			currentIter.pattern = 0;
			return;
		}
		currentIter.vectIter = currentIter.mapVIter->second->begin();
	}

	currentIter.key = currentIter.mapVIter->first;
	currentIter.pattern = *(currentIter.vectIter);
}

void CLMap::print(int counter)
{
	CLMap_Iterator iter = this->getFirstElement();

	while(true)
	{
		if(iter.pattern==0)
			break;

		cout<<counter<<":"<<endl;
		cout<<*(iter.pattern->getGraph())<<endl;
		cout<<"-----------"<<endl;
		this->advanceIterator(iter);
		counter++;
	}
}

void CLMap::clear()
{
	patterns_sig.clear();
	for(map<string, list<Pattern* >* >::iterator iter = patterns_nosig.begin();iter!=patterns_nosig.end();iter++)
	{
		delete iter->second;
	}
	patterns_nosig.clear();
	//size = 0;
}

void CLMap::deleteObjects()
{
	CLMap_Iterator iter = this->getFirstElement();
	while(iter.pattern!=0)
	{
		Pattern* pattern = iter.pattern;
		delete pattern;
		this->advanceIterator(iter);
	}
}
