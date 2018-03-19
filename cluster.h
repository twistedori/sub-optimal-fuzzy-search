#pragma once

#ifndef CLUSTER_H_
#define CLUSTER_H_

#include <map>
#include "Form.h"


extern bool** posibleConnections;

class ClusterClass
{
private:
	unsigned int seed;
	vector<unsigned int> connections;
	vector<unsigned int> possibilities;
	vector<unsigned int> offsets0;
	vector<unsigned int> offsets1;
	vector<unsigned int> offsets2;
public:
	ClusterClass(unsigned int location);

	unsigned int getSeed();
	vector<unsigned int> getPositions();
	vector<unsigned int> getPositionsLvl0();
	vector<unsigned int> getPositionsLvl1();
	vector<unsigned int> getPositionsLvl2();
	void addPosition(unsigned int position, unsigned int dist);
	void addPosibilitie(unsigned int posb);
	void addConnection(unsigned int conc);
	vector<unsigned int> getPosibilities();
	vector<unsigned int> getConnections();

};







#endif /* CLUSTER_H_ */
