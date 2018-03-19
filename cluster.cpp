#include "Cluster.h"
#include <vector>

ClusterClass::ClusterClass(unsigned int offset)
{
	seed = offset;
}


unsigned int ClusterClass::getSeed()
{
	return seed;
}

vector<unsigned int> ClusterClass::getPositions()
{
	vector<unsigned int> tmp;
	tmp.reserve(offsets0.size() + offsets1.size() + offsets2.size() + 1);
	tmp.insert(tmp.end(), seed);
	tmp.insert(tmp.end(), offsets0.begin(), offsets0.end());
	tmp.insert(tmp.end(), offsets1.begin(), offsets1.end());
	tmp.insert(tmp.end(), offsets2.begin(), offsets2.end());
	return tmp;
}

vector<unsigned int> ClusterClass::getPositionsLvl0()
{
	vector<unsigned int> tmp;
	tmp.reserve(offsets0.size() + 1);
	tmp.insert(tmp.end(), seed);
	tmp.insert(tmp.end(), offsets0.begin(), offsets0.end());
	return tmp;
}

vector<unsigned int> ClusterClass::getPositionsLvl1()
{
	vector<unsigned int> tmp;
	tmp.reserve(offsets1.size() + 1);
	tmp.insert(tmp.end(), seed);
	tmp.insert(tmp.end(), offsets1.begin(), offsets1.end());
	return tmp;
}

vector<unsigned int> ClusterClass::getPositionsLvl2()
{
	vector<unsigned int> tmp;
	tmp.reserve(offsets2.size() + 1);
	tmp.insert(tmp.end(), seed);
	tmp.insert(tmp.end(), offsets2.begin(), offsets2.end());
	return tmp;
}

void ClusterClass::addPosibilitie(unsigned int posb)
{
	posibleConnections[getSeed()][posb] = true;
	possibilities.push_back(posb);
}

void ClusterClass::addConnection(unsigned int conc)
{
	connections.push_back(conc);
}

vector<unsigned int> ClusterClass::getPosibilities()
{
	vector<unsigned int> tmp;
	tmp.reserve(possibilities.size() + 1);
	tmp.insert(tmp.end(), possibilities.begin(), possibilities.end());
	return tmp;
}

vector<unsigned int> ClusterClass::getConnections()
{
	vector<unsigned int> tmp;
	tmp.reserve(connections.size() + 1);
	tmp.insert(tmp.end(), connections.begin(), connections.end());
	return tmp;
}


void ClusterClass::addPosition(unsigned int position, unsigned int dist)
{
	if (dist <= 2) {
		switch (dist) {
		case 0:
			offsets0.push_back(position);
			break;
		case 1:
			offsets1.push_back(position);
			break;
		case 2:
			offsets2.push_back(position);
			break;
		default:
			printf_s(" ERROR! \n\n there is less than 90% accurcy");
			getchar();
			break;

		}


	}


}


