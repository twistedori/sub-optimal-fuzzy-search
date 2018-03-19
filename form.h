#pragma once


#ifndef FORM_H_
#define FORM_H_

#include <vector>
#include <ctime>
#include <map>


using namespace std;

const int NumberOfAminoAcids = 20;

const char AAcode[] = { 'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V' };





struct FormType
{
	string formType;
	vector<int> distances;
};

struct ProteinPosition
{
	unsigned int startIndex;
	unsigned int endIndex;
};

class Form
{
private:
	//Every form type has a vector of positions
	vector<vector<unsigned int>> *positions;
public:
	Form(unsigned char formTypesSize);
	Form(unsigned char formTypesSize, unsigned int proteinIndex);
	Form(unsigned char formTypesSize, unsigned char formType, unsigned int position);
	~Form();
	vector<vector<unsigned int>>* getPositions();
	void addPosition(unsigned char formType, unsigned int position);
};

#endif /* FORM_H_ */
