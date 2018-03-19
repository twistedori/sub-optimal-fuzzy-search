#include "Form.h"

Form::Form(unsigned char formTypesSize)
{
	positions = new vector<vector<unsigned int>>;
	positions->resize(formTypesSize);
}

Form::Form(unsigned char formTypesSize, unsigned int proteinIndex)
{
	positions = new vector<vector<unsigned int>>;
	positions->resize(formTypesSize);
}

Form::Form(unsigned char formTypesSize, unsigned char formType, unsigned int position)
{
	positions = new vector<vector<unsigned int>>;
	positions->resize(formTypesSize);
	addPosition(formType, position);
}

Form::~Form()
{
	delete positions;
}

vector<vector<unsigned int>>* Form::getPositions()
{
	return positions;
}

void Form::addPosition(unsigned char formType, unsigned int position)
{
	(*positions)[formType].push_back(position);
}