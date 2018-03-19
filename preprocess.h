#pragma once


#ifndef PREPROCESS_H_
#define PREPROCESS_H_


#include <string>
#include <unordered_map>
#include "Form.h"
#include "Cluster.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include <set>

using namespace std;

class PreprocessClass
{
private:
	int AATabQuick[255];
public:
	void PreProcess(std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>> &forms, vector<FormType> &formTypes, int form_size, int mismatches, string &db, multimap<int, ProteinPosition> &proteinsPositions, int seqSize);

	void PreProcess42(std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>> &forms, std::map<int, std::map<int, std::map<int, std::map<int, vector<Form*>>>>> &forms2, vector<FormType> &formTypes, int form_size, int mismatches, string &db, multimap<int, ProteinPosition> &proteinsPositions, int seqSize, vector<ClusterClass> &clusters);

	void PreProcess5(std::map<int, std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>>> &forms, vector<FormType> &formTypes, int form_size, int mismatches, string &db, multimap<int, ProteinPosition> &proteinsPositions, int seqSize);

	void PreProcess53(std::map<int, std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>>> &forms, vector<FormType> &formTypes, int form_size, int mismatches, string &db, multimap<int, ProteinPosition> &proteinsPositions, int seqSize, vector<ClusterClass> &clusters);

	void PreProcess6(std::map<int, std::map<int, std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>>>> &forms, vector<FormType> &formTypes, int form_size, int mismatches, string &db, multimap<int, ProteinPosition> &proteinsPositions, int seqSize);

	bool CheckIfPreProcessed(int form_size, int mismatches, int seqSize);

	void GetFormsHash(std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>> &forms, vector<FormType> &formTypes, string &db, multimap<int, ProteinPosition> &proteinsPositions, int seqSize, int form_size);

	void GetFormsHash42(std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>> &forms, std::map<int, std::map<int, std::map<int, std::map<int, vector<Form*>>>>> &forms2, vector<FormType> &formTypes, string &db, multimap<int, ProteinPosition> &proteinsPositions, int seqSize, int form_size, vector<ClusterClass> &clusters);

	void GetFormsHash5(std::map<int, std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>>> &forms, vector<FormType> &formTypes, string &db, multimap<int, ProteinPosition> &proteinsPositions, int seqSize, int form_size);

	void GetFormsHash53(std::map<int, std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>>> &forms, vector<FormType> &formTypes, string &db, multimap<int, ProteinPosition> &proteinsPositions, int seqSize, int form_size, vector<ClusterClass> &clusters);

	void GetFormsHash6(std::map<int, std::map<int, std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>>>> &forms, vector<FormType> &formTypes, string &db, multimap<int, ProteinPosition> &proteinsPositions, int seqSize, int form_size);

	void CreateFormTypes(vector<FormType> &formTypes, int seqSize, int form_size, int mismatches);

	void InitializeHash();

	inline int AminoToIndexHash(char amino)
	{
		return AATabQuick[amino];
	}


};

#endif /* PREPROCESS_H_ */
