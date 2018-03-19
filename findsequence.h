#pragma once


#ifndef FINDSEQUENCE_H_
#define FINDSEQUENCE_H_


#include <unordered_map>
#include <unordered_set>
#include <set>
#include <iostream>
#include "Form.h"
#include "cluster.h"
extern vector<ClusterClass> clusters;

extern bool** posibleConnections;
extern bool** Connections;
extern multimap<unsigned int, string> seedSequences;
extern multimap<unsigned int, string> seedSequences2;
extern multimap<unsigned int, unsigned int> seedHash;


using namespace std;

class FindSequenceClass
{
private:
	int AATabQuickFind[255];
public:
	void InitializeHashFind();

	inline int AminoToIndexHashFind(char amino)
	{
		return AATabQuickFind[amino];
	}

	void changeAminoIndex(char changeFrom, char changeTo);

	inline void GetFormsSequence(unordered_map<string, Form*> &sequenceForms, std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>> &forms, vector<FormType> &formTypes, string &text, int form_size)
	{
		int cnt;
		string formSeq;
		int distance;
		Form* form = nullptr;
		unordered_map<string, Form*>::iterator it;
		int i, j, k, b, limit;
		int index1, index2, index3, index4;
		//Loop through all form types
		int formTypesSize = formTypes.size();
		unsigned int textSize = text.size();
		int distancesSize;
		formSeq.resize(form_size);
	
		for (j = 0; j < formTypesSize; j++)
		{
			//Create the forms with their position in the text
			limit = textSize - formTypes[j].formType.size() + 1;
			for (i = 0; i < limit; i++)
			{
				//Get the form
				distance = 0;
				b = 0;
				distancesSize = formTypes[j].distances.size();
				for (k = 0; k < distancesSize; k++)
				{
					distance += formTypes[j].distances[k];
					formSeq[b] = text[i + distance];
					b++;
				}
				//Check if the form already exists in the form list, if not then create a new form, else add its appearence position in the text
				it = sequenceForms.find(formSeq);
			
	
				if (it == sequenceForms.end())
				{
					Form* newForm = new Form(formTypes.size(), j, i);
					pair<string, Form*> tempPair = make_pair(formSeq, newForm);
					sequenceForms.insert(tempPair);
				}
				else
					(*it).second->addPosition(j, i);
			}
		}
	}

	inline void GetFormsSequence5(unordered_map<string, Form*> &sequenceForms, std::map<int, std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>>> &forms, vector<FormType> &formTypes, string &text, int form_size)
	{
		string formSeq;
		int distance;
		Form* form;
		unordered_map<string, Form*>::iterator it;
		int i, j, k, b, limit;
		int index1, index2, index3, index4, index5;
		//Loop through all form types
		int formTypesSize = formTypes.size();
		unsigned int textSize = text.size();
		int distancesSize;
		formSeq.resize(form_size);
		for (j = 0; j < formTypesSize; j++)
		{
			//Create the forms with their position in the text
			limit = textSize - formTypes[j].formType.size() + 1;
			for (i = 0; i < limit; i++)
			{
				//Get the form
				distance = 0;
				b = 0;
				distancesSize = formTypes[j].distances.size();
				for (k = 0; k < distancesSize; k++)
				{
					distance += formTypes[j].distances[k];
					formSeq[b] = text[i + distance];
					b++;
				}
				//Check if the form already exists in the form list, if not then create a new form, else add its appearence position in the text
				it = sequenceForms.find(formSeq);
				if (it == sequenceForms.end())
				{
					Form *newForm = new Form(formTypes.size(), j, i);
					pair<string, Form*> tempPair = make_pair(formSeq, newForm);
					sequenceForms.insert(tempPair);
				}
				else
					(*it).second->addPosition(j, i);
			}
		}
	}

	inline void GetFormsSequence6(unordered_map<string, Form*> &sequenceForms, std::map<int, std::map<int, std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>>>> &forms, vector<FormType> &formTypes, string &text, int form_size)
	{
		int cnt;
		string formSeq;
		int distance;
		Form* form = nullptr;
		unordered_map<string, Form*>::iterator it;
		int i, j, k, b, limit;
		int index1, index2, index3, index4, index5, index6;
		//Loop through all form types
		int formTypesSize = formTypes.size();
		unsigned int textSize = text.size();
		int distancesSize;
		formSeq.resize(form_size);
	
		for (j = 0; j < formTypesSize; j++)
		{
			//Create the forms with their position in the text
			limit = textSize - formTypes[j].formType.size() + 1;
			for (i = 0; i < limit; i++)
			{
				//Get the form
				distance = 0;
				b = 0;
				distancesSize = formTypes[j].distances.size();
				for (k = 0; k < distancesSize; k++)
				{
					distance += formTypes[j].distances[k];
					formSeq[b] = text[i + distance];
					b++;
				}
				//Check if the form already exists in the form list, if not then create a new form, else add its appearence position in the text
				it = sequenceForms.find(formSeq);
			
				if (it == sequenceForms.end())
				{
					Form* newForm = new Form(formTypes.size(), j, i);
					pair<string, Form*> tempPair = make_pair(formSeq, newForm);
					sequenceForms.insert(tempPair);
				}
				else
					(*it).second->addPosition(j, i);
			}
		}
	}

	void FindSequence(string &db, unordered_set<unsigned int> &positions, vector<unsigned int> &proteinsIndexes, std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>> &forms, vector<FormType> &formTypes, string &sequence, int form_size, int mismatches, int seqSize, unsigned int seed, vector<ClusterClass> &clusters);

	void FindSequence5(string &db, unordered_set<unsigned int> &positions, vector<unsigned int> &proteinsIndexes, std::map<int, std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>>> &forms, vector<FormType> &formTypes, string &sequence, int form_size, int mismatches, int seqSize, unsigned int seed, vector<ClusterClass> &clusters);

	void FindSequence6(string &db, unordered_set<unsigned int> &positions, vector<unsigned int> &proteinsIndexes, std::map<int, std::map<int, std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>>>> &forms, vector<FormType> &formTypes, string &sequence, int form_size, int mismatches, int seqSize, unsigned int seed, vector<ClusterClass> &clusters);

	void FindPossibilities(string &db, unordered_set<unsigned int> &positions, vector<unsigned int> &proteinsIndexes, std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>> &forms, vector<FormType> &formTypes, string &sequence, int form_size, int mismatches, int seqSize, unsigned int seed, vector<ClusterClass> &clusters);

	void FindConnections(string &db, unordered_set<unsigned int> &positions, vector<unsigned int> &proteinsIndexes, std::map<int, std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>>> &forms, vector<FormType> &formTypes, string &sequence, int form_size, int mismatches, int seqSize, unsigned int seed, vector<ClusterClass> &clusters);
};

#endif /* FINDSEQUENCE_H_ */
