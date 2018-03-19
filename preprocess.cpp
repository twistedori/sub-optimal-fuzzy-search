#include "Preprocess.h"
#include "cluster.h"

bool PreprocessClass::CheckIfPreProcessed(int form_size, int mismatches, int seqSize)
{
	stringstream filename;
	filename << seqSize << "_" << form_size << "_" << mismatches << ".txt";
	ifstream myfile(filename.str());
	return myfile.good();
}

void PreprocessClass::InitializeHash()
{
	int i, j;
	for (i = 0; i < 255; i++)
		AATabQuick[i] = 20;

	for (i = 0; i < 255; i++)
		for (j = 0; j < 20; j++)
			if (AAcode[j] == i)
				AATabQuick[i] = j;
}

void PreprocessClass::CreateFormTypes(vector<FormType> &formTypes, int seqSize, int form_size, int mismatches)
{
	stringstream filename, ss;
	int i, j;
	int lastOne;
	filename << seqSize << "_" << form_size << "_" << mismatches << ".txt";
	ifstream myfile;
	myfile.open(filename.str());
	string line, form, formType;
	vector<string> formTypesStrings;
	while (getline(myfile, line) && line != "")
		formTypesStrings.push_back(line);
	myfile.close();

	for (i = 0; i < formTypesStrings.size(); i++)
	{
		FormType formType;
		formType.formType = formTypesStrings[i];
		lastOne = 0;
		for (j = 0; j < formTypesStrings[i].size(); j++)
		{
			if ((formTypesStrings[i])[j] == '1')
			{
				formType.distances.push_back(j - lastOne);
				lastOne = j;
			}
		}
		formTypes.push_back(formType);
	}
}



void PreprocessClass::PreProcess(std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>> &forms, vector<FormType> &formTypes, int form_size, int mismatches, string &db, multimap<int, ProteinPosition> &proteinsPositions, int seqSize)
{
	CreateFormTypes(formTypes, seqSize, form_size, mismatches);
	GetFormsHash(forms, formTypes, db, proteinsPositions, seqSize, form_size);
}

void PreprocessClass::PreProcess42(std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>> &forms, std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>> &forms2, vector<FormType> &formTypes, int form_size, int mismatches, string &db, multimap<int, ProteinPosition> &proteinsPositions, int seqSize, vector<ClusterClass> &clusters)
{
	CreateFormTypes(formTypes, seqSize, form_size, mismatches);

	GetFormsHash42(forms, forms2, formTypes, db, proteinsPositions, seqSize, form_size, clusters);


}

void PreprocessClass::PreProcess53(std::map<int, std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>>> &forms, vector<FormType> &formTypes, int form_size, int mismatches, string &db, multimap<int, ProteinPosition> &proteinsPositions, int seqSize, vector<ClusterClass> &clusters) {
	CreateFormTypes(formTypes, seqSize, form_size, mismatches);
	GetFormsHash53(forms, formTypes, db, proteinsPositions, seqSize, form_size, clusters);
}

void PreprocessClass::PreProcess5(std::map<int, std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>>> &forms, vector<FormType> &formTypes, int form_size, int mismatches, string &db, multimap<int, ProteinPosition> &proteinsPositions, int seqSize)
{
	CreateFormTypes(formTypes, seqSize, form_size, mismatches);
	GetFormsHash5(forms, formTypes, db, proteinsPositions, seqSize, form_size);
	
}

void PreprocessClass::PreProcess6(std::map<int, std::map<int, std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>>>> &forms, vector<FormType> &formTypes, int form_size, int mismatches, string &db, multimap<int, ProteinPosition> &proteinsPositions, int seqSize)
{
	CreateFormTypes(formTypes, seqSize, form_size, mismatches);
	GetFormsHash6(forms, formTypes, db, proteinsPositions, seqSize, form_size);

}



void PreprocessClass::GetFormsHash(std::map<int, std::map<int, std::map<int, std::map<int, vector<Form*>>>>> &forms, vector<FormType> &formTypes, string &db, multimap<int, ProteinPosition> &proteinsPositions, int seqSize, int form_size)
{
	int distance, index1, index2, index3, index4;
	multimap<int, ProteinPosition>::iterator it2;
	std::vector<Form *>::iterator it;
	bool skip = false;
	string formSeq;
	formSeq.resize(form_size);
	unsigned int i, j, k, b, limit;
	for (j = 0; j < formTypes.size(); j++)
	{
		//Create the forms with their position in the text according to the proteins positions
		it2 = proteinsPositions.begin();
		while (it2 != proteinsPositions.end())
		{
			if (((*it2).second.endIndex - (*it2).second.startIndex) < formTypes[j].formType.size())
			{
				//cout << protein << endl;
				it2++;
				continue;
			}
			limit = ((*it2).second.endIndex - formTypes[j].formType.size());
			for (i = (*it2).second.startIndex; i <= limit; i++)
			{
				//Get the form
				distance = 0;
				b = 0;
				for (k = 0; k < formTypes[j].distances.size(); k++)
				{
					distance += formTypes[j].distances[k];
					if (AminoToIndexHash(db[i + distance]) == 20)
					{
						skip = true;
						break;
					}
					else
					{
						formSeq[b] = db[i + distance];
						b++;
					}
				}
				if (!skip)
				{
					//Check if the form already exists in the form list, if not then create a new form, else add its appearance position in the text
					index1 = AminoToIndexHash(formSeq[0]);
					index2 = AminoToIndexHash(formSeq[1]);
					index3 = AminoToIndexHash(formSeq[2]);
					index4 = AminoToIndexHash(formSeq[3]);
					if (empty(forms[index1][index2][index3][index4])) {
						Form* temp = new Form(formTypes.size(), j, i);
						it = forms[index1][index2][index3][index4].begin();
						it = forms[index1][index2][index3][index4].insert(it, temp);
						it = forms[index1][index2][index3][index4].end();
					}
					else
						forms[index1][index2][index3][index4].at(0)->addPosition(j, i);

				}
				skip = false;
			}
			it2++;
		}
	}
}


void PreprocessClass::GetFormsHash42(std::map<int, std::map<int, std::map<int, std::map<int, vector<Form*>>>>> &forms, std::map<int, std::map<int, std::map<int, std::map<int, vector<Form*>>>>> &forms2, vector<FormType> &formTypes, string &db, multimap<int, ProteinPosition> &proteinsPositions, int seqSize, int form_size, vector<ClusterClass> &clusters)
{
	int distance, index1, index2, index3, index4;
	multimap<int, ProteinPosition>::iterator it2;
	std::vector<Form *>::iterator it;
	bool skip = false;
	string formSeq;
	formSeq.resize(form_size);
	unsigned int i, j, k, b, limit, clusterIndex = 0;
	for (j = 0; j < formTypes.size(); j++)
	{
		//Create the forms with their position in the text according to the proteins positions
		it2 = proteinsPositions.begin();
		while (clusterIndex < clusters.size())
		{
			printf("clusterNUM - %d,   out of  %d\n", clusterIndex, clusters.size());
			limit = (clusters[clusterIndex].getSeed() + seqSize - formTypes[j].formType.size());
			for (i = clusters[clusterIndex].getSeed(); i <= limit; i++)
			{
				//Get the form
				distance = 0;
				b = 0;
				for (k = 0; k < formTypes[j].distances.size(); k++)
				{
					distance += formTypes[j].distances[k];
					if (AminoToIndexHash(db[i + distance]) == 20)
					{
						skip = true;
						break;
					}
					else
					{
						formSeq[b] = db[i + distance];
						b++;
					}
				}
				if (!skip)
				{
					//Check if the form already exists in the form list, if not then create a new form, else add its appearance position in the text
					index1 = AminoToIndexHash(formSeq[0]);
					index2 = AminoToIndexHash(formSeq[1]);
					index3 = AminoToIndexHash(formSeq[2]);
					index4 = AminoToIndexHash(formSeq[3]);
					if (empty(forms[index1][index2][index3][index4])) {
						Form* temp = new Form(formTypes.size(), j, i);
						it = forms[index1][index2][index3][index4].begin();
						it = forms[index1][index2][index3][index4].insert(it, temp);
						it = forms[index1][index2][index3][index4].end();
					}
					else
					forms[index1][index2][index3][index4].at(0)->addPosition(j, i);
									
				}
				skip = false;
			}
		
			clusterIndex++;
		}
	}
}

void PreprocessClass::GetFormsHash5(std::map<int, std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>>> &forms, vector<FormType> &formTypes, string &db, multimap<int, ProteinPosition> &proteinsPositions, int seqSize, int form_size)
{
	int distance, index1, index2, index3, index4, index5;
	multimap<int, ProteinPosition>::iterator it2;
	std::vector<Form *>::iterator it;
	bool skip = false;
	string formSeq;
	formSeq.resize(form_size);
	unsigned int i, j, k, b, limit;
	for (j = 0; j < formTypes.size(); j++)
	{
		//Create the forms with their position in the text according to the proteins positions
		it2 = proteinsPositions.begin();
		while (it2 != proteinsPositions.end())
		{
			if (((*it2).second.endIndex - (*it2).second.startIndex) < formTypes[j].formType.size())
			{
				it2++;
				continue;
			}
			limit = ((*it2).second.endIndex - formTypes[j].formType.size());
			for (i = (*it2).second.startIndex; i <= limit; i++)
			{
				//Get the form
				distance = 0;
				b = 0;
				for (k = 0; k < formTypes[j].distances.size(); k++)
				{
					distance += formTypes[j].distances[k];
					if (AminoToIndexHash(db[i + distance]) == 20)
					{
						skip = true;
						break;
					}
					else
					{
						formSeq[b] = db[i + distance];
						b++;
					}
				}
				if (!skip)
				{
					//Check if the form already exists in the form list, if not then create a new form, else add its appearance position in the text
					
					index1 = AminoToIndexHash(formSeq[0]);
					index2 = AminoToIndexHash(formSeq[1]);
					index3 = AminoToIndexHash(formSeq[2]);
					index4 = AminoToIndexHash(formSeq[3]);
					index5 = AminoToIndexHash(formSeq[4]);
					if (empty(forms[index1][index2][index3][index4][index5])) {
						it = forms[index1][index2][index3][index4][index5].begin();
						it = forms[index1][index2][index3][index4][index5].insert(it, new Form(formTypes.size(), j, i));
						it = forms[index1][index2][index3][index4][index5].end();
					}
					else
						forms[index1][index2][index3][index4][index5].at(0)->addPosition(j, i);
				}
				skip = false;
			}
			it2++;
		}
	}
}


void PreprocessClass::GetFormsHash6(std::map<int, std::map<int, std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>>>> &forms, vector<FormType> &formTypes, string &db, multimap<int, ProteinPosition> &proteinsPositions, int seqSize, int form_size)
{
	int distance, index1, index2, index3, index4, index5, index6;
	multimap<int, ProteinPosition>::iterator it2;
	std::vector<Form *>::iterator it;
	bool skip = false;
	string formSeq;
	formSeq.resize(form_size);
	unsigned int i, j, k, b, limit, cnt = 0;
	for (j = 0; j < formTypes.size(); j++)
	{
		
		for (i = 0; i <= db.size() - 6; i++) {
			distance = 0;
			b = 0;
			for (k = 0; k < formTypes[j].distances.size(); k++)
			{
				distance += formTypes[j].distances[k];
				if (AminoToIndexHash(db[i + distance]) == 20)
				{
					skip = true;
					break;
				}
				else
				{
					formSeq[b] = db[i + distance];
					b++;
				}
			}
			if (!skip)
			{
				cnt++;
				if (!(cnt % 100000))
					printf("%d \n", cnt);
				//Check if the form already exists in the form list, if not then create a new form, else add its appearance position in the text
				index1 = AminoToIndexHash(formSeq[0]);
				index2 = AminoToIndexHash(formSeq[1]);
				index3 = AminoToIndexHash(formSeq[2]);
				index4 = AminoToIndexHash(formSeq[3]);
				index5 = AminoToIndexHash(formSeq[4]);
				index6 = AminoToIndexHash(formSeq[5]);
				if (empty(forms[index1][index2][index3][index4][index5][index6])) {
					it = forms[index1][index2][index3][index4][index5][index6].begin();
					it = forms[index1][index2][index3][index4][index5][index6].insert(it, new Form(formTypes.size(), j, i));
					it = forms[index1][index2][index3][index4][index5][index6].end();
				}
				else
					forms[index1][index2][index3][index4][index5][index6].at(0)->addPosition(j, i);
			}
			skip = false;
		}
	


	}
}






void PreprocessClass::GetFormsHash53(std::map<int, std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>>> &forms, vector<FormType> &formTypes, string &db, multimap<int, ProteinPosition> &proteinsPositions, int seqSize, int form_size, vector<ClusterClass> &clusters)
{
	int distance, index1, index2, index3, index4, index5;
	multimap<int, ProteinPosition>::iterator it2;
	std::vector<Form *>::iterator it;
	std::vector<unsigned int>::iterator it3;
	vector<unsigned int> dist;
	bool skip = false;
	string formSeq;
	formSeq.resize(form_size);
	unsigned int i, j, k, b, limit, limit2, clusterIndex = 0, distcnt;
	for (j = 0; j < formTypes.size(); j++)
	{
		//Create the forms with their position in the text according to the proteins positions
		while (clusterIndex < clusters.size())
		{
			distcnt = 0;
			dist = clusters[clusterIndex].getPositions();

			while (distcnt <= dist.size())
			{
				
				limit = (dist[distcnt] + seqSize - formTypes[j].formType.size());
				for (i = dist[distcnt]; i <= limit; i++)
				{
					//Get the form
					distance = 0;
					b = 0;
					for (k = 0; k < formTypes[j].distances.size(); k++)
					{
						distance += formTypes[j].distances[k];
						if (AminoToIndexHash(db[i + distance]) == 20)
						{
							skip = true;
							break;
						}
						else
						{
							formSeq[b] = db[i + distance];
							b++;
						}
					}
					if (!skip)
					{
						//Check if the form already exists in the form list, if not then create a new form, else add its appearance position in the text
						index1 = AminoToIndexHash(formSeq[0]);
						index2 = AminoToIndexHash(formSeq[1]);
						index3 = AminoToIndexHash(formSeq[2]);
						index4 = AminoToIndexHash(formSeq[3]);
						index5 = AminoToIndexHash(formSeq[4]);
						if (empty(forms[index1][index2][index3][index4][index5]))
						{
							it = forms[index1][index2][index3][index4][index5].begin();
							it = forms[index1][index2][index3][index4][index5].insert(it, new Form(formTypes.size(), j, i));
							it = forms[index1][index2][index3][index4][index5].end();
						}
						else
							forms[index1][index2][index3][index4][index5].at(0)->addPosition(j, i);
					}
					skip = false;


				}
				distcnt++;
			}


			clusterIndex++;
		}
	}
}
	