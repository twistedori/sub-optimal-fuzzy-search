

#include "FindSequence.h"

extern bool* in_cluster;

void FindSequenceClass::InitializeHashFind()
{
	int i, j;
	for (i = 0; i < 255; i++)
		AATabQuickFind[i] = 20;

	for (i = 0; i < 255; i++)
		for (j = 0; j < 20; j++)
			if (AAcode[j] == i)
				AATabQuickFind[i] = j;
}

void FindSequenceClass::changeAminoIndex(char changeFrom, char changeTo)
{
	AATabQuickFind[changeFrom] = AATabQuickFind[changeTo];
}

//Find sequence in a text
void FindSequenceClass::FindSequence(string &db, unordered_set<unsigned int> &positions, vector<unsigned int> &proteinsIndexes, std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>> &forms, vector<FormType> &formTypes, string &sequence, int form_size, int mismatches, int seqSize, unsigned int seed, vector<ClusterClass> &clusters)
{
	unordered_map<string, Form*> sequenceForms;

	//Geting forms from the sequence


	GetFormsSequence(sequenceForms, forms, formTypes, sequence, form_size);
	unordered_map<string, Form*>::iterator it = sequenceForms.begin();
	unordered_map<string, Form*>::iterator endIt = sequenceForms.end();
	Form* form;
	int index1, index2, index3, index4;
	int k, cluster_index = 0, count =0;
	unsigned int i, j, b, positionIndex, proteinNumber;
	unsigned char nonMatchCounter, matchCounter;
	unsigned char matches = seqSize - mismatches;
	int to_print = 0, pos_print = 0, timer = 0, deb =0;

	unsigned int positionsInSequenceSize, positionsInSequenceFormTypeSize, positionsInTextFormTypeSize, printsequence;
	while (clusters[cluster_index].getSeed() != seed) {
		cluster_index++;
	}

	while (it != endIt)
	{
		index1 = AminoToIndexHashFind(((*it).first)[0]);
		index2 = AminoToIndexHashFind(((*it).first)[1]);
		index3 = AminoToIndexHashFind(((*it).first)[2]);
		index4 = AminoToIndexHashFind(((*it).first)[3]);
		form = forms[index1][index2][index3][index4].at(0);
		if (form != NULL)
		{
	
			//Positions that the form appear in the db
			vector<vector<unsigned int>>* positionsInText = form->getPositions();
			//Positions that the form appear in the sequence
			vector<vector<unsigned int>>* positionsInSequence = (*it).second->getPositions();
			positionsInSequenceSize = positionsInSequence->size();
			//Go through all form types in the Form object
			for (unsigned int formTypeIndex = 0; formTypeIndex < positionsInSequenceSize; formTypeIndex++)
			{
				/*If the form types vector is not empty in both the db and the sequence
				it means that there is a possible match in the positions that appear
				in those vectors.*/
			
				if (!(*positionsInSequence)[formTypeIndex].empty())
					if(!(*positionsInText)[formTypeIndex].empty())
				{
					//Going through the positions of the form type
					positionsInSequenceFormTypeSize = (*positionsInSequence)[formTypeIndex].size();
					positionsInTextFormTypeSize = (*positionsInText)[formTypeIndex].size();
					//For each positions of the form type in the sequence
					for (i = 0; i < positionsInSequenceFormTypeSize; i++)
					{
						//For each positions of the form type in the db
						for (j = 0; j < positionsInTextFormTypeSize; j++)
						{
							/*Check the possibility that a positions in sequence is larger than the positions
							in the DB, skip this positions if it's true.      -----?? ?? ????*/
							if (((*positionsInText)[formTypeIndex])[j] >= ((*positionsInSequence)[formTypeIndex])[i])
							{
								/*Saving the position in the db that there could be a match starting from that position.*/
								positionIndex = ((*positionsInText)[formTypeIndex])[j] -((*positionsInSequence)[formTypeIndex])[i];
								nonMatchCounter = 0;
								matchCounter = 0;
								/*Go through the DB from the positionIndex and check if the substring there
								has no more than the requested mismatches.*/
								proteinNumber = proteinsIndexes[positionIndex];
								for (k = 0, b = positionIndex; k < seqSize && b < db.size() && proteinsIndexes[b] == proteinNumber; k++, b++)
								{
									if (db[b] != sequence[k])
										nonMatchCounter++;
									else
										matchCounter++;

									if (nonMatchCounter > mismatches)
										break;
								}
								if (matchCounter >= matches){
																	/*if we arrived to this location that means that we found a match within a threashold.
																	therefor we need to act accordingly. if the index//word is allready within a cluster:
																	y? continue
																	n? is it the first?
																
																	y? create a new cluster with this word as the seed.
																	n? add it to the cluster.*/
								
									if (!in_cluster[positionIndex]) {
										clusters[cluster_index].addPosition(positionIndex, mismatches/3);
										in_cluster[positionIndex] = true;
									

									
										

										}

									 }
								}
							}
						}
					}
				}
			}
		count++;

		it++;


		}
		
}


void FindSequenceClass::FindSequence5(string &db, unordered_set<unsigned int> &positions, vector<unsigned int> &proteinsIndexes, std::map<int, std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>>> &forms, vector<FormType> &formTypes, string &sequence, int form_size, int mismatches, int seqSize, unsigned int seed, vector<ClusterClass> &clusters)
{
	unordered_map<string, Form*> sequenceForms;
	
	//Geting forms from the sequence
	GetFormsSequence5(sequenceForms, forms, formTypes, sequence, form_size);
	unordered_map<string, Form*>::iterator it = sequenceForms.begin();
	unordered_map<string, Form*>::iterator endIt = sequenceForms.end();
	Form* form;
	int index1, index2, index3, index4, index5;
	int k, cluster_index = 0, count = 0;
	unsigned int i, j, b, positionIndex, proteinNumber;
	unsigned char nonMatchCounter, matchCounter;
	unsigned char matches = seqSize - mismatches;
	int to_print = 0, pos_print = 0, timer = 0, deb = 0;

	unsigned int positionsInSequenceSize, positionsInSequenceFormTypeSize, positionsInTextFormTypeSize, printsequence;
	while (clusters[cluster_index].getSeed() != seed) {
		cluster_index++;
	}

	while (it != endIt)
	{
		index1 = AminoToIndexHashFind(((*it).first)[0]);
		index2 = AminoToIndexHashFind(((*it).first)[1]);
		index3 = AminoToIndexHashFind(((*it).first)[2]);
		index4 = AminoToIndexHashFind(((*it).first)[3]);
		index5 = AminoToIndexHashFind(((*it).first)[4]);
		form = forms[index1][index2][index3][index4][index5].at(0);
		if (form != NULL)
		{
			
			
			//Positions that the form appear in the db
			vector<vector<unsigned int>>* positionsInText = form->getPositions();
			//Positions that the form appear in the sequence
			vector<vector<unsigned int>>* positionsInSequence = (*it).second->getPositions();
			positionsInSequenceSize = positionsInSequence->size();
			//Go through all form types in the Form object
			for (unsigned int formTypeIndex = 0; formTypeIndex < positionsInSequenceSize; formTypeIndex++)
			{
				/*If the form types vector is not empty in both the db and the sequence
				it means that there is a possible match in the positions that appear
				in those vectors.*/
			
				if (!(*positionsInSequence)[formTypeIndex].empty())
					if (!(*positionsInText)[formTypeIndex].empty())
					{
						//Going through the positions of the form type
						positionsInSequenceFormTypeSize = (*positionsInSequence)[formTypeIndex].size();
						positionsInTextFormTypeSize = (*positionsInText)[formTypeIndex].size();
						//For each positions of the form type in the sequence
						for (i = 0; i < positionsInSequenceFormTypeSize; i++)
						{
							//For each positions of the form type in the db
							for (j = 0; j < positionsInTextFormTypeSize; j++)
							{
								/*Check the possibility that a positions in sequence is larger than the positions
								in the DB, skip this positions if it's true.      -----?? ?? ????*/
								if (((*positionsInText)[formTypeIndex])[j] >= ((*positionsInSequence)[formTypeIndex])[i])
								{
									/*Saving the position in the db that there could be a match starting from that position.*/
									positionIndex = ((*positionsInText)[formTypeIndex])[j] - ((*positionsInSequence)[formTypeIndex])[i];
									nonMatchCounter = 0;
									matchCounter = 0;
									/*Go through the DB from the positionIndex and check if the substring there
									has no more than the requested mismatches.*/
									proteinNumber = proteinsIndexes[positionIndex];
									for (k = 0, b = positionIndex; k < seqSize && b < db.size() && proteinsIndexes[b] == proteinNumber; k++, b++)
									{
										if (db[b] != sequence[k])
											nonMatchCounter++;
										else
											matchCounter++;

										if (nonMatchCounter > mismatches)
											break;
									}
									if (matchCounter >= matches) {
										/*if we arrived to this location that means that we found a match within a threashold.
										therefor we need to act accordingly. if the index//word is allready within a cluster:
										y? continue
										n? is it the first?
										y? create a new cluster with this word as the seed.
										n? add it to the cluster.*/
										
										if (!in_cluster[positionIndex]) {
											clusters[cluster_index].addPosition(positionIndex, mismatches / 3);
											in_cluster[positionIndex] = true;
										





										}

									}
								}
							}
						}
					}
			}
		}
		count++;
		
		it++;

	}
}

void FindSequenceClass::FindSequence6(string &db, unordered_set<unsigned int> &positions, vector<unsigned int> &proteinsIndexes, std::map<int, std::map<int, std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>>>> &forms, vector<FormType> &formTypes, string &sequence, int form_size, int mismatches, int seqSize, unsigned int seed, vector<ClusterClass> &clusters)
{
	unordered_map<string, Form*> sequenceForms;
	//Geting forms from the sequence

	GetFormsSequence6(sequenceForms, forms, formTypes, sequence, form_size);
	unordered_map<string, Form*>::iterator it = sequenceForms.begin();
	unordered_map<string, Form*>::iterator endIt = sequenceForms.end();
	Form* form;
	int index1, index2, index3, index4, index5, index6;
	int k, cluster_index = 0, count = 0, times,c1,c2,c3;
	unsigned int i, j, b, positionIndex, proteinNumber;
	unsigned char nonMatchCounter, matchCounter;
	unsigned char matches = seqSize - mismatches;
	int to_print = 0, pos_print = 0, timer = 0, deb = 0;

	unsigned int positionsInSequenceSize, positionsInSequenceFormTypeSize, positionsInTextFormTypeSize, printsequence;
	while (clusters[cluster_index].getSeed() != seed) {
		cluster_index++;
	}

	while (it != endIt)
	{
		index1 = AminoToIndexHashFind(((*it).first)[0]);
		index2 = AminoToIndexHashFind(((*it).first)[1]);
		index3 = AminoToIndexHashFind(((*it).first)[2]);
		index4 = AminoToIndexHashFind(((*it).first)[3]);
		index5 = AminoToIndexHashFind(((*it).first)[4]);
		index6 = AminoToIndexHashFind(((*it).first)[5]);
		form = forms[index1][index2][index3][index4][index5][index6].at(0);
		if (form != NULL)
		{
			
			//Positions that the form appear in the db
			vector<vector<unsigned int>>* positionsInText = form->getPositions();
			//Positions that the form appear in the sequence
			vector<vector<unsigned int>>* positionsInSequence = (*it).second->getPositions();
			positionsInSequenceSize = positionsInSequence->size();
			//Go through all form types in the Form object
			for (unsigned int formTypeIndex = 0,c1=0; formTypeIndex < positionsInSequenceSize; formTypeIndex++,c1++)
			{
				/*If the form types vector is not empty in both the db and the sequence
				it means that there is a possible match in the positions that appear
				in those vectors.*/
				
				if (!(*positionsInSequence)[formTypeIndex].empty())
					if (!(*positionsInText)[formTypeIndex].empty())
					{
						//Going through the positions of the form type
						positionsInSequenceFormTypeSize = (*positionsInSequence)[formTypeIndex].size();
						positionsInTextFormTypeSize = (*positionsInText)[formTypeIndex].size();
						//For each positions of the form type in the sequence
						for (c2=0,i = 0; i < positionsInSequenceFormTypeSize; i++,c2++)
						{
							//For each positions of the form type in the db
							for (c3=0,j = 0; j < positionsInTextFormTypeSize; j++,c3++)
							{
								/*Check the possibility that a positions in sequence is larger than the positions
								in the DB, skip this positions if it's true.      -----?? ?? ????*/
								if (((*positionsInText)[formTypeIndex])[j] >= ((*positionsInSequence)[formTypeIndex])[i])
								{
									/*Saving the position in the db that there could be a match starting from that position.*/
									positionIndex = ((*positionsInText)[formTypeIndex])[j] - ((*positionsInSequence)[formTypeIndex])[i];
									nonMatchCounter = 0;
									matchCounter = 0;
									/*Go through the DB from the positionIndex and check if the substring there
									has no more than the requested mismatches.*/
					
									for (times=0,k = 0, b = positionIndex; k < seqSize && b < db.size(); k++, b++,times++)
									{
								
										if (db[b] != sequence[k])
											nonMatchCounter++;
										else
											matchCounter++;

										if (nonMatchCounter > mismatches)
											break;
									}
									if (matchCounter >= matches) {
											positions.insert(positionIndex);/*if we arrived to this location that means that we found a match within a threashold.
										therefor we need to act accordingly. if the index//word is allready within a cluster:
										y? continue
										n? is it the first?

										y? create a new cluster with this word as the seed.
										n? add it to the cluster.*/
									
										if (!in_cluster[positionIndex]) {
											clusters[cluster_index].addPosition(positionIndex, mismatches);
											in_cluster[positionIndex] = true;
																		
										}

									}
								}
							}
						}
					}
			}
		}
		count++;
		
		it++;
	

	}

}

void FindSequenceClass::FindPossibilities(string & db, unordered_set<unsigned int>& positions, vector<unsigned int>& proteinsIndexes, std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>>& forms, vector<FormType>& formTypes, string & sequence, int form_size, int mismatches, int seqSize, unsigned int seed, vector<ClusterClass>& clusters)
{
	unordered_map<string, Form*> sequenceForms;
	//Geting forms from the sequence

	GetFormsSequence(sequenceForms, forms, formTypes, sequence, form_size);
	unordered_map<string, Form*>::iterator it = sequenceForms.begin();
	unordered_map<string, Form*>::iterator endIt = sequenceForms.end();
	Form* form;
	int index1, index2, index3, index4;
	int printer[15];
	static int stat = 0;
	int k, cluster_index = 0, count = 0;
	unsigned int i, j, b, positionIndex, proteinNumber;
	unsigned char nonMatchCounter, matchCounter;
	unsigned char matches = seqSize - mismatches;
	int to_print = 0, pos_print = 0, timer = 0, deb = 0;
		unsigned int positionsInSequenceSize, positionsInSequenceFormTypeSize, positionsInTextFormTypeSize, printsequence;
	while (clusters[cluster_index].getSeed() != seed) {
		cluster_index++;
	}

	while (count < sequenceForms.size())
	{
		index1 = AminoToIndexHashFind(((*it).first)[0]);
		index2 = AminoToIndexHashFind(((*it).first)[1]);
		index3 = AminoToIndexHashFind(((*it).first)[2]);
		index4 = AminoToIndexHashFind(((*it).first)[3]);
		if( forms[index1][index2][index3][index4].empty())
			form = NULL;
		else(form = forms[index1][index2][index3][index4].at(0));
		if (form != NULL)
		{
						
			//Positions that the form appear in the db
			vector<vector<unsigned int>>* positionsInText = form->getPositions();
			//Positions that the form appear in the sequence
			vector<vector<unsigned int>>* positionsInSequence = (*it).second->getPositions();
			positionsInSequenceSize = positionsInSequence->size();
			//Go through all form types in the Form object
			for (unsigned int formTypeIndex = 0; formTypeIndex < positionsInSequenceSize; formTypeIndex++)
			{
				/*If the form types vector is not empty in both the db and the sequence
				it means that there is a possible match in the positions that appear
				in those vectors.*/
					if (!(*positionsInSequence)[formTypeIndex].empty())
					if (!(*positionsInText)[formTypeIndex].empty())
					{
						//Going through the positions of the form type
						positionsInSequenceFormTypeSize = (*positionsInSequence)[formTypeIndex].size();
						positionsInTextFormTypeSize = (*positionsInText)[formTypeIndex].size();
//For each positions of the form type in the sequence
for (i = 0; i < positionsInSequenceFormTypeSize; i++)
{
	//For each positions of the form type in the db
	for (j = 0; j < positionsInTextFormTypeSize; j++)
	{
		
		/*Check the possibility that a positions in sequence is larger than the positions
		in the DB, skip this positions if it's true.      -----?? ?? ????*/
		if (((*positionsInText)[formTypeIndex])[j] >= ((*positionsInSequence)[formTypeIndex])[i])
		{
			
			/*Saving the position in the db that there could be a match starting from that position.*/
			positionIndex = ((*positionsInText)[formTypeIndex])[j] - ((*positionsInSequence)[formTypeIndex])[i];
			if (positionIndex >= clusters.size()) continue;
			nonMatchCounter = 0;
			matchCounter = 0;
			/*Go through the DB from the positionIndex and check if the substring there
			has no more than the requested mismatches.*/
			proteinNumber = proteinsIndexes[positionIndex];
			auto search = seedHash.find(positionIndex);

			if (posibleConnections[cluster_index][search->second] || posibleConnections[search->second][cluster_index] || search->second == cluster_index) continue;
			for (k = 0, b = positionIndex; k < seqSize && b < db.size() ; k++, b++)
			{
				posibleConnections[cluster_index][search->second] = true;
				posibleConnections[search->second][cluster_index] = true;

				if (db[b] != sequence[k])
					nonMatchCounter++;
				else
					matchCounter++;

				if (nonMatchCounter > mismatches)
					break;
			}
			if (nonMatchCounter < mismatches) {
				clusters[cluster_index].addPosibilitie(positionIndex);
				clusters[(search->second)].addPosibilitie(seed);
			
			}

		}
	}
}
						}
					}
			}
			count++;
			
			it++;


}

	}



	void FindSequenceClass::FindConnections(string & db, unordered_set<unsigned int>& positions, vector<unsigned int>& proteinsIndexes, std::map<int, std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>>>& forms, vector<FormType>& formTypes, string & sequence, int form_size, int mismatches, int seqSize, unsigned int seed, vector<ClusterClass>& clusters)
	{

		

		vector<unsigned int> posits;
		vector< unsigned int > posib_posits;
		vector<unsigned int> possibilities;
		int p_index, p_index2, matchCounter;
		int indexes[15], i;
		unsigned int tmp1, tmp2;
		for (i = 0; i < 15; i++)indexes[i] = 0;
		auto search = seedHash.find(seed);
		tmp1 = search->second;
		posits = clusters[tmp1].getPositions();
		possibilities = clusters[tmp1].getPosibilities();
		for (p_index = 0; p_index < possibilities.size(); p_index++) {

			auto search2 = seedHash.find(possibilities[p_index]);
			tmp2 = search2->second;
			if (Connections[tmp1][tmp2] || Connections[tmp2][tmp1]) continue;
			posib_posits = clusters[tmp2].getPositions();
			for (indexes[0]; indexes[0] < posits.size(); indexes[0]++)
				for (indexes[1]; indexes[1] < posib_posits.size(); indexes[1]++)
				{
					matchCounter = 0;
					for (indexes[2]; indexes[2] < seqSize; indexes[2]++)
					{
						if (db[indexes[0] + indexes[2]] == db[indexes[1] + indexes[2]]) matchCounter++;
					}
					if (matchCounter >= seqSize - mismatches) {
						Connections[tmp1][tmp2];
						Connections[tmp2][tmp1];
						clusters[tmp1].addConnection(tmp2);
						clusters[tmp2].addConnection(tmp1);

				}


			}
		}
	}



