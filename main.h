#pragma once
#ifndef MAIN_H_
#define MAIN_H_

#include <iostream>
#include <string>
#include <unordered_map>
#include <fstream>
#include <sstream>

#include <chrono>
#include "Form.h"
#include "Preprocess.h"
#include "FindSequence.h"
#include "cluster.h"
#include <windows.h>
#include <direct.h>
#include <Shlwapi.h>

using namespace std;



struct SequenceResults
{
	string sequence;
	unordered_set<unsigned int> positions;
};


int SEQ_SIZE = 20;
const int SRQ_SIZE2 = 10;

void startFindAlgorithm(string parameters);

void ReadDB(string path, string fileName);

void InitializePreprocess(int mismatches, int form_size);

void InitializePreprocess2(int mismatches, int form_size);

void InitializePreprocess3(int mismatches, int form_size);

void GetSequences(string sequencesPath, string sequencesFileName, int mismatches, int form_size, vector<ClusterClass> &clusters);

void PossibleMatchAlgorithm(string parameters);

void PossibleClusters40(string sequencesPath, string sequencesFileName, int mismatches, int form_size, vector<ClusterClass> &clusters);

void connectionsAlgorithm(string parameters);

void clustersConnection(string sequencesPath, string sequencesFileName, int mismatches, int form_size, vector<ClusterClass> &clusters);

int main(int argc, char *argv[]);


#endif /* MAIN_H_ */
