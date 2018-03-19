/*Created  by Ori Erel ------- main function starts at line 504*/
#include "Main.h"

string db;
multimap<int, ProteinPosition> proteinsPositions;
vector<unsigned int> proteinsIndexes;
vector<string> proteinsNames;
int mismatches, form_size;
vector<string> sequences;
multimap<unsigned int, string> seedSequences;
multimap<unsigned int, string> seedSequences2;

bool* in_cluster; //this will be the array that simbolizes weather a certain index is allready within a cluster
bool* in_cluster2;
PreprocessClass preprocessClass;
FindSequenceClass findSeqClass;

multimap<unsigned int, unsigned int> seedHash;

bool** posibleConnections;
bool** Connections;
vector<ClusterClass> clusters;



std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>> forms4;

std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>> forms42;

std::map<int, std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>>> forms5;

std::map<int, std::map<int, std::map<int, std::map<int, std::map<int, std::map<int, std::vector<Form*>>>>>>> forms6;




vector<FormType> formTypes;


void ReadDB(string path, string fileName)
{
	unsigned int i;
	ifstream myfile;
	stringstream buffer, directory;
	string line;
	directory << "ProteinDB4.txt";
	//myfile.open(path + fileName);
	myfile.open(directory.str());
	if (!myfile.good())
	{
		cout << "DB file does not exist" << endl;
		cin.get();
		exit(1);
	}

	cout << "DB read started" << endl;
	int proteinIndex = 0;
	int proteinPosition = 0;
	int proteinSize = 0;
	vector <string> name;
	bool reachedSize = false;
	string temp;
	
	//Read the database and add each occurence and its name to the vectors
	while (myfile >> line)
	{

		if (line == "sz")
		reachedSize = true;
		if (reachedSize)
		{
			myfile >> proteinSize;
			myfile >> temp;
			i = 0;
			while (temp.size() < proteinSize || temp.substr(temp.size() - 3, temp.size()) == "fas")
			{
				myfile >> temp;
				if (temp.size() < proteinSize || i == 0)
				{
					name.push_back(temp);
					if (i == 0 && temp.size() >= proteinSize)
						myfile >> temp;
				}
				i++;
			}
			proteinsNames.push_back(name[0]);
			for (i = 0; i < name.size(); i++)
			{
				proteinsNames[proteinsNames.size() - 1].append(" ");
				proteinsNames[proteinsNames.size() - 1].append(name[i]);
			}
			name.clear();
			proteinPosition = db.size();
			db.append(temp);
			reachedSize = false;



			ProteinPosition newPositions;
			newPositions.startIndex = proteinPosition;
			newPositions.endIndex = db.size() - 1;
			pair<int, ProteinPosition> tempPair2 = make_pair(proteinIndex, newPositions);
			proteinsPositions.insert(tempPair2);

			for (i = 0; i < proteinSize; i++)
				proteinsIndexes.push_back(proteinIndex);

			proteinIndex++;

		}
	}
	myfile.close();
	cout << "DB read finished" << endl;
}

void InitializePreprocess(int mismatches, int form_size)
{

	if (form_size != 4 && form_size != 5 && form_size != 6)
	{
		cout << "Form size is invalid" << endl;
		cin.get();
		exit(1);
	}

	/*Initialize ASCII hash, both functions are doing exactly the same thing
	but they are divided between 2 classes and are used by inline functions.*/
	preprocessClass.InitializeHash();
	findSeqClass.InitializeHashFind();

	//Check if there is already a file that contain the preprocessed form types results
	cout << "Preprocess started" << endl;
	//Check if a form types set file exist
	bool preprocessed = preprocessClass.CheckIfPreProcessed(form_size, mismatches, SEQ_SIZE);
	if (!preprocessed)
	{
		cout << "Form types file is missing" << endl;
		cin.get();
		exit(1);
	}
	else
	{
		if (form_size == 4)
			preprocessClass.PreProcess(forms4, formTypes, form_size, mismatches, db, proteinsPositions, SEQ_SIZE);
		if (form_size == 5)
			preprocessClass.PreProcess5(forms5, formTypes, form_size, mismatches, db, proteinsPositions, SEQ_SIZE);
		if (form_size == 6);
			preprocessClass.PreProcess6(forms6, formTypes, form_size, mismatches, db, proteinsPositions, SEQ_SIZE);
	}

	cout << "Preprocess finished" << endl;
}


void InitializePreprocess2(int mismatches, int form_size)
{
	unsigned int clusterSize, i,j;
	pair<unsigned int, unsigned int> hashPair;
	vector<unsigned int> temp;
	
	if (form_size != 4 && form_size != 5 && form_size != 6)
	{
		cout << "Form size is invalid" << endl;
		cin.get();
		exit(1);
	}

	/*Initialize ASCII hash, both functions are doing exactly the same thing
	but they are divided between 2 classes and are used by inline functions.*/
	preprocessClass.InitializeHash();
	findSeqClass.InitializeHashFind();

	
	//Check if there is already a file that contain the preprocessed form types results
	cout << "Preprocess2 started" << endl;
	//Check if a form types set file exist
	bool preprocessed = preprocessClass.CheckIfPreProcessed(form_size, mismatches, SEQ_SIZE);
	if (!preprocessed)
	{
		cout << "Form types file is missing" << endl;
		cin.get();
		exit(1);
	}
	else
	{
		if (form_size == 4)
			preprocessClass.PreProcess42(forms4, forms42, formTypes, form_size, mismatches, db, proteinsPositions, SEQ_SIZE, clusters);
		
	}

	clusterSize = clusters.size();
	in_cluster2 = new bool[clusters.size()];
	posibleConnections = new bool*[clusterSize];
	for (i = 0; i < clusterSize; i++) {
		posibleConnections[i] = new bool[clusterSize]();
		pair<unsigned int, unsigned int> hashPair = make_pair(clusters[i].getSeed(), i);

		seedHash.begin();
		seedHash.insert(hashPair);
		in_cluster2[i] = true;
		temp = clusters[i].getPositions();
		
		seedHash.end();
		for (j = 1; j < temp.size(); j++) {
			pair<unsigned int, unsigned int> hashPair = make_pair(temp[j], i);
			seedHash.begin();
			seedHash.insert(hashPair);
			seedHash.end();


		}
	}

cout << "Preprocess2 finished" << endl;
}







void InitializePreprocess3(int mismatches, int form_size)
{
	unsigned int clusterSize, i;
	pair<unsigned int, unsigned int> hashPair;
	//Currently only form sizes of 4 and 5 are implemented
	if (form_size != 4 && form_size != 5 && form_size != 6)
	{
		cout << "Form size is invalid" << endl;
		cin.get();
		exit(1);
	}

	/*Initialize ASCII hash, both functions are doing exactly the same thing
	but they are divided between 2 classes and are used by inline functions.*/
	preprocessClass.InitializeHash();
	findSeqClass.InitializeHashFind();


	//Check if there is already a file that contain the preprocessed form types results
	cout << "Preprocess3 started" << endl;
	//Check if a form types set file exist
	bool preprocessed = preprocessClass.CheckIfPreProcessed(form_size, mismatches, SEQ_SIZE);
	if (!preprocessed)
	{
		cout << "Form types file is missing" << endl;
		cin.get();
		exit(1);
	}
	else
	{
		if (form_size == 4)
			preprocessClass.PreProcess53(forms5, formTypes, form_size, mismatches, db, proteinsPositions, SEQ_SIZE, clusters);

	}

	clusterSize = clusters.size();
	Connections = new bool*[clusterSize];
	seedHash.clear();
	for (i = 0; i < clusterSize; i++) {
		Connections[i] = new bool[clusterSize]();
		pair<unsigned int, unsigned int> hashPair = make_pair(clusters[i].getSeed(), i);
		seedHash.begin();
		seedHash.insert(hashPair);
		seedHash.end();
	}

	cout << "Preprocess3 finished" << endl;
}




void GetSequences(string sequencesPath, string sequencesFileName, int mismatches, int form_size, vector<ClusterClass> &clusters)
{
	double algorithmRunTime, naiveRunTime;
	unsigned int i, threadIndex;
	int to_print = 0, pos_print = 0;

	for (i = 0; i < db.size() - 20; i++) {
		sequences.push_back(db.substr(i, 20));
	}
	unsigned int startSequence, endSequence;
	
	unsigned int numOfSeq = sequences.size();

											
	for (i = 0; i < numOfSeq; i++)
	{
		startSequence = i;
		if (!in_cluster[startSequence]) {

			ClusterClass nCluster(startSequence);
			clusters.push_back(nCluster);
			in_cluster[startSequence] = true;
			stringstream ss;
			ss << startSequence << " ";
			startFindAlgorithm(ss.str());
		}
	}

	for (to_print = 0; to_print < clusters.size(); to_print++) {
		vector<unsigned int> tmp = clusters[to_print].getPositions();
		printf_s("cluster number  -  %d \n", clusters[to_print].getSeed());
		for (vector<unsigned int>::iterator it = tmp.begin(); it != tmp.end(); it++, pos_print++) {
			cout << *it;
			printf_s(" ");

		}
		printf_s("\n");
	}
}


void startFindAlgorithm(string parameters)
{
	unordered_set<unsigned int> positions;
	vector<SequenceResults> sequencesResults;
	stringstream ss;
	unsigned int startSequence, endSequence, i;
	ss << parameters;
	ss >> startSequence;
	i = startSequence;

	/*A function that finds the sequence in the database, the function will
	write to the positions vector the positions were there is a match*/
	
		if (form_size == 4)
			findSeqClass.FindSequence(db, positions, proteinsIndexes, forms4, formTypes, sequences[i], form_size, mismatches, SEQ_SIZE, startSequence, clusters);
		else if (form_size == 5)
			findSeqClass.FindSequence5(db, positions, proteinsIndexes, forms5, formTypes, sequences[i], form_size, mismatches, SEQ_SIZE, startSequence, clusters);
		else if (form_size == 6)
			findSeqClass.FindSequence6(db, positions, proteinsIndexes, forms6, formTypes, sequences[i], form_size, mismatches, SEQ_SIZE, startSequence, clusters);
		if (positions.size() > 0)
		{
			SequenceResults newResults;
			newResults.sequence = sequences[i];
			newResults.positions = positions;
			sequencesResults.push_back(newResults);
		}
		positions.clear();

}



void PossibleClusters40(string sequencesPath, string sequencesFileName, int mismatches, int form_size, vector<ClusterClass> &clusters)
{
	unsigned int i;
	int to_print = 0, pos_print = 0;
	pair<unsigned int, string> tempPair, tempPair2;

	for (i = 0; i < clusters.size(); i++) {

		pair<unsigned int, string> tempPair = make_pair(i, db.substr(clusters[i].getSeed(), 20));
	
		seedSequences.insert(tempPair);
	
		printf("seq num %d \n", seedSequences.size());
	}
	unsigned int startSequence, endSequence;
	
	unsigned int numOfSeq = seedSequences.size();


	for (i = 0; i < clusters.size() ; i++)
	{
		startSequence = i;			
					
		auto search = seedHash.find(i);
		if (in_cluster2[search->second])continue;
		in_cluster2[search->second] = true;
		stringstream ss;
		ss << startSequence << " ";
		PossibleMatchAlgorithm(ss.str());
	}


	for (to_print = 0; to_print < clusters.size(); to_print++) {
		vector<unsigned int> tmp = clusters[to_print].getPosibilities();
		printf_s("possible connections of cluster number  -  %d \n", clusters[to_print].getSeed());
		for (vector<unsigned int>::iterator it = tmp.begin(); it != tmp.end(); it++, pos_print++) {
			cout << *it;
			printf_s(" ");

		}
		printf_s("\n");

	}
}

	void PossibleMatchAlgorithm(string parameters)
	{
		unordered_set<unsigned int> positions;
		vector<SequenceResults> sequencesResults;
		stringstream ss;
		string searchSeq, searchSeq2;

		unsigned int startSequence, endSequence, i;
		ss << parameters;
		ss >> startSequence;
		i = startSequence;
		auto search = seedSequences.find(i);
		searchSeq = search->second;
		if (form_size == 4) {
			findSeqClass.FindPossibilities(db, positions, proteinsIndexes, forms4, formTypes, searchSeq, form_size, mismatches, SEQ_SIZE, startSequence, clusters);
			}
		if (positions.size() > 0)
		{
			SequenceResults newResults;
			newResults.sequence = sequences[i];
			newResults.positions = positions;
			sequencesResults.push_back(newResults);
		}
		positions.clear();
	
	}


	void clustersConnection(string sequencesPath, string sequencesFileName, int mismatches, int form_size, vector<ClusterClass> &clusters)
	{
		unsigned int i,j;
		int to_print = 0, pos_print = 0;
		vector<unsigned int> dists;
		pair<unsigned int, string> tempPair, tempPair2;

		for (i = 0; i < clusters.size(); i++) {
			dists = clusters[i].getPositions();
			for (j = 0; j < dists.size(); j++) {
				pair<unsigned int, string> tempPair = make_pair(i, db.substr(dists[j], 20));
				
				seedSequences2.insert(tempPair);
			
				printf("seq num %d \n", seedSequences2.size());
				if (seedSequences2.size() == 659)
					putchar('s');
			}
		}
		unsigned int startSequence, endSequence;
	
		unsigned int numOfSeq = seedSequences2.size();


		for (i = 0; i < clusters.size(); i++)
		{
			startSequence = i;
			stringstream ss;
			ss << startSequence << " ";
			connectionsAlgorithm(ss.str());
		}


		for (to_print = 0; to_print < clusters.size(); to_print++) {
			vector<unsigned int> tmp = clusters[to_print].getConnections();
			printf_s("real connections of cluster number  -  %d \n", clusters[to_print].getSeed());
			for (vector<unsigned int>::iterator it = tmp.begin(); it != tmp.end(); it++, pos_print++) {
				cout << *it;
				printf_s(" ");

			}
			printf_s("\n");
			
		}
	}



	void connectionsAlgorithm(string parameters)
	{
		unordered_set<unsigned int> positions;
		vector<SequenceResults> sequencesResults;
		stringstream ss;
		string searchSeq, searchSeq2;	
		unsigned int startSequence, endSequence, i;
		ss << parameters;
		ss >> startSequence;
		i = startSequence;
		auto search = seedSequences.find(i);
		searchSeq = search->second;
		
		/*A function that finds the sequence in the database, the function will
		write to the positions vector the positions were there is a match*/
	
		if (form_size == 4) {
			findSeqClass.FindPossibilities(db, positions, proteinsIndexes, forms4, formTypes, searchSeq, form_size, mismatches, SEQ_SIZE, startSequence, clusters);
			
		}
			else if (form_size == 5)
		findSeqClass.FindConnections(db, positions, proteinsIndexes, forms5, formTypes, sequences[i], form_size, mismatches, SEQ_SIZE, startSequence, clusters);
		
		if (positions.size() > 0)
		{
			SequenceResults newResults;
			newResults.sequence = sequences[i];
			newResults.positions = positions;
			sequencesResults.push_back(newResults);
		}
		positions.clear();
	
	}







int main(int argc, char *argv[])
{
	string DBFileName, DBpath, sequencesFileName, sequencesPath;
	int flag = 0, break_flag = 1;
	
	/*If the program is executed like this for example:
	file.exe db.txt 6 2
	.*/

	//Comment the following if statement if you test the program in Visual studio.
	
	if (argc <= 1)
	{
		cout << "For help use -help when running the application:\nfile.exe -help" << endl;
		return 0;
	}
	if (argc >= 2 && (argv[1])[0] == '-' && (argv[1])[1] == 'h' && (argv[1])[2] == 'e' && (argv[1])[3] == 'l' && (argv[1])[4] == 'p')
	{
		cout << "1) If you want to search for individual sequences in DB:\nOn the command line, run the exe file like this:\nfile.exe [db file] [mismatches] [form_size]\nFor example : file.exe C:\\db.txt 8 4\nPress q to exit the program\n\n2) If you want to search a list of sequences in DB:\nOn the command line, run the exe file like this:\nfile.exe [db file] [mismatches] [form_size]\nFor example: file.exe -l C:\\db.txt C:\\SequencesList.txt 8 4" << endl;
		return 0;
	}
	else if (argc >= 2 && (argv[1])[0] == '-' && (argv[1])[1] == 'l')
	{
		flag = 2;
		DBFileName = argv[flag];
		sequencesFileName = argv[1 + flag];
	}

	else if (argc >= 2)
		DBFileName = argv[1 + flag];
	else
		DBFileName = "C:\\Users\\ProBook\\ProteinDB4.txt";

	//Divide the path from the filename
	for (int i = DBFileName.size(); i >= 0 && break_flag; i--)
	{
		if (DBFileName[i] == '\\')
		{
			DBpath = DBFileName.substr(0, i) + '\\';
			DBFileName = DBFileName.substr(i + 1, DBFileName.size());
			break_flag = 0;
		}
	}
	if (flag == 2 && (argv[1])[0] == '-' && (argv[1])[1] == 'l')
	{
		for (int i = sequencesFileName.size(); i >= 0; i--)
		{
			if (sequencesFileName[i] == '\\')
			{
				sequencesPath = sequencesFileName.substr(0, i) + '\\';
				sequencesFileName = sequencesFileName.substr(i + 1, sequencesFileName.size());
			}
		}
	}

	
	ReadDB(DBpath, DBFileName);

							   //Get the number of mismatches and form size from the command line or set them to 8 and 4 as default
	mismatches = (argc >= 3) ? atoi(argv[2 + flag]) : 8;
	form_size = (argc >= 4) ? atoi(argv[3 + flag]) : 4;

	in_cluster = new bool[db.size()]();
	

	InitializePreprocess(mismatches, form_size);

	
		//Search for a list of sequences
		GetSequences(sequencesPath, sequencesFileName, mismatches, form_size, clusters);
		system("pause");


		// **********up until here weve done preprocess-I********


		
		forms6.clear();
		formTypes.clear();

		//********preprocess-II*********

		form_size = 4;
		mismatches = 12;
		
		InitializePreprocess2(mismatches, form_size);
		

		PossibleClusters40(sequencesPath, sequencesFileName, 12, 4, clusters);
		system("pause");
		forms4.clear();

		form_size = 5;
		mismatches = 8;

		InitializePreprocess3(mismatches, form_size);


		clustersConnection(sequencesPath, sequencesFileName, mismatches, form_size, clusters);
		system("pause");
	cout << "End of program" << endl;

	if (argc <= 1)
		cin.get();
	return 0;
}