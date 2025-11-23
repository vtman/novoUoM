//D:\NOVO\2025_02_Cf252\in\det_000153.root D:\NOVO\2025_02_Cf252\test\temp\153_coinTwo.bin D:\NOVO\2025_02_Cf252\test\temp\153_twoHits.root

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TBasket.h>
#include <TKey.h>
#include <TLeaf.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include <ctime>

const long long LEN_D = 22;
const int MAX_TREES = 16;
const int N_PREVIOUS = 20;

class Bin2Root {
public:
	Bin2Root();
	~Bin2Root();

	int getBinaryData();
	int openRootFile();
	int getKeys();
	int countEntries();
	int initOutputFile();
	int startProcessing();

	char inputRootName[1000], outputRootName[1000], binaryFileName[1000];
	
	FILE* fiBin;

	TFile *inFile, *outFile;
	TList* listOfKeys;
	TKey* key;
	TTree* tree, *outTree, **vTrees;
	TBranch* branch;

	char* vBinaryData, * treeNames;

	int* vActiveTrees, *vReverse;
	int* vStatus, * vSame, * vIndDTkeys;

	int nKeys, nDTkeys, nActiveTrees, nTrees;

	long long * vActiveSelectionCounter, * vTreeCounter, **pvSelection;
	long long fileSize;

	Long64_t* vNumEntries;
};

Bin2Root::Bin2Root() {
	vReverse = nullptr;
	pvSelection = nullptr;
	vActiveSelectionCounter = nullptr;
	vBinaryData = nullptr;
	inFile = nullptr;
	outFile = nullptr;
	vIndDTkeys = nullptr;
	vNumEntries = nullptr;
	treeNames = nullptr;
	vStatus = nullptr;
	vSame = nullptr;

	vTrees = nullptr;

	nTrees = 0;

	vActiveTrees = nullptr;
	vTreeCounter = nullptr;

	fiBin = nullptr;
}

Bin2Root::~Bin2Root() {
	if (vReverse != nullptr) {free(vReverse); vReverse = nullptr;}
	if (vActiveSelectionCounter != nullptr) {
		free(vActiveSelectionCounter); vActiveSelectionCounter = nullptr;
	}
	if (pvSelection != nullptr) {
		for (int i = 0; i < nActiveTrees; i++) {
			if (pvSelection[i] != nullptr) {free(pvSelection[i]); pvSelection[i] = nullptr;}
		}
		free(pvSelection); pvSelection = nullptr;
	}
	if (vBinaryData != nullptr) {free(vBinaryData); vBinaryData = nullptr;}
	if (fiBin != nullptr) {fclose(fiBin); fiBin = nullptr;}
	if (vActiveTrees != nullptr) {free(vActiveTrees); vActiveTrees = nullptr;}
	if (vTreeCounter != nullptr) {free(vTreeCounter); vTreeCounter = nullptr;}
	if (vTrees != nullptr) {free(vTrees); vTrees = nullptr;}
	if (inFile != nullptr) {inFile->Close(); inFile = nullptr;}
	if (outFile != nullptr) {outFile->Close(); outFile = nullptr;}
	if (vIndDTkeys != nullptr) {free(vIndDTkeys); vIndDTkeys = nullptr;}
	if (vNumEntries != nullptr) {free(vNumEntries); vNumEntries = nullptr;}
	if (treeNames != nullptr) {free(treeNames); treeNames = nullptr;}
	if (vStatus != nullptr) {free(vStatus); vStatus = nullptr;}
	if (vSame != nullptr) {free(vSame); vSame = nullptr;}
}

int Bin2Root::getBinaryData() {
	long long nRecords, pos;
	long long iPrevious, iCurrent;
	bool isFound;
	int iB;

	fiBin = fopen(binaryFileName, "rb");
	printf("Binary file: %s\n", binaryFileName);
	if (fiBin == nullptr) {
		printf("Error: cannot open the binary file %s\n", binaryFileName);
		return -1;
	}

	vActiveTrees = (int*)malloc(sizeof(int) * MAX_TREES);
	vReverse = (int*)malloc(sizeof(int) * MAX_TREES);

	for (int i = 0; i < MAX_TREES; i++) {
		vReverse[i] = -1;
	}

	vTreeCounter = (long long*)malloc(sizeof(long long) * MAX_TREES);
	for (int i = 0; i < MAX_TREES; i++) {
		vTreeCounter[i] = 0;
	}

#ifdef _WIN32
	_fseeki64(fiBin, 0, SEEK_END);
	fileSize = _ftelli64(fiBin);
	_fseeki64(fiBin, 0, SEEK_SET);
#else
	fseeko64(fiBin, 0, SEEK_END);
	fileSize = ftello64(fiBin);
	fseeko64(fiBin, 0, SEEK_SET);
#endif

	printf("File size: %lli\n", fileSize);

	vBinaryData = (char*)malloc(sizeof(char) * fileSize);
	if (vBinaryData == nullptr) {
		printf("Error: memory for the binary data\n");
		return -2;
	}

	fread(vBinaryData, sizeof(char), fileSize, fiBin);
	fclose(fiBin); fiBin = nullptr;

	nRecords = fileSize / LEN_D;
	printf("Number of records: %lli\n", nRecords);

	for (long long i = 0; i < nRecords; i++) {
		iB = vBinaryData[i * LEN_D + 16] >> 4;
		vTreeCounter[iB]++;
	}

	printf("\nSelected entries:\n");
	nActiveTrees = 0;
	for (int i = 0; i < MAX_TREES; i++) {
		if (vTreeCounter[i] == 0) continue;
		vActiveTrees[nActiveTrees] = i;
		vReverse[i] = nActiveTrees;
		printf("%i) dt%c: %lli\n", nActiveTrees + 1, 'a' + i, vTreeCounter[i]);
		nActiveTrees++;
	}
	printf("\n");

	if (nActiveTrees == 0) return 0;
	vActiveSelectionCounter = (long long*)malloc(sizeof(long long) * nActiveTrees);

	nActiveTrees = 0;
	for (int i = 0; i < MAX_TREES; i++) {
		if (vTreeCounter[i] == 0) continue;
		vActiveSelectionCounter[nActiveTrees] = vTreeCounter[i];
		nActiveTrees++;
	}

	pvSelection = (long long**)malloc(sizeof(long long*) * nActiveTrees);

	for (int i = 0; i < nActiveTrees; i++) {
		pvSelection[i] = (long long*)malloc(sizeof(long long) * vActiveSelectionCounter[i]);
		vActiveSelectionCounter[i] = 0;
	}

	for (long long i = 0; i < nRecords; i++) {
		iB = vReverse[vBinaryData[i * LEN_D + 16] >> 4];
		pos = *(long long*)(vBinaryData + i * LEN_D + 8);
		isFound = false;
		iCurrent = vActiveSelectionCounter[iB];
		iPrevious = iCurrent - N_PREVIOUS;
		if (iPrevious < 0) iPrevious = 0;
		for (long long j = iPrevious; j < iCurrent; j++) {
			if (pvSelection[iB][j] == pos) {
				isFound = true;
				break;
			}
		}
		if (!isFound) {
			pvSelection[iB][vActiveSelectionCounter[iB]] = pos;
			vActiveSelectionCounter[iB]++;
		}
	}

	free(vBinaryData); vBinaryData = nullptr;

	return 0;
}

int Bin2Root::openRootFile() {
	printf("Input File: %s\n", inputRootName);
	inFile = TFile::Open(inputRootName);
	if (!inFile || inFile->IsZombie()) {
		printf("Error: cannot open file %s\n", inputRootName);
		return -1;
	}
	printf("The file is openned\n");
	return 0;
}

int Bin2Root::getKeys() {
	int nSame, indCycle, nBest;
	listOfKeys = inFile->GetListOfKeys();
	nKeys = listOfKeys->GetEntries();
	printf("Number of keys: %i\n", nKeys);

	if (nKeys == 0) return -1;

	nDTkeys = 0;

	vStatus = (int*)malloc(sizeof(int) * nKeys);
	vSame = (int*)malloc(sizeof(int) * nKeys);

	nDTkeys = 0;

	for (int i = 0; i < nKeys; i++) {
		vStatus[i] = 0;
		key = (TKey *)listOfKeys->At(i);
		if (strcmp(key->GetClassName(), "TTree") != 0) continue;
		if (strncmp("dt", key->GetName(), 2) != 0) continue;
		if (strlen(key->GetName()) != 3) continue;
		vStatus[i] = 1;
		nDTkeys++;
	}

	printf("Number of keys (TTree, dt...): %i\n", nDTkeys);

	nDTkeys = 0;
	for (int i = 0; i < nKeys; i++) {
		if (vStatus[i] == 0) continue;
		nSame = 0;
		vSame[nSame] = i;
		nSame++;
		for (int j = i + 1; j < nKeys; j++) {
			if (vStatus[j] == 0) continue;
			if (strcmp(listOfKeys->At(j)->GetName(), listOfKeys->At(i)->GetName()) != 0) continue;
			vSame[nSame] = j;
			nSame++;
		}
		if (nSame == 1) {
			nDTkeys++;
			continue;
		}
		indCycle = -1;
		nBest = -1;
		for (int j = 0; j < nSame; j++) {
			key = (TKey*)listOfKeys->At(vSame[j]);
			if (key->GetCycle() > indCycle) {
				indCycle = key->GetCycle();
				nBest = j;
			}
		}

		for (int j = 0; j < nSame; j++) {
			if (j != nBest) vStatus[vSame[j]] = 0;
		}
		nDTkeys++;
	}

	free(vSame); vSame = nullptr;

	printf("Number of keys (TTree, dt..., single cycle): %i\n", nDTkeys);
		
	vIndDTkeys = (int*)malloc(sizeof(int) * nDTkeys);
	nDTkeys = 0;

	nTrees = 0;

	for (int i = 0; i < nKeys; i++) {
		if (vStatus[i] == 0)continue;
		key = (TKey*)listOfKeys->At(i);
		tree = (TTree*)inFile->Get(key->GetName());
		vIndDTkeys[nDTkeys] = i;
		nTrees++;
		nDTkeys++;
	}

	vTrees = (TTree**)malloc(sizeof(TTree*) * nTrees);

	nTrees = 0;
	for (int i = 0; i < nKeys; i++) {
		if (vStatus[i] == 0)continue;
		key = (TKey*)listOfKeys->At(i);
		vTrees[nTrees] = (TTree*)inFile->Get(key->GetName());
		nTrees++;
	}

	free(vStatus); vStatus = nullptr;

	if (nDTkeys == 0) return -1;

	return 0;
}

int Bin2Root::countEntries() {
	vNumEntries = (Long64_t*)malloc(sizeof(Long64_t) * nDTkeys);
	treeNames = (char*)malloc(sizeof(char) * 50 * nDTkeys);

	for (int i = 0; i < nDTkeys; i++) {
		key = (TKey*)listOfKeys->At(vIndDTkeys[i]);
		sprintf(treeNames + 50 * i, "%s", key->GetName());
		tree = (TTree*)inFile->Get(key->GetName());
		vNumEntries[i] = tree->GetEntries();
		printf("%i) %s\t%lli\n", i + 1, treeNames + 50 * i, vNumEntries[i]);
	}
	printf("\n");

	return 0;
}

int Bin2Root::initOutputFile() {
	char locTreeName[1000];
	int indActive, indTree;
	long long nC, * vSelection;

	outFile = TFile::Open(outputRootName, "RECREATE");
	if (outFile == nullptr) {
		printf("Error: cannot create output file %s\n", outputRootName);
		return -1;
	}
	if (outFile->IsZombie()) {
		printf("Error: cannot create output file %s\n", outputRootName);
		return -2;
	}

	for (int i = 0; i < nTrees; i++) {
		outTree = vTrees[i]->CloneTree(0);
		sprintf(locTreeName, "%s", vTrees[i]->GetName());
		indTree = locTreeName[2] - 'a';
		indActive = vReverse[indTree];
		vSelection = pvSelection[indActive];
		nC = vActiveSelectionCounter[indActive];
		vTrees[i]->SetBranchStatus("*", 1);

		for (long long j = 0; j < nC; j++) {
			vTrees[i]->GetEntry(vSelection[j]);
			outTree->Fill();
		}
		outTree->Write(locTreeName);

		delete outTree; outTree = nullptr;
	}

	outFile->Close();

	return 0;
}


int Bin2Root::startProcessing() {
	if (getBinaryData() != 0) return -1;
	if (openRootFile() != 0) return -2;
	if (getKeys() != 0) return -3;
	if (countEntries() != 0) return -4;
	if (initOutputFile() != 0) return -5;
	
	return 0;
}

int printUsage() {
	printf("Input parameters:\n");
	printf("\t1) Input ROOT file (string)\n");
	printf("\t2) Binary file-manager (string)\n");
	printf("\t3) Output ROOT file (string)\n\n");
	return 0;
}

int main(int argc, char** argv) {
	int iResult;
	Bin2Root* br;
	
	if (argc != 4) {
		printUsage();
		return -99;
	}

/*	TFile* tf = TFile::Open("D:\\NOVO\\2025_02_Cf252\\test\\temp\\153_twoHits.root", "READ");
	TTree* tree = (TTree*)tf->Get("dta");
	tree->Scan("*", "", "", 22);
	printf("\n");
	TTree* tree2 = (TTree*)tf->Get("dtb");
	tree2->Scan("*", "", "", 22);
	return 0;
	*/

	auto start = std::chrono::high_resolution_clock::now();

	br = new Bin2Root();

	sprintf(br->inputRootName, "%s", argv[1]);
	sprintf(br->binaryFileName, "%s", argv[2]);
	sprintf(br->outputRootName, "%s", argv[3]);
	
	iResult = br->startProcessing();
	delete br; br = nullptr;

	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsed = finish - start;
	printf("Elapsed Time: %10.2f s", 0.001 * elapsed.count());

	return iResult;
}
