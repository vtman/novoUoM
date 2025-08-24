//D:\NOVO\conData\in\det_000206.root D:\NOVO\conData\out\ts_000206_cbs.bin D:\NOVO\conData\out\ts_000206_out.root

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

const int BUFF_LENGTH = 1000000;


class CoinRoot {
public:
	CoinRoot();
	~CoinRoot();

	int openRootFile();
	int getKeys();
	int countEntries();
	int getBranches();
	int initOutputFile();

	int mergeTrees();
	int startProcessing();

	char inputFileName[1000], outputFileName[1000], manFileName[1000];
	char dtNames[100];

	
	TFile *inFile, *outFile;
	TList* listOfKeys;
	TKey* key;
	TTree* tree, *outTree, **vTrees;
	TBranch* branch;

	int nTrees;

	char* treeNames, * branchNames;
	char* vTypeName, * vTempBuffer;
	char* vControl;

	int* vStatus, * vSame;
	int* vIndDTkeys, * ivTemp;
	int* vTypeSize, * vArraySize;
	int nKeys, nDTkeys, lenD, nBranches;

	int bitsPerEntry;
	int* offsetBranch;
	char* localVariables;
	
	Long64_t* vNumEntries;

	char *ccc;
};

CoinRoot::CoinRoot() {
	vControl = nullptr;
	inFile = nullptr;
	outFile = nullptr;
	ivTemp = nullptr;
	vIndDTkeys = nullptr;
	vNumEntries = nullptr;
	treeNames = nullptr;
	vStatus = nullptr;
	vSame = nullptr;
	branchNames = nullptr;
	offsetBranch = nullptr;
	localVariables = nullptr;

	vTypeSize = nullptr;
	vArraySize = nullptr;
	vTypeName = nullptr;

	vTrees = nullptr;

	lenD = 18;
	nTrees = 0;

	vTempBuffer = nullptr;
}

CoinRoot::~CoinRoot() {
	if (vControl != nullptr) {
		free(vControl); vControl = nullptr;
	}
	if (vTrees != nullptr) {
		free(vTrees); vTrees = nullptr;
	}
	if (offsetBranch != nullptr) {
		free(offsetBranch); offsetBranch = nullptr;
	}
	if (localVariables != nullptr) {
		free(localVariables); localVariables = nullptr;
	}
	if (inFile != nullptr) {
		inFile->Close(); inFile = nullptr;
	}
	if (outFile != nullptr) {
		outFile->Close(); outFile = nullptr;
	}

	if (vIndDTkeys != nullptr) {
		free(vIndDTkeys); vIndDTkeys = nullptr;
	}
	if (ivTemp != nullptr) {
		free(ivTemp); ivTemp = nullptr;
	}
	if (vNumEntries != nullptr) {
		free(vNumEntries); vNumEntries = nullptr;
	}
	if (treeNames != nullptr) {
		free(treeNames); treeNames = nullptr;
	}
	if (vStatus != nullptr) {
		free(vStatus); vStatus = nullptr;
	}
	if (vSame != nullptr) {
		free(vSame); vSame = nullptr;
	}
	if (branchNames != nullptr) {
		free(branchNames); branchNames = nullptr;
	}
	if (vTypeSize != nullptr) {
		free(vTypeSize); vTypeSize = nullptr;
	}
	if (vArraySize != nullptr) {
		free(vArraySize); vArraySize = nullptr;
	}
	if (vTypeName != nullptr) {
		free(vTypeName); vTypeName = nullptr;
	}
	if (vTempBuffer != nullptr) {
		free(vTempBuffer); vTempBuffer = nullptr;
	}
}

int CoinRoot::openRootFile() {
	printf("Input File: %s\n", inputFileName);
	inFile = TFile::Open(inputFileName);
	if (!inFile || inFile->IsZombie()) {
		printf("Error: cannot open file %s\n", inputFileName);
		return -1;
	}
	printf("The file is openned\n");
	return 0;
}


int CoinRoot::getKeys() {
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

int CoinRoot::countEntries() {
	vNumEntries = (Long64_t*)malloc(sizeof(Long64_t) * nDTkeys);
	treeNames = (char*)malloc(sizeof(char) * 50 * nDTkeys);

	for (int i = 0; i < nDTkeys; i++) {
		key = (TKey*)listOfKeys->At(vIndDTkeys[i]);
		sprintf(treeNames + 50 * i, "%s", key->GetName());
		dtNames[i] = treeNames[50 * i + 2];
		tree = (TTree*)inFile->Get(key->GetName());
		vNumEntries[i] = tree->GetEntries();
		printf("%i) %s\t%lli\n", i + 1, treeNames + 50 * i, vNumEntries[i]);
	}
	printf("\n");

	return 0;
}

int CoinRoot::getBranches() {
	int n, ntot;
	TBranch* locBranch;
	TLeaf* leaf;
	TObjArray* allBranches;
	char* tempNames, * tempType;
	int* vtempSize, * vtempArray, * vStatus;

	key = (TKey*)listOfKeys->At(vIndDTkeys[0]);
	tree = (TTree*)inFile->Get(key->GetName());

	allBranches = tree->GetListOfBranches();

	n = allBranches->GetEntries();
	printf("Number of branches (initial): %i\n", n);
	if (n == 0) return -1;
	
	tempNames = (char*)malloc(sizeof(char) * 100 * n);

	vtempSize = (int*)malloc(sizeof(int) * n);
	vtempArray = (int*)malloc(sizeof(int) * n);
	tempType = (char*)malloc(sizeof(char) * n * 20);

	ntot = 0;

	for (Int_t i = 0; i < n; i++) {
		locBranch = (TBranch*)allBranches->At(i);
		leaf = locBranch->GetLeaf(locBranch->GetName());

		if (!leaf) {
			printf("Warning: cannot get leaf for branch %s\n", locBranch->GetName());
			continue;
		}

		vtempArray[ntot] = leaf->GetNdata();
		sprintf(tempType + 20 * ntot, "%s", leaf->GetTypeName());
		vtempSize[ntot] = leaf->GetLenType();

		sprintf(tempNames + 100 * ntot, "%s", locBranch->GetName());

		ntot++;
	}
	printf("Number of branches (med): %i\n", ntot);
	if (ntot == 0) return 0;

	vStatus = (int*)malloc(sizeof(int) * ntot);

	for (int i = 0; i < ntot; i++) {
		vStatus[i] = 1;
	}

	for (int k = 0; k < nDTkeys; k++) {
		key = (TKey*)listOfKeys->At(vIndDTkeys[k]);
		tree = (TTree*)inFile->Get(key->GetName());

		for (int i = 0; i < ntot; i++) {
			locBranch = tree->GetBranch(tempNames + 100 * i);
			if (locBranch == nullptr) {
				vStatus[i] = 0;
				continue;
			}
			leaf = locBranch->GetLeaf(tempNames + 100 * i);

			if (vtempArray[i] != leaf->GetNdata()) {
				vStatus[i] = 0;
				continue;
			}
			if (strcmp(tempType + 20 * i, leaf->GetTypeName()) != 0) {
				vStatus[i] = 0;
				continue;
			}
			if (vtempSize[i] != leaf->GetLenType()) {
				vStatus[i] = 0;
				continue;
			}
		}
	}
	nBranches = 0;
	for (int i = 0; i < ntot; i++) {
		if (vStatus[i] == 1) nBranches++;
	}

	printf("Number of branches (final): %i\n", nBranches);

	vTypeSize = (int*)malloc(sizeof(int) * nBranches);
	vArraySize = (int*)malloc(sizeof(int) * nBranches);
	vTypeName = (char*)malloc(sizeof(char) * nBranches * 20);
	branchNames = (char*)malloc(sizeof(char) * nBranches * 100);


	nBranches = 0;
	for (int i = 0; i < ntot; i++) {
		if (vStatus[i] == 0) continue;
		vTypeSize[nBranches] = vtempSize[i];
		vArraySize[nBranches] = vtempArray[i];
		strncpy(branchNames + 100 * nBranches, tempNames + i * 100, 100);
		strncpy(vTypeName + 20 * nBranches, tempType + 20 * i, 20);
		nBranches++;
	}

	free(vtempSize); vtempSize = nullptr;
	free(vtempArray); vtempArray = nullptr;
	free(tempType); tempType = nullptr;
	free(tempNames); tempNames = nullptr;

	free(vStatus); vStatus = nullptr;

	for (int i = 0; i < nBranches; i++) {
		printf("%i) %s\t%s\t%i\t%i\n", i + 1, branchNames + 100 * i, vTypeName + 20 * i, vTypeSize[i], vArraySize[i]);
	}

	printf("\n");

	return 0;
}

int CoinRoot::initOutputFile() {
	TTree* templateTree;
	TBranch* locBranch;
	TLeaf* leaf;
	char branchName[100], leafList[100];

	outFile = TFile::Open(outputFileName, "RECREATE");
	if (outFile == nullptr) {
		printf("Error: cannot create output file %s\n", outputFileName);
		return -1;
	}
	if (outFile->IsZombie()) {
		printf("Error: cannot create output file %s\n", outputFileName);
		return -2;
	}
	
	outTree = new TTree("dtAll", "selected entries");

	templateTree = vTrees[0];


	offsetBranch = (int*)malloc(sizeof(int) * (nBranches + 1));
	offsetBranch[0] = 0;

	bitsPerEntry = 0;

	for (int i = 0; i < nBranches; i++) {
		locBranch = templateTree->GetBranch(branchNames + 100 * i);
		leaf = locBranch->GetLeaf(branchNames + 100 * i);
		offsetBranch[i+ 1] = offsetBranch[i] + (leaf->GetNdata() * leaf->GetLenType());
		printf("%i\t%i\n", i, offsetBranch[i + 1]);
	}

	bitsPerEntry = offsetBranch[nBranches];

	localVariables = (char*)malloc(sizeof(char) * (bitsPerEntry + 1));

	ccc = localVariables + bitsPerEntry;

	for (int i = 0; i < nBranches; i++) {
		locBranch = templateTree->GetBranch(branchNames + 100 * i);
		
		strcpy(branchName, locBranch->GetName());
		strcpy(leafList, locBranch->GetTitle());
		for (int j = 0; j < nTrees; j++) {
			vTrees[j]->SetMakeClass(1);
			vTrees[j]->SetBranchAddress(branchName, localVariables + offsetBranch[i]);
		}
		outTree->Branch(branchName, localVariables + offsetBranch[i], locBranch->GetTitle());
	}
	
	outTree->Branch("board", ccc, "board/b");

	return 0;
}

int CoinRoot::mergeTrees() {
	long long fileSize, nRead, n2r, m2r;
	Long64_t nCEntries, indEntry;
	FILE* fc;
	int indTree;
	char* vloc;
	vControl = (char*)malloc(sizeof(char) * lenD * BUFF_LENGTH);
	if (vControl == nullptr) {
		printf("Error: memory\n");
		return -1;
	}


	fc = fopen(manFileName, "rb");
	if (fc == nullptr) {
		printf("Error: cannot open file %s\n", manFileName);
		return -2;
	}

#ifdef _WIN32
	_fseeki64(fc, 0, SEEK_END);
	fileSize = _ftelli64(fc);
	_fseeki64(fc, 0, SEEK_SET);
	#else
	fseeko64(fc, 0, SEEK_END);
	fileSize = ftello64(fc);
	fseeko64(fc, 0, SEEK_SET);
	#endif
	
	nCEntries = fileSize / ((long long)lenD);
	n2r = BUFF_LENGTH / lenD;
	m2r = n2r * lenD;

	printf("Number of output entries: %lli\n", nCEntries);




	do {
		nRead = fread(vControl, sizeof(char), m2r, fc);
		if (nRead <= 0) break;
		nRead /= lenD;

		for (long long i = 0; i < nRead; i++) {
			vloc = vControl + i * lenD;
			indTree = (int)(*(vloc + 17));
			indEntry = *(long long*)(vloc + 8);
			if (indTree < 0 || indTree >= nTrees)continue;
			if (indEntry < 0 || indEntry >= vNumEntries[indTree])continue;
			*ccc = dtNames[indTree];
			vTrees[indTree]->GetEntry(indEntry);
			outTree->Fill();
		}
	} while (true);

	

	fclose(fc); fc = nullptr;

	outTree->Write();

	return 0;
}

int CoinRoot::startProcessing() {
	if (openRootFile() != 0) return -1;
	if (getKeys() != 0) return -2;
	if (countEntries() != 0) return -3;
	if (getBranches() != 0) return -4;
	if (initOutputFile() != 0) return -5;
	if (mergeTrees() != 0) return -6;
	
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
	CoinRoot* se;
	
	if (argc != 4) {
		printUsage();
		return -99;
	}

	auto start = std::chrono::high_resolution_clock::now();

	se = new CoinRoot();

	sprintf(se->inputFileName, "%s", argv[1]);
	sprintf(se->manFileName, "%s", argv[2]);
	sprintf(se->outputFileName, "%s", argv[3]);
	
	iResult = se->startProcessing();
	delete se; se = nullptr;

	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsed = finish - start;
	printf("Elapsed Time: %10.2f s", 0.001 * elapsed.count());

	return iResult;
}
