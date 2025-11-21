//D:\NOVO\conData\in\det_000206.root D:\NOVO\conData\out test

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TKey.h>
#include <TLeaf.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <chrono>

const int64_t BUFF_LENGTH = 1000000;
const int nFTbits = 10ULL;
const int nScale = (1ULL << nFTbits);
const uint64_t maskFT = ((1ULL << nFTbits) - 1ULL);
const float fScale = float(nScale);
const int64_t LEN_D = 22;

class ExtractTime {
public:
	ExtractTime();
	~ExtractTime();

	int openRootFile();
	int getKeys();
	int countEntries();
	int getBranches();
	int processNew();
	int startProcessing();

	char inputFileName[1000], outputFileName[1000], prefix[1000], outputFolder[1000];
	
	TFile *inFile;
	TList* listOfKeys;
	TKey* key;
	TTree* tree;
	TBranch* branch;

	char* treeNames, * branchNames, * vTypeName, * vTempBuffer;

	int* vStatus, * vSame, * vIndDTkeys, * ivTemp, * vTypeSize, * vArraySize;
	int nKeys, nDTkeys, nBranches;
	
	Long64_t* vNumEntries;
};

ExtractTime::ExtractTime() {
	inFile = nullptr;
	ivTemp = nullptr;
	vIndDTkeys = nullptr;
	vNumEntries = nullptr;
	treeNames = nullptr;
	vStatus = nullptr;
	vSame = nullptr;
	branchNames = nullptr;

	vTypeSize = nullptr;
	vArraySize = nullptr;
	vTypeName = nullptr;

	vTempBuffer = nullptr;
}

ExtractTime::~ExtractTime() {
	if (inFile != nullptr) {inFile->Close(); inFile = nullptr;}
	if (vIndDTkeys != nullptr) {free(vIndDTkeys); vIndDTkeys = nullptr;}
	if (ivTemp != nullptr) {free(ivTemp); ivTemp = nullptr;}
	if (vNumEntries != nullptr) {free(vNumEntries); vNumEntries = nullptr;}
	if (treeNames != nullptr) {free(treeNames); treeNames = nullptr;}
	if (vStatus != nullptr) {free(vStatus); vStatus = nullptr;}
	if (vSame != nullptr) {free(vSame); vSame = nullptr;}
	if (branchNames != nullptr) {free(branchNames); branchNames = nullptr;}
	if (vTypeSize != nullptr) {free(vTypeSize); vTypeSize = nullptr;}
	if (vArraySize != nullptr) {free(vArraySize); vArraySize = nullptr;}
	if (vTypeName != nullptr) {free(vTypeName); vTypeName = nullptr;}
	if (vTempBuffer != nullptr) {free(vTempBuffer); vTempBuffer = nullptr;}
}

int ExtractTime::openRootFile() {
	printf("Input File: %s\n", inputFileName);
	inFile = TFile::Open(inputFileName);
	if (!inFile || inFile->IsZombie()) {
		printf("Error: cannot open file %s\n", inputFileName);
		return -1;
	}
	printf("The file is openned\n");
	return 0;
}

void correct3t(char* vc, uint64_t utemp) {
	uint64_t u1, u2, u3;
	u1 = *(uint64_t*)vc;
	u2 = utemp & maskFT;
	u3 = (u1 << nFTbits) | u2;
	*(uint64_t*)vc = u3;
}

int ExtractTime::processNew() {
	Long64_t nTot, nLeft, ij;
	uint16_t uBase, uMin, uMax, uU, uBase2;
	char cb, locBuffer[20], boardName[100], B;
	int icb, sw;
	bool isCurrent, isPrevious;
	FILE* fo;

	vTempBuffer = (char*)malloc(sizeof(char) * LEN_D * BUFF_LENGTH);

	//timestamp (8)
	//position (8)
	//CB (1)
	//Eshort (2)
	//Elong (2)
	//flags (1)

	int ipCurrent, ipNew;

	int maxBoards;

	maxBoards = -1;
	

	for (int it = 0; it < nDTkeys; it++) {
		key = (TKey*)listOfKeys->At(vIndDTkeys[it]);
		tree = (TTree*)inFile->Get(key->GetName());
		nTot = tree->GetEntries();
		sprintf(boardName, "%s", key->GetName());
		printf("Tree: %i (%s), %lli\n", it, boardName, nTot);
		icb = boardName[2] - 'a';

		if (icb > maxBoards) maxBoards = icb;

		sprintf(outputFileName, "%s/%s_tree_%i.bin", outputFolder, prefix, icb);
		fo = fopen(outputFileName, "wb");
		if (fo == nullptr) {
			printf("Error: cannot open output file %s\n", outputFileName);
			return -1;
		}

		ipCurrent = -1;

		cb = icb;
		B = (icb << 4);

		memset(vTempBuffer, 0, sizeof(char) * LEN_D * BUFF_LENGTH);
		for (int k = 0; k < BUFF_LENGTH; k++) {
			*(vTempBuffer + LEN_D * (k + 1) - 1) = cb;
		}

		uBase = 0;
		isPrevious = false;

		for (Long64_t i = 0; i < nTot; i += BUFF_LENGTH) {
			ipNew = (int)floor(double(i) * 100.0 / double(nTot));

			if (ipNew != ipCurrent) {
				printf("\t%i%%\t%lli\n", ipNew, i);
				ipCurrent = ipNew;
			}
			
			nLeft = nTot - i;
			if (nLeft > BUFF_LENGTH) nLeft = BUFF_LENGTH;

			for (Long64_t j = 0; j < nLeft; j++) {
				ij = i + j;
				*(Long64_t *)(vTempBuffer + LEN_D * j + 8) = ij;
				*(int32_t*)(vTempBuffer + LEN_D * j + 4) = 0;
			}

			//timestamp
			branch = tree->GetBranch("timestamp");
			branch->SetAddress(locBuffer);
			sw = 4;
			for (Long64_t j = 0; j < nLeft; j++) {
				ij = i + j;
				branch->GetEntry(ij);
				memcpy(vTempBuffer + j * LEN_D + 0, locBuffer, sw * sizeof(char));
			}

			//timestampExtended
			branch = tree->GetBranch("timestampExtended");
			branch->SetAddress(locBuffer);
			sw = 2;

			ij = i;
			branch->GetEntry(ij);
			memcpy(vTempBuffer + 4, locBuffer, sw * sizeof(char));
			uU = *(uint16_t*)locBuffer;
			uMin = uU;
			uMax = uU;

			for (Long64_t j = 0; j < nLeft; j++) {
				ij = i + j;
				branch->GetEntry(ij);
				memcpy(vTempBuffer + j * LEN_D + 4, locBuffer, sw * sizeof(char));
				uU = *(uint16_t*)locBuffer;
				if (uU < uMin) uMin = uU;
				if (uU > uMax) uMax = uU;
			}

			if (uMax - uMin > 30000) {
				isCurrent = true;
			}
			else {
				isCurrent = false;
			}
			if (isCurrent) {
				uBase2 = uBase + 1;
				for (Long64_t j = 0; j < nLeft; j++) {
					uU = *(uint16_t*)(vTempBuffer + j * LEN_D + 4);

					if (uU > 30000) {
						*(uint16_t*)(vTempBuffer + j * LEN_D + 6) = uBase;
					}
					else {
						*(uint16_t*)(vTempBuffer + j * LEN_D + 6) = uBase2;
					}
				}
			}else if (isPrevious) {
				uBase++;
				for (Long64_t j = 0; j < nLeft; j++) {
					*(uint16_t*)(vTempBuffer + j * LEN_D + 6) = uBase;
				}
			}
			else {
				for (Long64_t j = 0; j < nLeft; j++) {
					*(uint16_t*)(vTempBuffer + j * LEN_D + 6) = uBase;
				}
			}
			isPrevious = isCurrent;

			//time
			branch = tree->GetBranch("time");
			branch->SetAddress(locBuffer);
			sw = 4;
			uint64_t utemp;
			for (Long64_t j = 0; j < nLeft; j++) {
				ij = i + j;
				branch->GetEntry(ij);
				utemp = (uint64_t)(0.5f + *(float*)(locBuffer)*fScale);
				correct3t(vTempBuffer + j * LEN_D, utemp);
			}

			//channel
			branch = tree->GetBranch("channel");
			branch->SetAddress(locBuffer);
			sw = 1;
			for (Long64_t j = 0; j < nLeft; j++) {
				ij = i + j;
				branch->GetEntry(ij);
				locBuffer[0] |= B;
				memcpy(vTempBuffer + j * LEN_D + 16, locBuffer, sw * sizeof(char));
			}

			//Eshort
			branch = tree->GetBranch("Eshort");
			branch->SetAddress(locBuffer);
			sw = 2;
			for (Long64_t j = 0; j < nLeft; j++) {
				ij = i + j;
				branch->GetEntry(ij);
				memcpy(vTempBuffer + j * LEN_D + 17, locBuffer, sw * sizeof(char));
			}

			//Elong
			branch = tree->GetBranch("Elong");
			branch->SetAddress(locBuffer);
			sw = 2;
			for (Long64_t j = 0; j < nLeft; j++) {
				ij = i + j;
				branch->GetEntry(ij);
				memcpy(vTempBuffer + j * LEN_D + 19, locBuffer, sw * sizeof(char));
			}

			//flags
			branch = tree->GetBranch("flags");
			branch->SetAddress(locBuffer);
			sw = 1;
			for (Long64_t j = 0; j < nLeft; j++) {
				ij = i + j;
				branch->GetEntry(ij);
				memcpy(vTempBuffer + j * LEN_D + 20, locBuffer, sw * sizeof(char));
			}

			fwrite(vTempBuffer, sizeof(char), LEN_D * nLeft, fo);
			
		}
		printf("Done\n\n");
		fclose(fo); fo = nullptr;
	}

	maxBoards++;
	sprintf(outputFileName, "%s/%s_maxBoards.txt", outputFolder, prefix);
	fo = fopen(outputFileName, "w");
	if (fo == nullptr) {
		printf("Error: cannot open output file %s\n", outputFileName);
		return -1;
	}
	fprintf(fo, "%i\n", maxBoards);
	fclose(fo); fo = nullptr;

	return 0;
}

int ExtractTime::getKeys() {
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
	for (int i = 0; i < nKeys; i++) {
		if (vStatus[i] == 0)continue;
		vIndDTkeys[nDTkeys] = i;
		nDTkeys++;
	}

	free(vStatus); vStatus = nullptr;

	if (nDTkeys == 0) return -1;

	return 0;
}

int ExtractTime::countEntries() {
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

int ExtractTime::getBranches() {
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

int ExtractTime::startProcessing() {
	if (openRootFile() != 0) return -1;
	if (getKeys() != 0) return -2;
	if (countEntries() != 0) return -3;
	if (getBranches() != 0) return -4;
	if (processNew() != 0) return -5;
	
	return 0;
}

int printUsage() {
	printf("Input parameters:\n");
	printf("\t1) Input ROOT file (string)\n");
	printf("\t2) Output folder (string)\n");
	printf("\t3) prefix (string)\n\n");
	return 0;
}

int main(int argc, char** argv) {
	int iResult;
	ExtractTime* se;
	
	if (argc != 4) {
		printUsage();
		return -99;
	}

	auto start = std::chrono::high_resolution_clock::now();

	se = new ExtractTime();

	sprintf(se->inputFileName, "%s", argv[1]);
	sprintf(se->outputFolder, "%s", argv[2]);
	sprintf(se->prefix, "%s", argv[3]);
	
	iResult = se->startProcessing();
	delete se; se = nullptr;

	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsed = finish - start;
	printf("Elapsed Time: %10.2f s", 0.001 * elapsed.count());

	return iResult;
}
