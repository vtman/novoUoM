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

#include <omp.h>

const int64_t LONG_LENGTH = 1000000;
const int64_t SHORT_LENGTH = 10000;
const int64_t LEN_D = 18;

class SortTreeData {
public:
	SortTreeData();
	~SortTreeData();

	int openRootFile();
	int getKeys();
	int countEntries();
	int getBranches();
	int radixSort(char* vSortBufferMain, char* vSortBuffer0, char* vSortBuffer1, int n);

	int processNew();
	int startProcessing();
	int copyBinaryFile(const char* source, const char* destination);

	char inputFileName[1000], ioFileName[1000], prefix[1000], outputFolder[1000];
	char inFileName[1000], outFileName[1000];

	
	TFile *inFile;
	TList* listOfKeys;
	TKey* key;
	TTree* tree;
	TBranch* branch;

	char* treeNames, * branchNames;
	char* vTypeName;

	char** pvBuffer;

	int nThreads;

	int* vStatus, * vSame;
	int* vIndDTkeys, * ivTemp;
	int* vTypeSize, * vArraySize;
	int nKeys, nDTkeys, nBranches;
	
	Long64_t* vNumEntries;
};


int SortTreeData::copyBinaryFile(const char* source, const char* destination) {
	char* buf;
	long long nBuf;
	size_t bytes;
	FILE* fi, * fo;

	fi = fopen(source, "rb");
	if (fi == nullptr) {
		printf("Error: cannot open file %s\n", source);
		return -1;
	}

	fo = fopen(destination, "wb");
	if (fo == nullptr) {
		printf("Error: cannot open file %s\n", destination);
		return -2;
	}

	nBuf = LEN_D * LONG_LENGTH;
	buf = pvBuffer[0];

	while ((bytes = fread(buf, 1, nBuf, fi)) > 0) {
		fwrite(buf, 1, bytes, fo);
	}

	fclose(fi); fi = nullptr;
	fclose(fo); fo = nullptr;

	return 0;
}

SortTreeData::SortTreeData() {
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

	pvBuffer = nullptr;

	nThreads = 0;
}

SortTreeData::~SortTreeData() {
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
	if (pvBuffer != nullptr) {
		for (int i = 0; i < 3*nThreads; i++) {
			if (pvBuffer[i] != nullptr) {free(pvBuffer[i]); pvBuffer[i] = nullptr;}
		}
		free(pvBuffer); pvBuffer = nullptr;
	}
	
}

int SortTreeData::openRootFile() {
	printf("Input File: %s\n", inputFileName);
	inFile = TFile::Open(inputFileName);
	if (!inFile || inFile->IsZombie()) {
		printf("Error: cannot open file %s\n", inputFileName);
		return -1;
	}
	printf("The file is openned\n");
	return 0;
}

int SortTreeData::radixSort(char * vSortBufferMain, char *vSortBuffer0, char * vSortBuffer1, int n) {
	int n0, n1, nbits;
	long long Tmin, Tmax, T, dT;

	Tmin = *(long long*)(vSortBufferMain);
	Tmax = Tmin;

	for (int i = 1; i < n; i++) {
		T = *(long long*)(vSortBufferMain + LEN_D * i);
		if (T < Tmin) Tmin = T;
		if (T > Tmax) Tmax = T;
	}
	
	dT = Tmax - Tmin;
	nbits = 1;
	while (dT > 0) {
		dT /= 2;
		nbits++;
	}

	for (int i = 0; i < n; i++) {
		T = *(long long*)(vSortBufferMain + LEN_D * i) - Tmin;
		*(long long*)(vSortBufferMain + LEN_D * i) = T;
	}

	for (int i = 0; i < 4; i++) {
		n0 = 0;
		n1 = 0;
		for (int k = 0; k < n; k++) {
			if ((( (int)(*(vSortBufferMain + k * LEN_D + 16)) >> i) & 1) == 0) {
				memcpy(vSortBuffer0 + n0 * LEN_D, vSortBufferMain + k * LEN_D, LEN_D);
				n0++;
			}
			else {
				memcpy(vSortBuffer1 + n1 * LEN_D, vSortBufferMain + k * LEN_D, LEN_D);
				n1++;
			}
		}
		if (n0 == 0 || n1 == 0)continue;
		memcpy(vSortBufferMain, vSortBuffer0, n0 * LEN_D);
		memcpy(vSortBufferMain + n0 * LEN_D, vSortBuffer1, n1 * LEN_D);
	}

	for (int i = 0; i < nbits; i++) {
		n0 = 0;
		n1 = 0;
		for (int k = 0; k < n; k++) {
			if ((((*(Long64_t*)(vSortBufferMain + k * LEN_D)) >> i) & 1) == 0) {
				memcpy(vSortBuffer0 + n0 * LEN_D, vSortBufferMain + k * LEN_D, LEN_D);
				n0++;
			}
			else {
				memcpy(vSortBuffer1 + n1 * LEN_D, vSortBufferMain + k * LEN_D, LEN_D);
				n1++;
			}
		}
		if (n0 == 0 || n1 == 0)continue;
		memcpy(vSortBufferMain, vSortBuffer0, n0 * LEN_D);
		memcpy(vSortBufferMain + n0 * LEN_D, vSortBuffer1, n1 * LEN_D);
	}

	for (int i = 0; i < n; i++) {
		T = *(long long*)(vSortBufferMain + LEN_D * i) + Tmin;
		*(long long*)(vSortBufferMain + LEN_D * i) = T;
	}

	return 0;
}


int SortTreeData::processNew() {
	Long64_t nTot;
	FILE *fio;
	char boardName[100];

	long long maxRowsInChunk, M;
	long long* vStart, * vEnd;
	int nBlocks, iPercOld, iPerc, nSmallRows, icb;
	int* vDone;

	maxRowsInChunk = LONG_LENGTH / LEN_D;
	nSmallRows = SHORT_LENGTH / LEN_D;

	pvBuffer = (char**)malloc(sizeof(char*) * 3* nThreads);

	for (int i = 0; i < 3 * nThreads; i++) {
		pvBuffer[i] = nullptr;
		pvBuffer[i] = (char*)malloc(sizeof(char) * LEN_D * LONG_LENGTH);
		if (pvBuffer[i] == nullptr) {
			printf("Error: cannot allocate memory\n");
			return -1;
		}
	}

	for (int it = 0; it < nDTkeys; it++) {
		key = (TKey*)listOfKeys->At(vIndDTkeys[it]);
		sprintf(boardName, "%s", key->GetName());
		tree = (TTree*)inFile->Get(boardName);
		nTot = tree->GetEntries();
		printf("Tree: %i (%s), %lli\n", it, boardName, nTot);
		
		icb = boardName[2] - 'a';
		sprintf(inFileName, "%s/%s_tree_%i.bin", outputFolder, prefix, icb);
		sprintf(outFileName, "%s/%s_tree_sorted_%i.bin", outputFolder, prefix, icb);
		if (copyBinaryFile(inFileName, outFileName) != 0) return -1;

		printf("File #%i: %s\n", it + 1, outFileName);
		fio = fopen(outFileName, "rb+");
		if (fio == nullptr) {
			printf("Error: cannot open output file %s\n", outFileName);
			return -2;
		}

		M = vNumEntries[it];

		nBlocks = M / maxRowsInChunk;
		if (M % maxRowsInChunk > 0) nBlocks++;
		printf("Number of blocks: %i\n", nBlocks);

		vStart = (long long*)malloc(sizeof(long long) * nBlocks);
		vEnd = (long long*)malloc(sizeof(long long) * nBlocks);

		vDone = (int*)malloc(sizeof(int) * nBlocks);

		for (int i = 0; i < nBlocks; i++) {
			vDone[i] = 0;
			vStart[i] = maxRowsInChunk * (long long)(i);
			vEnd[i] = maxRowsInChunk * (long long)(i + 1);
		}
		if (vEnd[nBlocks - 1] > M) vEnd[nBlocks - 1] = M;

		iPercOld = -1;

#pragma omp parallel
		{
			int tid;
			bool isContinue;
			long long nSize, nRows;
			char* locBuf, *locBuf0, *locBuf1;
			tid = omp_get_thread_num();
			locBuf = pvBuffer[3 * tid];
			locBuf0 = pvBuffer[3 * tid + 1];
			locBuf1 = pvBuffer[3 * tid + 2];
			for (int ib = 0; ib < nBlocks; ib++) {
#pragma omp critical
				{
					isContinue = true;
					if (vDone[ib] == 1) {
						isContinue = false;
					}
					else {
						vDone[ib] = 1;
						nRows = vEnd[ib] - vStart[ib];
						nSize = nRows * LEN_D;
						#ifdef _WIN32
						_fseeki64(fio, vStart[ib] * LEN_D, SEEK_SET);
						#else
						fseeko64(fio, vStart[ib] * LEN_D, SEEK_SET);
						#endif
						fread(locBuf, sizeof(char), nSize, fio);

						iPerc = (int)(double(vStart[ib]) * 100.0 / double(M));
						if (iPerc != iPercOld) {
							printf("Long %i %%\n", iPerc);
							iPercOld = iPerc;
						}
					}
				}
				if (!isContinue) continue;

				radixSort(locBuf, locBuf0, locBuf1, nRows);

#pragma omp critical
				{
					#ifdef _WIN32
					_fseeki64(fio, vStart[ib] * LEN_D, SEEK_SET);
					#else
					fseeko64(fio, vStart[ib] * LEN_D, SEEK_SET);
					#endif
					fwrite(locBuf, sizeof(char), nSize, fio);
				}
			}
		}

		for (int i = 0; i < nBlocks - 1; i++) {
			vDone[i] = 0;
			vStart[i] = maxRowsInChunk * (long long)(i + 1) - nSmallRows;
			vEnd[i] = vStart[i] + 2 * nSmallRows;
		}

		vDone[nBlocks - 1] = 0;
		vStart[nBlocks - 1] = M - nSmallRows;
		vEnd[nBlocks - 1] = M;

		if (vEnd[nBlocks - 2] >= vStart[nBlocks - 1]) {
			vEnd[nBlocks - 2] = M;
			nBlocks--;
		}

		iPercOld = -1;

#pragma omp parallel
		{
			int tid;
			bool isContinue;
			long long nSize, nRows;
			char* locBuf, * locBuf0, * locBuf1;
			tid = omp_get_thread_num();
			locBuf = pvBuffer[3 * tid];
			locBuf0 = pvBuffer[3 * tid + 1];
			locBuf1 = pvBuffer[3 * tid + 2];
			for (int ib = 0; ib < nBlocks; ib++) {
#pragma omp critical
				{
					isContinue = true;
					if (vDone[ib] == 1) {
						isContinue = false;
					}
					else {
						vDone[ib] = 1;
						nRows = vEnd[ib] - vStart[ib];
						nSize = nRows * LEN_D;
						#ifdef _WIN32
						_fseeki64(fio, vStart[ib] * LEN_D, SEEK_SET);
						#else
						fseeko64(fio, vStart[ib] * LEN_D, SEEK_SET);
						#endif
						fread(locBuf, sizeof(char), nSize, fio);
						iPerc = (int)(double(vStart[ib]) * 100.0 / double(M));
						if (iPerc != iPercOld) {
							printf("Short %i %%\n", iPerc);
							iPercOld = iPerc;
						}
					}
				}
				if (!isContinue) continue;

				radixSort(locBuf, locBuf0, locBuf1, nRows);

#pragma omp critical
				{
					#ifdef _WIN32
					_fseeki64(fio, vStart[ib] * LEN_D, SEEK_SET);
					#else
					fseeko64(fio, vStart[ib] * LEN_D, SEEK_SET);
					#endif
					fwrite(locBuf, sizeof(char), nSize, fio);
				}
			}
		}
		fclose(fio); fio = nullptr;
		

		free(vDone); vDone = nullptr;
		free(vStart); vStart = nullptr;
		free(vEnd); vEnd = nullptr;

		printf("\n");
		
	}
	return 0;
}

int SortTreeData::getKeys() {
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

int SortTreeData::countEntries() {
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

int SortTreeData::getBranches() {
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

int SortTreeData::startProcessing() {
	if (openRootFile() != 0) return -1;
	if (getKeys() != 0) return -2;
	if (countEntries() != 0) return -3;
	if (getBranches() != 0) return -4;

#pragma omp parallel
	{
		nThreads = omp_get_num_threads();
	}

	printf("\nNumber of threads: %i\n", nThreads);

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
	SortTreeData* sd;
	
	if (argc != 4) {
		printUsage();
		return -99;
	}

	auto start = std::chrono::high_resolution_clock::now();

	sd = new SortTreeData();

	sprintf(sd->inputFileName, "%s", argv[1]);
	sprintf(sd->outputFolder, "%s", argv[2]);
	sprintf(sd->prefix, "%s", argv[3]);
	
	iResult = sd->startProcessing();
	delete sd; sd = nullptr;

	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsed = finish - start;
	printf("Elapsed Time: %10.2f s", 0.001 * elapsed.count());

	return iResult;
}
