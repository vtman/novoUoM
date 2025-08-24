//D:\NOVO\conData\in\det_000206.root D:\NOVO\conData\out test

#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <stdio.h>
#include <fstream>
#include <chrono>

#include <omp.h>

const int BUFF_LENGTH = 100000;
const int CHUNK_LENGTH = 10000;

class MergeSort {
public:
	MergeSort();
	~MergeSort();

	int openRootFile();
	int getKeys();
	int startProcessing();
	int allocateMemory();
	int getChunks();
	int prepareMT();
	int runMT();

	char inputFileName[1000], outputFileName[1000], outputFolder[1000], prefix[1000];
	char ioFileName[1000];
	
	TFile *inFile;
	TList* listOfKeys;
	TKey* key;
	TTree* tree;

	FILE** vFI, *fo;

	int* vStatus, * vSame, *ivTemp;
	int nKeys, lenD, nLists, nChunks, nShift, nThreads;

	char** pvBufferInput, ** pvBufferOutput;
	char* vTempBuffer, * treeNames;
	long long** pvIndex;

	long long* vIndStartChunk, * vIndStartOutput, scMax;
	long long** pvStartInput;
	
	Long64_t* vNumberOfElements;
	Long64_t nElements;
};

MergeSort::MergeSort() {
	inFile = nullptr;
	ivTemp = nullptr;
	treeNames = nullptr;
	vStatus = nullptr;
	vSame = nullptr;

	vIndStartChunk = nullptr;
	vIndStartOutput = nullptr;

	vFI = nullptr;
	fo = nullptr;

	pvBufferInput = nullptr;
	pvBufferOutput = nullptr;
	pvStartInput = nullptr;
	pvIndex = nullptr;

	vNumberOfElements = nullptr;
	vTempBuffer = nullptr;

	nLists = 0;
	lenD = 18;
	nThreads = 0;

}

MergeSort::~MergeSort() {
	if (pvIndex != nullptr) {
		for (int i = 0; i < nThreads; i++) {
			if (pvIndex[i] != nullptr) {
				free(pvIndex[i]); pvIndex[i] = nullptr;
			}
		}
		free(pvIndex); pvIndex = nullptr;
	}

	if (pvBufferInput != nullptr) {
		for (int i = 0; i < nThreads; i++) {
			if (pvBufferInput[i] != nullptr) {
				free(pvBufferInput[i]); pvBufferInput[i] = nullptr;
			}
		}
		free(pvBufferInput); pvBufferInput = nullptr;
	}

	if (pvBufferOutput != nullptr) {
		for (int i = 0; i < nThreads; i++) {
			if (pvBufferOutput[i] != nullptr) {
				free(pvBufferOutput[i]); pvBufferOutput[i] = nullptr;
			}
		}
		free(pvBufferOutput); pvBufferOutput = nullptr;
	}

	if (pvStartInput != nullptr) {
		for (int i = 0; i < nThreads; i++) {
			if (pvStartInput[i] != nullptr) {
				free(pvStartInput[i]); pvStartInput[i] = nullptr;
			}
		}
		free(pvStartInput); pvStartInput = nullptr;
	}

	if (vIndStartChunk != nullptr) {
		free(vIndStartChunk); vIndStartChunk = nullptr;
	}
	if (vIndStartOutput != nullptr) {
		free(vIndStartOutput); vIndStartOutput = nullptr;
	}
	if (inFile != nullptr) {
		inFile->Close(); inFile = nullptr;
	}
	if (fo != nullptr) {
		fclose(fo); fo = nullptr;
	}
	if (ivTemp != nullptr) {
		free(ivTemp); ivTemp = nullptr;
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
	if (vNumberOfElements != nullptr) {
		free(vNumberOfElements); vNumberOfElements = nullptr;
	}
	if (vTempBuffer != nullptr) {
		free(vTempBuffer); vTempBuffer = nullptr;
	}
	if (vFI != nullptr) {
		for (int i = 0; i < nLists; i++) {
			if (vFI[i] != nullptr) {
				fclose(vFI[i]); vFI[i] = nullptr;
			}
		}
		free(vFI); vFI = nullptr;
	}
}

int MergeSort::openRootFile() {
	printf("Input File: %s\n", inputFileName);
	inFile = TFile::Open(inputFileName);
	if (!inFile || inFile->IsZombie()) {
		printf("Error: cannot open file %s\n", inputFileName);
		return -1;
	}
	printf("The file is openned\n");
	return 0;
}

int MergeSort::allocateMemory() {
	vNumberOfElements = (Long64_t*)malloc(sizeof(Long64_t) * nLists);
	if (vNumberOfElements == nullptr) return -1;

	nElements = 0;

	sprintf(ioFileName, "%s/%s_merged.bin", outputFolder, prefix);
	fo = fopen(ioFileName, "wb");
	if (fo == nullptr) {
		printf("Error: cannot open file %s\n", ioFileName);
		return -7;
	}

	vFI = (FILE**)malloc(sizeof(FILE*) * nLists);
	for (int i = 0; i < nLists; i++) {
		sprintf(ioFileName, "%s/%s_tree_sorted_%i.bin", outputFolder, prefix, i);
		vFI[i] = fopen(ioFileName, "rb");
		if (vFI[i] == nullptr) {
			printf("Error: cannot open file %s\n", ioFileName);
			return -6;
		}
		#ifdef _WIN32
		_fseeki64(vFI[i], 0, SEEK_END);
		vNumberOfElements[i] = _ftelli64(vFI[i]) / (long long)(lenD);
		_fseeki64(vFI[i], 0, SEEK_SET);
		#else
		fseeko64(vFI[i], 0, SEEK_END);
		vNumberOfElements[i] = ftello64(vFI[i]) / (long long)(lenD);
		fseeko64(vFI[i], 0, SEEK_SET);
		#endif
	}

	for (int i = 0; i < nLists; i++) {
		nElements += vNumberOfElements[i];
	}
	printf("Total number of entries: %lli\n", nElements);

	vTempBuffer = (char*)malloc(sizeof(char) * lenD * BUFF_LENGTH);

	return 0;
}

int MergeSort::getKeys() {
	int nSame, indCycle, nBest;
	listOfKeys = inFile->GetListOfKeys();
	nKeys = listOfKeys->GetEntries();
	printf("Number of keys: %i\n", nKeys);

	if (nKeys == 0) return -1;

	vStatus = (int*)malloc(sizeof(int) * nKeys);
	vSame = (int*)malloc(sizeof(int) * nKeys);

	nLists = 0;

	for (int i = 0; i < nKeys; i++) {
		vStatus[i] = 0;
		key = (TKey *)listOfKeys->At(i);
		if (strcmp(key->GetClassName(), "TTree") != 0) continue;
		if (strncmp("dt", key->GetName(), 2) != 0) continue;
		vStatus[i] = 1;
		nLists++;
	}

	printf("Number of keys (TTree, dt...): %i\n", nLists);

	nLists = 0;
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
			nLists++;
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
		nLists++;
	}
	free(vSame); vSame = nullptr;

	printf("Number of keys (TTree, dt..., single cycle): %i\n", nLists);

	nLists = 0;
	for (int i = 0; i < nKeys; i++) {
		if (vStatus[i] == 0)continue;
		nLists++;
	}

	free(vStatus); vStatus = nullptr;
	return 0;
}

int getBits(long long A) {
	int n;
	long long B;
	if (A == 0) return 1;
	B = A;
	n = 0;
	while (B > 0) {
		n++;
		B = B >> 1;
	}
	return n;
}

int MergeSort::getChunks() {
	long long maxElements, maxValue, iValue;
	long long nChunks1;
	int nBitsMax;
	long long indPrev, indCurrent, n, nRead, ij;
	long long s1, scLoc;

	maxElements = 0;
	maxValue = 0;
	for (int i = 0; i < nLists; i++) {
	        #ifdef _WIN32
		_fseeki64(vFI[i], lenD * (vNumberOfElements[i] - 1), SEEK_SET);
		fread(&iValue, sizeof(long long), 1, vFI[i]);
		_fseeki64(vFI[i], 0, SEEK_SET);
		#else
		fseeko64(vFI[i], lenD * (vNumberOfElements[i] - 1), SEEK_SET);
		fread(&iValue, sizeof(long long), 1, vFI[i]);
		fseeko64(vFI[i], 0, SEEK_SET);
		#endif
		printf("%lli\n", iValue);
		if (iValue > maxValue) maxValue = iValue;
		if (vNumberOfElements[i] > maxElements) maxElements = vNumberOfElements[i];
	}
	printf("Maximum number of elements: %lli\n", maxElements);
	printf("Maximum timestamp value: %lli\n", maxValue);

	nChunks1 = maxElements / CHUNK_LENGTH;
	printf("Number of chunks (estimated): %lli\n", nChunks1);

	nBitsMax = getBits(maxValue);

	nShift = nBitsMax - getBits(nChunks1);
	printf("Shift: %i\n", nShift);

	nChunks = (maxValue >> nShift) + 2;
	printf("Number of chunks (final): %i\n", nChunks);

	vIndStartChunk = (long long*)malloc(sizeof(long long) * nChunks * nLists);
	vIndStartOutput = (long long*)malloc(sizeof(long long) * nChunks);

	for (int i = 0; i < nLists; i++) {
		indPrev = -1;
		n = vNumberOfElements[i];
		#ifdef _WIN32
		_fseeki64(vFI[i], 0, SEEK_SET);
		#else
		fseeko64(vFI[i], 0, SEEK_SET);
		#endif
		ij = 0;
		for (long long j = 0; j < n; j += BUFF_LENGTH) {
			nRead = n - j;
			if (nRead > BUFF_LENGTH) nRead = BUFF_LENGTH;
			fread(vTempBuffer, sizeof(char), lenD * nRead, vFI[i]);

			for (long long k = 0; k < nRead; k++) {
				indCurrent = (*(long long*)(vTempBuffer + k * lenD)) >> nShift;
				if (indCurrent == indPrev)continue;
				ij = j + k;
				for (long long m = indPrev + 1; m <= indCurrent; m++) {
					vIndStartChunk[m * nLists + i] = ij;
				}
				indPrev = indCurrent;
			}
		}
		for (long long m = indPrev + 1; m < nChunks; m++) {
			vIndStartChunk[m * nLists + i] = vNumberOfElements[i];
		}
	}

	scMax = 0;

	for (int j = 0; j < nChunks; j++) {
		s1 = 0;
		for (int i = 0; i < nLists; i++) {
			s1 += vIndStartChunk[j * nLists + i];
		}
		vIndStartOutput[j] = s1;
	}
	for (int i = 0; i < nChunks - 1; i++) {
		scLoc = vIndStartOutput[i + 1] - vIndStartOutput[i];
		if (scLoc > scMax) scMax = scLoc;
	}

	printf("Maximum number of elements in output chuncks: %lli\n", scMax);
	return 0;
}

int MergeSort::prepareMT() {
#pragma omp parallel
	{
		nThreads = omp_get_num_threads();
	}

	printf("\nNumber of threads: %i\n", nThreads);

	pvBufferInput = (char**)malloc(sizeof(char*) * nThreads);
	pvBufferOutput = (char**)malloc(sizeof(char*) * nThreads);
	pvIndex = (long long**)malloc(sizeof(long long*) * nThreads);

	for (int i = 0; i < nThreads; i++) {
		pvBufferInput[i] = nullptr;
		pvBufferInput[i] = (char*)malloc(sizeof(char) * lenD * scMax);
		if (pvBufferInput[i] == nullptr) {
			printf("Error: buffer\n");
			return -1;
		}
	}
	for (int i = 0; i < nThreads; i++) {
		pvBufferOutput[i] = nullptr;
		pvBufferOutput[i] = (char*)malloc(sizeof(char) * lenD * scMax);
		if (pvBufferOutput[i] == nullptr) {
			printf("Error: buffer\n");
			return -2;
		}
	}

	for (int i = 0; i < nThreads; i++) {
		pvIndex[i] = nullptr;
		pvIndex[i] = (long long*)malloc(sizeof(long long) * nLists);
		if (pvIndex[i] == nullptr) {
			printf("Error: buffer\n");
			return -2;
		}
	}

	pvStartInput = (long long**)malloc(sizeof(long long*) * nThreads);

	for (int i = 0; i < nThreads; i++) {
		pvStartInput[i] = (long long*)malloc(sizeof(long long) * (nLists + 1));
	}

	for (int i = 0; i < nLists; i++) {
	#ifdef _WIN32
		_fseeki64(vFI[i], 0, SEEK_SET);
		#else
		fseeko64(vFI[i], 0, SEEK_SET);
		#endif
	}

	return 0;
}

int mSort(char* vIn, char* vOut, long long* vStart, long long *vPos, int n, long long lenD) {
	long long M, valBest, val;
	int indBest;
	M = vStart[n];

	for (int i = 0; i < n; i++) {
		vPos[i] = vStart[i];
	}
	for (long long i = 0; i < M; i++) {
		indBest = -1;
		valBest = 0x7FFFFFFFFFFFFFFF;
		for (int j = 0; j < n; j++) {
			if (vPos[j] >= vStart[j + 1])continue;
			val = *(long long*)(vIn + lenD * vPos[j]);
			if (val >= valBest) continue;
			valBest = val;
			indBest = j;
		}
		memcpy(vOut + lenD * i, vIn + vPos[indBest] * lenD, lenD);
		vPos[indBest]++;
	}
	return 0;
}

int MergeSort::runMT() {
	int* vDone;
	int indPerc, indPercOld;
	vDone = (int*)malloc(sizeof(int) * (nChunks - 1));

	indPercOld = -1;

	for (int i = 0; i < nChunks - 1; i++) {
		vDone[i] = 0;
	}

#pragma omp parallel
	{
		int tid;
		bool isContinue;
		char* vLocIn, * vLocOut;
		long long* vLocStart, *vLocIndex;
		long long nRows;

		tid = omp_get_thread_num();
		vLocIn = pvBufferInput[tid];
		vLocOut = pvBufferOutput[tid];
		vLocStart = pvStartInput[tid];
		vLocIndex = pvIndex[tid];
		for (int ic = 0; ic < nChunks - 1; ic++) {
#pragma omp critical
			{
				isContinue = true;
				if (vDone[ic] == 1) {
					isContinue = false;
				}
				else {
					vDone[ic] = 1;
					vLocStart[0] = 0;
					for (int k = 0; k < nLists; k++) {
						nRows = vIndStartChunk[(ic + 1) * nLists + k] - vIndStartChunk[ic * nLists + k];
						vLocStart[k + 1] = vLocStart[k] + nRows;
						#ifdef _WIN32
						_fseeki64(vFI[k], lenD * vIndStartChunk[ic * nLists + k], SEEK_SET);
						#else
						fseeko64(vFI[k], lenD * vIndStartChunk[ic * nLists + k], SEEK_SET);
						#endif
						fread(vLocIn + vLocStart[k]*lenD, sizeof(char), lenD * nRows, vFI[k]);
					}
				}
			}
			if (!isContinue) continue;

			mSort(vLocIn, vLocOut, vLocStart, vLocIndex, nLists, lenD);

#pragma omp critical
			{
				#ifdef _WIN32
				_fseeki64(fo, vIndStartOutput[ic] * lenD, SEEK_SET);
				#else
				fseeko64(fo, vIndStartOutput[ic] * lenD, SEEK_SET);
				#endif
				fwrite(vLocOut, sizeof(char), lenD* vLocStart[nLists], fo);
				indPerc = (int)(double(ic) * 100.0 / double(nChunks - 1));
				if (indPerc != indPercOld) {
					printf("%i %%\n", indPerc);
					indPercOld = indPerc;
				}
			}
		}
	}

	free(vDone); vDone = nullptr;

	return 0;
}

int MergeSort::startProcessing() {
	if (openRootFile() != 0) return -1;
	if (getKeys() != 0) return -2;
	if (allocateMemory() != 0) return -3;
	if (getChunks() != 0) return -4;
	if (prepareMT() != 0) return -5;
	if (runMT() != 0) return -6;
	
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
	MergeSort* se;
	
	if (argc != 4) {
		printUsage();
		return -99;
	}

	auto start = std::chrono::high_resolution_clock::now();

	se = new MergeSort();

	sprintf(se->inputFileName, "%s", argv[1]);
	sprintf(se->outputFolder, "%s", argv[2]);
	sprintf(se->prefix, "%s", argv[3]);
	
	iResult = se->startProcessing();
	delete se; se = nullptr;

	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsed = finish - start;
	printf("Elapsed Time: %10.2f s", 0.001*elapsed.count());

	return iResult;
}
