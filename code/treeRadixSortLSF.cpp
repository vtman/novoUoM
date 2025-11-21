//D:\NOVO\conData\out test

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <chrono>

#include <omp.h>

const int64_t LONG_LENGTH = 1000000;
const int64_t SHORT_LENGTH = 10000;
const int64_t LEN_D = 22;

class SortTreeData {
public:
	SortTreeData();
	~SortTreeData();

	int getMaxBoards();
	int radixSort(char* vSortBufferMain, char* vSortBuffer0, char* vSortBuffer1, int n);

	int processNew();
	int startProcessing();
	int copyBinaryFile(const char* source, const char* destination, long long *M);

	char inputFileName[1000], ioFileName[1000], prefix[1000], outputFolder[1000];
	char inFileName[1000], outFileName[1000];

	char* treeNames, * branchNames;
	char* vTypeName;

	char** pvBuffer;

	int nThreads, maxBoards;

	int* vStatus, * vSame;
	int* vIndDTkeys, * ivTemp;
	int* vTypeSize, * vArraySize;
	int nKeys, nDTkeys, nBranches;

	long long* vNumEntries;
};


int SortTreeData::copyBinaryFile(const char* source, const char* destination, long long* M) {
	char* buf;
	long long nBuf, m;
	size_t bytes;
	FILE* fi, * fo;

	fi = fopen(source, "rb");
	if (fi == nullptr) {
		printf("Warning: cannot open file %s\n", source);
		return -1;
	}

	fo = fopen(destination, "wb");
	if (fo == nullptr) {
		printf("Error: cannot open file %s\n", destination);
		return -2;
	}

	m = 0;

	nBuf = LEN_D * LONG_LENGTH;
	buf = pvBuffer[0];

	while ((bytes = fread(buf, 1, nBuf, fi)) > 0) {
		fwrite(buf, 1, bytes, fo);
		m += bytes;
	}

	fclose(fi); fi = nullptr;
	fclose(fo); fo = nullptr;

	*M = m / LEN_D;

	return 0;
}

SortTreeData::SortTreeData() {
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
	if (vIndDTkeys != nullptr) { free(vIndDTkeys); vIndDTkeys = nullptr; }
	if (ivTemp != nullptr) { free(ivTemp); ivTemp = nullptr; }
	if (vNumEntries != nullptr) { free(vNumEntries); vNumEntries = nullptr; }
	if (treeNames != nullptr) { free(treeNames); treeNames = nullptr; }
	if (vStatus != nullptr) { free(vStatus); vStatus = nullptr; }
	if (vSame != nullptr) { free(vSame); vSame = nullptr; }
	if (branchNames != nullptr) { free(branchNames); branchNames = nullptr; }
	if (vTypeSize != nullptr) { free(vTypeSize); vTypeSize = nullptr; }
	if (vArraySize != nullptr) { free(vArraySize); vArraySize = nullptr; }
	if (vTypeName != nullptr) { free(vTypeName); vTypeName = nullptr; }
	if (pvBuffer != nullptr) {
		for (int i = 0; i < 3 * nThreads; i++) {
			if (pvBuffer[i] != nullptr) { free(pvBuffer[i]); pvBuffer[i] = nullptr; }
		}
		free(pvBuffer); pvBuffer = nullptr;
	}

}

int SortTreeData::radixSort(char* vSortBufferMain, char* vSortBuffer0, char* vSortBuffer1, int n) {
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
			if ((((int)(*(vSortBufferMain + k * LEN_D + 16)) >> i) & 1) == 0) {
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
			if ((((*(long long*)(vSortBufferMain + k * LEN_D)) >> i) & 1) == 0) {
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

int SortTreeData::getMaxBoards() {
	FILE* fi;
	sprintf(inFileName, "%s/%s_maxBoards.txt", outputFolder, prefix);
	fi = fopen(inFileName, "r");
	if (fi == nullptr) {
		printf("Error: cannot open file %s\n", inFileName);
		return -1;
	}
	fscanf(fi, "%i", &maxBoards);
	fclose(fi); fi = nullptr;
	return 0;
}

int SortTreeData::processNew() {
	FILE* fio;

	long long maxRowsInChunk, M;
	long long* vStart, * vEnd;
	int nBlocks, iPercOld, iPerc, nSmallRows;
	int* vDone;

	maxRowsInChunk = LONG_LENGTH / LEN_D;
	nSmallRows = SHORT_LENGTH / LEN_D;

	pvBuffer = (char**)malloc(sizeof(char*) * 3 * nThreads);

	for (int i = 0; i < 3 * nThreads; i++) {
		pvBuffer[i] = nullptr;
		pvBuffer[i] = (char*)malloc(sizeof(char) * LEN_D * LONG_LENGTH);
		if (pvBuffer[i] == nullptr) {
			printf("Error: cannot allocate memory\n");
			return -1;
		}
	}

	vNumEntries = (long long*)malloc(sizeof(long long) * maxBoards);

	for (int it = 0; it < maxBoards; it++) {
		sprintf(inFileName, "%s/%s_tree_%i.bin", outputFolder, prefix, it);
		sprintf(outFileName, "%s/%s_tree_sorted_%i.bin", outputFolder, prefix, it);
		if (copyBinaryFile(inFileName, outFileName, vNumEntries + it) != 0) continue;

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

int SortTreeData::startProcessing() {
	if (getMaxBoards()!= 0) return -1;

#pragma omp parallel
	{
		nThreads = omp_get_num_threads();
	}

	printf("\nNumber of threads: %i\n", nThreads);

	if (processNew() != 0) return -2;

	return 0;
}

int printUsage() {
	printf("Input parameters:\n");
	printf("\t1) IO folder (string)\n");
	printf("\t2) prefix (string)\n\n");
	return 0;
}

int main(int argc, char** argv) {
	int iResult;
	SortTreeData* sd;

	if (argc != 3) {
		printUsage();
		return -99;
	}

	auto start = std::chrono::high_resolution_clock::now();

	sd = new SortTreeData();

	sprintf(sd->outputFolder, "%s", argv[1]);
	sprintf(sd->prefix, "%s", argv[2]);

	iResult = sd->startProcessing();
	delete sd; sd = nullptr;

	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsed = finish - start;
	printf("Elapsed Time: %10.2f s", 0.001 * elapsed.count());

	return iResult;
}
