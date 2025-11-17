//D:\NOVO\conData\in\det_000206.root D:\NOVO\conData\out test
//D:\NOVO\2025_02_Cf252\test\temp 211 2 4 10

#include <stdio.h>
#include <fstream>
#include <chrono>
#include <omp.h>

const int64_t LINES_PER_BLOCK = 1000000;
const int LINES_PER_GAP = 100;
const int MAX_BOARDS = 2;
const int LEN_D = 18;

class BlockData {
public:
	BlockData();
	~BlockData();

	int allocateMemory();
	int processBlock();

	int iState;
	int nRecords, nChannels, nPairs, nLen;
	int iChannel, iBoard, ind, indPair;
	int im;
	int64_t iTime, iDiff, iValue;
	int64_t minShort, minLong, maxShort, maxLong;
	int nBitsShort, nBitsLong;

	char* vBuffer;
	int64_t* vLastPosition, *vDistribution;
};

class TimeDiff {
public:
	TimeDiff();
	~TimeDiff();

	int checkParameters();
	int startProcessing();
	int blockProcessing();
	int writeShort();
	int writeLong();
	int allocateMemory();
	int countPairs();
	int getGapData();
	int prepareMT();
	int mergeData();

	char inputFileName[1000], outputFileName[1000], outputFolder[1000], prefix[1000];
	char ioFileName[1000];
	char* vBuf;

	FILE* fi, *foShort, *foLong;

	int nBlocks;
	int nChannels, nPairs;
	int nBitsShort, nBitsLong, nSamples;
	int iMinShort, iMinLong, nLen;
	int nStepShort, nStepLong;
	int maxGapTime, nThreads;

	bool* vBool;

	BlockData** vBD;

	int64_t maxShort, maxLong;

	int64_t* vLastPosition;
	int64_t* vDistribution, *vPairCount;
	int64_t nRecords;
	int* vActiveShort, * vActiveLong, *vState;
	int nActiveShort, nActiveLong;
};

BlockData::BlockData() {
	vLastPosition = nullptr;
	vDistribution = nullptr;
	vBuffer = nullptr;

	iState = 0;
}

BlockData::~BlockData() {
	if (vLastPosition != nullptr) { free(vLastPosition); vLastPosition = nullptr; }
	if (vDistribution != nullptr) { free(vDistribution); vDistribution = nullptr; }
	if (vBuffer != nullptr) { free(vBuffer); vBuffer = nullptr; }
}

int BlockData::allocateMemory() {
	nPairs = nChannels * nChannels;
	vBuffer = (char*)malloc(sizeof(char) * LEN_D * LINES_PER_BLOCK);
	if (vBuffer == nullptr) return -1;
	vLastPosition = (int64_t*)malloc(sizeof(int64_t) * nChannels);
	vDistribution = (int64_t*)malloc(sizeof(int64_t) * nPairs * nLen);

	for (int i = 0; i < nPairs * nLen; i++) {
		vDistribution[i] = 0;
	}
	return 0;
}

TimeDiff::TimeDiff() {
	nThreads = 0;
	fi = nullptr;
	foShort = nullptr;
	foLong = nullptr;

	vLastPosition = nullptr;
	vDistribution = nullptr;
	vActiveLong = nullptr;
	vActiveShort = nullptr;
	vBuf = nullptr;
	vState = nullptr;
	vBool = nullptr;
}

TimeDiff::~TimeDiff() {
	if (vBD != nullptr) {
		for (int i = 0; i < nThreads; i++) {
			if (vBD[i] != nullptr) {
				delete vBD[i]; vBD[i] = nullptr;
			}
		}
		free(vBD); vBD = nullptr;
	}
	if (fi != nullptr) { fclose(fi); fi = nullptr; }
	if (foShort != nullptr) { fclose(foShort); foShort = nullptr; }
	if (foLong != nullptr) { fclose(foLong); foLong = nullptr; }

	if (vBuf != nullptr) { free(vBuf); vBuf = nullptr; }
	if (vLastPosition != nullptr) { free(vLastPosition); vLastPosition = nullptr; }
	if (vDistribution != nullptr) { free(vDistribution); vDistribution = nullptr; }
	if (vActiveLong != nullptr) { free(vActiveLong); vActiveLong = nullptr; }
	if (vActiveShort != nullptr) { free(vActiveShort); vActiveShort = nullptr; }
	if (vState != nullptr) { free(vState); vState = nullptr; }
	if (vBool != nullptr) { free(vBool); vBool = nullptr; }
}




int TimeDiff::checkParameters() {
	long long fileSize;
	sprintf(ioFileName, "%s/%s_merged.bin", outputFolder, prefix);
	fi = fopen(ioFileName, "rb");
	if (fi == nullptr) {
		printf("Error: cannot open file %s\n", ioFileName);
		return -1;
	}

#ifdef _WIN32
	_fseeki64(fi, 0, SEEK_END);
	fileSize = _ftelli64(fi);
	_fseeki64(fi, 0, SEEK_SET);
#else
	fseeko64(fi, 0, SEEK_END);
	fileSize = ftello64(fi);
	fseeko64(fi, 0, SEEK_SET);
#endif

	printf("File size: %lli\n", fileSize);
	nRecords = fileSize / LEN_D;
	printf("Number of records: %lli\n", nRecords);

	nBlocks = nRecords / LINES_PER_BLOCK;
	if (nRecords % LINES_PER_BLOCK > 0) nBlocks++;

	printf("Number of blocks: %i\n", nBlocks);

	printf("Bits (short): %i\n", nBitsShort);
	if (nBitsShort < 0) {
		printf("Error: wrong value");
		return -2;
	}

	nStepShort = 1 << nBitsShort;
	printf("Step size (short): %i\n", nStepShort);

	printf("Bits (long): %i\n", nBitsLong);
	if (nBitsLong < 0) {
		printf("Error: wrong value");
		return -3;
	}

	nStepLong = 1 << nBitsLong;
	printf("Step size (long): %i\n", nStepLong);

	nLen = 2 * nSamples + 1;

	iMinShort = -nSamples * nStepShort - nStepShort/2;
	printf("Range (short): [%i, %i]\n", iMinShort, iMinShort + nLen * nStepShort - 1);

	iMinLong = -nSamples * nStepLong - nStepLong/2;
	printf("Range (long): [%i, %i]\n", iMinLong, iMinLong + nLen * nStepLong - 1);


	maxShort = nLen * nStepShort;
	maxLong = nLen * nStepLong;

	if (maxShort < maxLong) {
		maxGapTime = maxLong;
	}
	else {
		maxGapTime = maxShort;
	}

	return 0;
}

int TimeDiff::allocateMemory() {
	nChannels = 16 * MAX_BOARDS;
	nPairs = nChannels * nChannels;

	vActiveLong = (int*)malloc(sizeof(int) * nPairs);
	vActiveShort = (int*)malloc(sizeof(int) * nChannels);

	vPairCount = (int64_t*)malloc(sizeof(int64_t) * nPairs);
	vDistribution = (int64_t*)malloc(sizeof(int64_t) * nPairs * nLen);
	vLastPosition = (int64_t*)malloc(sizeof(int64_t) * nChannels * nBlocks);

	vBuf = (char*)malloc(sizeof(char) * LEN_D * LINES_PER_GAP);

	vState = (int*)malloc(sizeof(int) * nBlocks);
	vBool = (bool*)malloc(sizeof(bool) * nChannels);

	return 0;
}

int TimeDiff::countPairs() {
	int64_t totCount;
	int64_t* vLoc;
	int indP;
	for (int i = 0; i < nChannels; i++) {
		for (int j = 0; j < nChannels; j++) {
			totCount = 0;
			indP = i * nChannels + j;
			vLoc = vDistribution + indP * nLen;
			for (int k = 0; k < nLen; k++) {
				totCount += vLoc[k];
			}
			vPairCount[indP] = totCount;
		}
	}
	return 0;
}

int BlockData::processBlock() {
	char* vLoc;

	for (int i = 0; i < nRecords; i++) {
		vLoc = vBuffer + i * LEN_D;
		iChannel = vLoc[16];
		iBoard = vLoc[17];
		if (iBoard >= MAX_BOARDS) {
			iState |= 1;
			continue;
		}
		if (iChannel >= 16) {
			iState |= 2;
			continue;
		}
		iTime = *(int64_t*)vLoc;

		ind = iBoard * 16 + iChannel;

		if ((iChannel & 15) < 15) {
			indPair = ind + 1 - 2 * (iChannel & 1);

			iDiff = iTime - vLastPosition[indPair];

			iValue = iDiff - minShort;
			if (iValue >= 0 && iValue < maxShort) {
				im = ((int32_t)iValue) >> nBitsShort;
				vDistribution[(ind * nChannels + indPair) * nLen + im]++;
			}
			iValue = -iDiff - minShort;
			if (iValue >= 0 && iValue < maxShort) {
				im = ((int32_t)iValue) >> nBitsShort;
				vDistribution[(indPair * nChannels + ind) * nLen + im]++;
			}
		}
		for (int j = 0; j < nChannels; j++) {
			if (ind == j)continue;
			if (ind == indPair)continue;
			iDiff = iTime - vLastPosition[j];
			iValue = iDiff - minLong;
			if (iValue >= 0 && iValue < maxLong) {
				im = ((int32_t)iValue) >> nBitsLong;
				vDistribution[(ind * nChannels + j) * nLen + im]++;
			}
			iValue = -iDiff - minLong;
			if (iValue >= 0 && iValue < maxLong) {
				im = ((int32_t)iValue) >> nBitsLong;
				vDistribution[(j * nChannels + ind) * nLen + im]++;
			}

		}
		vLastPosition[ind] = iTime;
	}

	return 0;
}

int TimeDiff::blockProcessing() {
	int indPerc, indPercOld;

	indPercOld = -1;


	for (int i = 0; i < nBlocks; i++) {
		vState[i] = 0;
	}

#pragma omp parallel
	{
		int tid;
		bool isContinue;
		BlockData* bd;
		int64_t nLeft, startRecord, startPosition;

		tid = omp_get_thread_num();
		bd = vBD[tid];


		for (int ib = 0; ib < nBlocks; ib++) {
#pragma omp critical
			{
				isContinue = true;
				if (vState[ib] == 1) {
					isContinue = false;
				}
				else {
					vState[ib] = 1;
					indPerc = (int)(double(ib) * 100.0 / double(nBlocks - 1));
					if (indPerc != indPercOld) {
						printf("%i %%\n", indPerc);
						indPercOld = indPerc;
					}
				}
			}
			if (!isContinue) continue;

#pragma omp critical
			{
				
				startRecord = ((int64_t)ib) * LINES_PER_BLOCK;
				nLeft = nRecords - startRecord;
				if (nLeft > LINES_PER_BLOCK) nLeft = LINES_PER_BLOCK;
				startPosition = startRecord * LEN_D;
#ifdef _WIN32
				_fseeki64(fi, startPosition, SEEK_SET);
#else
				fseeko64(fi, startPosition, SEEK_SET);
#endif
				fread(bd->vBuffer, sizeof(char), LEN_D * nLeft, fi);
			}
			

			for (int i = 0; i < nChannels; i++) {
				bd->vLastPosition[i] = vLastPosition[ib * nChannels + i];
			}
			bd->nRecords = nLeft;
			bd->processBlock();
		}
	}

	return 0;
}

int TimeDiff::writeShort() {
	int i1, i2;
	nActiveShort = 0;
	for (int i = 0; i < nChannels; i += 2) {
		if (vPairCount[i * nChannels + i + 1] == 0)continue;
		vActiveShort[nActiveShort] = i * nChannels + i + 1;
		nActiveShort++;
		vActiveShort[nActiveShort] = (i + 1) * nChannels + i;
		nActiveShort++;
	}

	sprintf(ioFileName, "%s/%s_timeDifferenceShort.txt", outputFolder, prefix);
	foShort = fopen(ioFileName, "w");
	if (foShort == nullptr) {
		printf("Error: cannot open file %s\n", ioFileName);
		return -1;
	}

	fprintf(foShort, "timeDiff");

	for (int i = 0; i < nActiveShort; i++) {
		i1 = vActiveShort[i] / nChannels;
		i2 = vActiveShort[i] % nChannels;
		fprintf(foShort, "\t%03i-%03i", i1, i2);
	}
	fprintf(foShort, "\n");
	for (int k = 0; k < nLen; k++) {
		fprintf(foShort, "%i", iMinShort + k * nStepShort + nStepShort/2);
		for (int i = 0; i < nActiveShort; i++) {
			fprintf(foShort, "\t%lli", vDistribution[vActiveShort[i] * nLen + k]);
		}
		fprintf(foShort, "\n");
	}

	fclose(foShort); foShort = nullptr;

	/*
	for (int i = 0; i < nActiveShort; i++) {
		i1 = vActiveShort[i] / nChannels;
		i2 = vActiveShort[i] % nChannels;
		sprintf(ioFileName, "%s/data/%s_timeDifference_%i-%i.txt", outputFolder, prefix, i1, i2);
		foShort = fopen(ioFileName, "w");
		if (foShort == nullptr) {
			printf("Error: cannot open file %s\n", ioFileName);
			return -1;
		}

		for (int k = 0; k < nLen; k++) {
			fprintf(foShort, "%i", iMinShort + k * nStepShort + nStepShort / 2);
			fprintf(foShort, "\t%lli", vDistribution[vActiveShort[i] * nLen + k]);
			fprintf(foShort, "\n");
		}
		fclose(foShort); foShort = nullptr;
	}
	*/
	return 0;
}

int TimeDiff::writeLong() {
	int i1, i2;
	nActiveLong = 0;
	for (int i = 0; i < nChannels; i++) {
		if ((i & 15) > 13)continue;
		for (int j = 0; j < nChannels; j++) {
			if (i == j)continue;
			if ((j & 15) > 13)continue;
			if ((i >> 1) == (j >> 1)) continue;
			if (vPairCount[i * nChannels + j] == 0)continue;
			vActiveLong[nActiveLong] = i * nChannels + j;
			nActiveLong++;
		}
	}

	sprintf(ioFileName, "%s/%s_timeDifferenceLong.txt", outputFolder, prefix);
	foLong = fopen(ioFileName, "w");
	if (foLong == nullptr) {
		printf("Error: cannot open file %s\n", ioFileName);
		return -1;
	}

	fprintf(foLong, "timeDiff");

	for (int i = 0; i < nActiveLong; i++) {
		i1 = vActiveLong[i] / nChannels;
		i2 = vActiveLong[i] % nChannels;
		fprintf(foLong, "\t%03i-%03i", i1, i2);
	}
	fprintf(foLong, "\n");
	for (int k = 0; k < nLen; k++) {
		fprintf(foLong, "%i", iMinLong + k * nStepLong + nStepLong/2);
		for (int i = 0; i < nActiveLong; i++) {
			fprintf(foLong, "\t%lli", vDistribution[vActiveLong[i] * nLen + k]);
		}
		fprintf(foLong, "\n");
	}

	fclose(foLong); foLong = nullptr;

	/*
	for (int i = 0; i < nActiveLong; i++) {
		i1 = vActiveLong[i] / nChannels;
		i2 = vActiveLong[i] % nChannels;
		sprintf(ioFileName, "%s/data/%s_timeDifference_%i-%i.txt", outputFolder, prefix, i1, i2);
		foLong = fopen(ioFileName, "w");
		if (foLong == nullptr) {
			printf("Error: cannot open file %s\n", ioFileName);
			return -1;
		}

		for (int k = 0; k < nLen; k++) {
			fprintf(foLong, "%i", iMinLong + k * nStepLong + nStepLong / 2);
			fprintf(foLong, "\t%lli", vDistribution[vActiveLong[i] * nLen + k]);
			fprintf(foLong, "\n");
		}
		fclose(foLong); foLong = nullptr;
	}
	*/

	return 0;
}

int TimeDiff::prepareMT() {
#pragma omp parallel
	{
		nThreads = omp_get_num_threads();
	}

	printf("\nNumber of threads: %i\n", nThreads);

	vBD = (BlockData**)malloc(sizeof(BlockData*) * nThreads);
	for (int i = 0; i < nThreads; i++) {
		vBD[i] = new BlockData();
		vBD[i]->nLen = nLen;
		vBD[i]->iState = 0;
		vBD[i]->nBitsShort = nBitsShort;
		vBD[i]->nBitsLong = nBitsLong;
		vBD[i]->minLong = iMinLong;
		vBD[i]->minShort = iMinShort;
		vBD[i]->nChannels = nChannels;
		vBD[i]->maxShort = maxShort;
		vBD[i]->maxLong = maxLong;

		if (vBD[i]->allocateMemory() != 0) {
			printf("Error: block data\n");
			return -1;
		}
	}
	return 0;
}

int TimeDiff::mergeData() {
	for (int i = 0; i < nPairs * nLen; i++) {
		vDistribution[i] = 0;
	}

	for (int k = 0; k < nThreads; k++) {
		for (int i = 0; i < nPairs * nLen; i++) {
			vDistribution[i] += vBD[k]->vDistribution[i];
		}
	}
	return 0;
}

int TimeDiff::getGapData() {
	int64_t iRecord, startPosition, iTime, iTimeTop, iTimeBottom;
	int nGap, nGapMax;
	int iChannel, iBoard, ind;
	char *vLoc;
	for (int i = 0; i < nChannels * nBlocks; i++) {
		vLastPosition[i] = 0;
	}

	nGapMax = -1;
	
	for (int ib = 1; ib < nBlocks; ib++) {
		for (int i = 0; i < nChannels; i++) {
			vBool[i] = true;
		}


		iRecord = ((int64_t)ib) * LINES_PER_BLOCK - LINES_PER_GAP;
		startPosition = iRecord * LEN_D;
#ifdef _WIN32
		_fseeki64(fi, startPosition, SEEK_SET);
#else
		fseeko64(fi, startPosition, SEEK_SET);
#endif
		fread(vBuf, sizeof(char), LEN_D * LINES_PER_GAP, fi);

		iTimeTop = *(int64_t*)(vBuf + (LINES_PER_GAP - 1) * LEN_D);
		iTimeBottom = iTimeTop - maxGapTime;

		nGap = LINES_PER_GAP;

		for (int i = LINES_PER_GAP - 1; i >= 0; i--) {
			vLoc = vBuf + i * LEN_D;
			iChannel = vLoc[16];
			iBoard = vLoc[17];

			if (iBoard >= MAX_BOARDS) continue;
			if (iChannel >= 16) continue;
			iTime = *(int64_t*)vLoc;
			if (iTime < iTimeBottom) {
				nGap = LINES_PER_GAP - i;
				break;
			}

			ind = iBoard * 16 + iChannel;
			if (vBool[ind]) {
				vLastPosition[ind] = iTime;
				vBool[ind] = false;
			}
		}

		if (nGap > nGapMax) nGapMax = nGap;
	}

	printf("Max gap size: %i\n", nGapMax);
	return 0;
}

int TimeDiff::startProcessing() {
	if (checkParameters() != 0) return -1;
	if (allocateMemory() != 0) return -2;
	if (getGapData() != 0) return -3;
	if (prepareMT() != 0) return -4;
	if (blockProcessing() != 0) return -5;
	if (mergeData() != 0) return -6;
	if (countPairs() != 0) return -7;
	if (writeShort() != 0) return -8;
	if (writeLong() != 0) return -9;
	
	return 0;
}

int printUsage() {
	printf("Input parameters:\n");
	printf("\t1) IO folder (string)\n");
	printf("\t2) prefix (string)\n");
	printf("\t3) bits short (int, >= 0)\n");
	printf("\t4) bits long (int, >= 0)\n");
	printf("\t5) samples (int, > 0)\n\n");
	return 0;
}

int main(int argc, char** argv) {
	int iResult;
	TimeDiff* td;
	
	if (argc != 6) {
		printUsage();
		return -99;
	}

	auto start = std::chrono::high_resolution_clock::now();

	td = new TimeDiff();

	sprintf(td->outputFolder, "%s", argv[1]);
	sprintf(td->prefix, "%s", argv[2]);
	td->nBitsShort = atoi(argv[3]);
	td->nBitsLong = atoi(argv[4]);
	td->nSamples = atoi(argv[5]);
	
	iResult = td->startProcessing();
	delete td; td = nullptr;

	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsed = finish - start;
	printf("Elapsed Time: %10.2f s", 0.001*elapsed.count());

	return iResult;
}
