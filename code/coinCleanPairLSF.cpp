//D:\NOVO\conData\out test 5

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <chrono>

const long long BUFF_LENGTH = 100000;
const long long LEN_D = 22;

class CoinPair {
public:
	CoinPair();
	~CoinPair();
	
	int processNew();
	int startProcessing();
	int allocateMemory();
	int printInfo();

	char inputFileName[1000], outputFileName[1000], ioFolder[1000], prefix[1000];

	FILE *fo, *fi, *fInfo;
	int iDiffSame;
	char* bufTemp, * bufInput, * bufOutput;
	int* vStatus;
	long long* vCounter;
};

CoinPair::CoinPair() {
	bufInput = nullptr;
	bufOutput = nullptr;
	bufTemp = nullptr;
	vStatus = nullptr;
	vCounter = nullptr;

	fo = nullptr;
	fi = nullptr;
	fInfo = nullptr;
}

CoinPair::~CoinPair() {
	if (bufInput != nullptr) { free(bufInput); bufInput = nullptr;}
	if (bufOutput != nullptr) { free(bufOutput); bufOutput = nullptr;}
	if (bufTemp != nullptr) { free(bufTemp); bufTemp = nullptr;}
	if (vStatus != nullptr) { free(vStatus); vStatus = nullptr;}
	if (vCounter != nullptr) { free(vCounter); vCounter = nullptr; }
	if (fo != nullptr) { fclose(fo); fo = nullptr;}
	if (fi != nullptr) { fclose(fi); fi = nullptr;}
	if (fInfo != nullptr) { fclose(fInfo); fInfo = nullptr; }
}

int CoinPair::printInfo() {
	int iBoard, iChannel;
	long long tot;
	tot = 0;
	for (int i = 0; i < 256; i+=2) {
		if (vCounter[i] == 0)continue;
		iBoard = i / 16;
		iChannel = i;
		fprintf(fInfo, "dt%c\t%i\t%lli\n", 'a' + iBoard, iChannel, vCounter[i]);
		tot += vCounter[i];
	}
	fprintf(fInfo, "\nTotal_pairs: %lli\n", tot);
	return 0;
}

int CoinPair::allocateMemory() {
	bufTemp = (char*)malloc(sizeof(char) * LEN_D * BUFF_LENGTH);
	bufInput = (char*)malloc(sizeof(char) * LEN_D * BUFF_LENGTH);
	bufOutput = (char*)malloc(sizeof(char) * LEN_D * BUFF_LENGTH);
	vStatus = (int*)malloc(sizeof(int) * BUFF_LENGTH);
	vCounter = (long long*)malloc(sizeof(long long) * 256);

	for (int i = 0; i < 256; i++) {
		vCounter[i] = 0;
	}

	if (bufInput == nullptr || bufOutput == nullptr || bufTemp == nullptr || vStatus == nullptr) {
		printf("Error: buffer\n");
		return -1;
	}

	sprintf(inputFileName, "%s/%s_merged.bin", ioFolder, prefix);
	sprintf(outputFileName, "%s/%s_cleanPairs.bin", ioFolder, prefix);

	fi = fopen(inputFileName, "rb");
	if (fi == nullptr) {
		printf("Error: cannot open input file %s\n", inputFileName);
		return -2;
	}

	fo = fopen(outputFileName, "wb");
	if (fo == nullptr) {
		printf("Error: cannot open output file %s\n", outputFileName);
		return -3;
	}

	sprintf(outputFileName, "%s/%s_pairsInfo.txt", ioFolder, prefix);
	fInfo = fopen(outputFileName, "w");
	if (fInfo == nullptr) {
		printf("Error: cannot open output file %s\n", outputFileName);
		return -4;
	}

	return 0;
}


int CoinPair::processNew() {
	char* vrow, * vrow2;
	long long fileSize, startPosition, nDone, nBefore, nRead;
	long long nRows, nEntries, nLeft;
	long long Tlast, T, T1, T2, TT;
	long long nPairs;

	int iCB1, iCB2, iCBodd, iCBeven;
	int iPerc, iPercOld, iChannel2;
	int ind1, ind2, indBuf;
	int iDiff;
	bool isFound, isLast;

	iDiff = 1024 * iDiffSame;

	#ifdef _WIN32
	_fseeki64(fi, 0, SEEK_END);
	fileSize = _ftelli64(fi);
	nEntries = fileSize / LEN_D;
	_fseeki64(fi, 0, SEEK_SET);
	#else
	fseeko64(fi, 0, SEEK_END);
	fileSize = ftello64(fi);
	nEntries = fileSize / LEN_D;
	fseeko64(fi, 0, SEEK_SET);
	#endif

	printf("Number of entries (input): %lli\n", nEntries);

	startPosition = 0;
	nBefore = 0;
	nRows = 0;
	indBuf = 0;

	memset(vStatus, 0, sizeof(int) * BUFF_LENGTH);

	iPercOld = -1;

	nPairs = 0;

	do {
		nLeft = nEntries - startPosition - nBefore;

		iPerc = (int)(100.0 * double(startPosition) / double(nEntries));
		if (iPerc != iPercOld) {
			iPercOld = iPerc;
			printf("%i %%\n", iPerc);
		}
		if (nLeft <= 0) break;
		
		nRows = nBefore;
		nRead = nLeft;
		isLast = true;
		if (nRead > BUFF_LENGTH - nRows) {
			nRead = BUFF_LENGTH - nRows;
			isLast = false;
		}
		fread(bufInput + nRows * LEN_D, sizeof(char), LEN_D * nRead, fi);
		memset(vStatus + nRows, 0, sizeof(int) * nRead);
		nRows += nRead;
		nDone = nRows - 1;
		if (!isLast) {
			Tlast = *(long long*)(bufInput + nDone * LEN_D) - iDiff;
			do {
				nDone--;
				T = *(long long*)(bufInput + nDone * LEN_D);
			} while (T > Tlast);
		}
		nDone++;

		for (long long i = 0; i < nDone; i++) {
			vrow = bufInput + i * LEN_D;
			iCB1 = (int)(*(vrow + 16));

			if (((iCB1 & 15) >> 1) == 7) continue;

			iCBeven = (iCB1 >> 1) << 1;
			iCBodd = iCBeven + 1;

			if (iCB1 == iCBeven) {
				iCB2 = iCBodd;
			}
			else {
				iCB2 = iCBeven;
			}

			T1 = *(long long*)vrow;

			ind1 = i;
			TT = T1 + iDiff;
			isFound = false;
			for (long long j = i + 1; j < nRows; j++) {
				vrow2 = bufInput + j * LEN_D;
				iChannel2 = (int)(*(vrow2 + 16));
				if (iChannel2 != iCB2) continue;
				T2 = *(long long*)vrow2;
				if (T2 > TT) break;
				ind2 = j;
				isFound = true;
				break;
			}
			if (!isFound) continue;

			vStatus[ind1] = 1;
			vStatus[ind2] = 1;

			vCounter[iCBeven]++;

			vrow = bufInput + iCBeven * LEN_D;
			memcpy(bufOutput + indBuf * LEN_D, vrow, LEN_D * sizeof(char));
			indBuf++;
			vrow = bufInput + iCBodd * LEN_D;
			memcpy(bufOutput + indBuf * LEN_D, vrow, LEN_D * sizeof(char));
			indBuf++;
			nPairs++;

			if (indBuf == BUFF_LENGTH) {
				fwrite(bufOutput, sizeof(char), BUFF_LENGTH * LEN_D, fo);
				indBuf = 0;
			}
		}

		nBefore = nRows - nDone;
		if (nBefore > 0) {
			for (int i = 0; i < nBefore; i++) {
				vStatus[i] = vStatus[i + nDone];
			}
			memcpy(bufTemp, bufInput + LEN_D * nDone, sizeof(char) * LEN_D * nBefore);
			memcpy(bufInput, bufTemp, sizeof(char) * LEN_D * nBefore);
		}

		startPosition += nDone;

	} while (true);

	if (indBuf > 0) {
		fwrite(bufOutput, sizeof(char), indBuf * LEN_D, fo);
	}

	printf("Number of pairs: %lli\n", nPairs);

	return 0;
}

int CoinPair::startProcessing() {
	if (allocateMemory() != 0) return -1;
	if (processNew() != 0) return -2;
	if (printInfo() != 0) return -3;
	
	return 0;
}

int printUsage() {
	printf("Input parameters:\n");
	printf("\t1) IO folder (string)\n");
	printf("\t2) prefix (string)\n");
	printf("\t3) max difference for time stamps (in samples, integer, positive)\n");
	return 0;
}

int main(int argc, char** argv) {
	int iResult;
	CoinPair* cbs;
	
	if (argc != 4) {
		printUsage();
		return -99;
	}

	auto start = std::chrono::high_resolution_clock::now();

	cbs = new CoinPair();

	sprintf(cbs->ioFolder, "%s", argv[1]);
	sprintf(cbs->prefix, "%s", argv[2]);
	cbs->iDiffSame = atoi(argv[3]);
	
	iResult = cbs->startProcessing();
	delete cbs; cbs = nullptr;

	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsed = finish - start;
	printf("Elapsed Time: %10.2f s", 0.001*elapsed.count());

	return iResult;
}
