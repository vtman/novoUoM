//D:\NOVO\conData\out test 5

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <chrono>

const long long BUFF_LENGTH = 100000;
const long long LEN_D = 22 * 2;

class CoinPair {
public:
	CoinPair();
	~CoinPair();
	
	int processTwo();
	int startProcessing();
	int allocateMemory();
	int printInfo();

	char inputFileName[1000], outputFileName[1000], ioFolder[1000], prefix[1000];

	FILE *fo, *fi, *fInfo;
	int iDiffOther;
	char* bufTemp, * bufInput, * bufOutput;
	int* vActive, nActive;
	long long* matCounter;
};

CoinPair::CoinPair() {
	bufInput = nullptr;
	bufOutput = nullptr;
	bufTemp = nullptr;
	matCounter = nullptr;
	vActive = nullptr;

	fo = nullptr;
	fi = nullptr;
	fInfo = nullptr;
}

CoinPair::~CoinPair() {
	if (bufInput != nullptr) { free(bufInput); bufInput = nullptr;}
	if (bufOutput != nullptr) { free(bufOutput); bufOutput = nullptr;}
	if (bufTemp != nullptr) { free(bufTemp); bufTemp = nullptr;}
	if (matCounter != nullptr) { free(matCounter); matCounter = nullptr; }
	if (vActive != nullptr) { free(vActive); vActive = nullptr; }
	if (fo != nullptr) { fclose(fo); fo = nullptr;}
	if (fi != nullptr) { fclose(fi); fi = nullptr;}
	if (fInfo != nullptr) { fclose(fInfo); fInfo = nullptr; }
}

int CoinPair::printInfo() {
	int iBoard, iChannel;
	int ind, ind2;
	long long tot;
	tot = 0;

	nActive = 0;
	for (int i = 0; i < 128; i++) {
		tot = 0;
		for (int j = 0; j < 128; j++) {
			tot += matCounter[i * 128 + j];
		}
		if (tot == 0) continue;
		vActive[nActive] = i;
		nActive++;
	}

	printf("Number of active channels: %i\n", nActive);

	fprintf(fInfo, "Channel");
	for (int i = 0; i < nActive; i++) {
		ind = vActive[i];
		iBoard = ind >> 3;
		iChannel = (ind & 7) << 1;
		fprintf(fInfo, "\tdt%c_%i", 'a' + iBoard, iChannel);
	}
	fprintf(fInfo, "\n");

	for (int i = 0; i < nActive; i++) {
		ind = vActive[i];
		iBoard = ind >> 3;
		iChannel = (ind & 7) << 1;
		fprintf(fInfo, "dt%c_%i", 'a' + iBoard, iChannel);

		for (int j = 0; j < nActive; j++) {
			ind2 = vActive[j];
			fprintf(fInfo, "\t%lli", matCounter[ind * 128 + ind2]);
		}
		fprintf(fInfo, "\n");
	}
	return 0;
}

int CoinPair::allocateMemory() {
	bufTemp = (char*)malloc(sizeof(char) * LEN_D * BUFF_LENGTH);
	bufInput = (char*)malloc(sizeof(char) * LEN_D * BUFF_LENGTH);
	bufOutput = (char*)malloc(sizeof(char) * LEN_D * BUFF_LENGTH);
	matCounter = (long long*)malloc(sizeof(long long) * 128 * 128);
	vActive = (int*)malloc(sizeof(int) * 128);

	for (int i = 0; i < 128*128; i++) {
		matCounter[i] = 0;
	}

	if (bufInput == nullptr || bufOutput == nullptr || bufTemp == nullptr) {
		printf("Error: buffer\n");
		return -1;
	}

	sprintf(inputFileName, "%s/%s_cleanPairs.bin", ioFolder, prefix);
	sprintf(outputFileName, "%s/%s_coinTwo.bin", ioFolder, prefix);
	

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

	sprintf(outputFileName, "%s/%s_matrixCoincidence.txt", ioFolder, prefix);
	fInfo = fopen(outputFileName, "w");
	if (fInfo == nullptr) {
		printf("Error: cannot open output file %s\n", outputFileName);
		return -4;
	}

	return 0;
}


int CoinPair::processTwo() {
	char* vrow, * vrow2;
	long long fileSize, startPosition, nDone, nBefore, nRead;
	long long nRows, nEntries, nLeft;
	long long Tlast, T, T1, T2, TT;
	long long nPairs;

	int iCB1, iCB2, iPerc, iPercOld, indBuf, iDiff;
	bool isLast;

	iDiff = 1024 * iDiffOther;

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

		for (int i = 0; i < 30; i++) {
			vrow = bufInput + i * LEN_D;
			T1 = *(long long*)vrow;
			printf("%lli\n", T1);
		}

		for (long long i = 0; i < nDone; i++) {
			vrow = bufInput + i * LEN_D;
			iCB1 = (int)(*(vrow + 16)) >> 1;

			T1 = *(long long*)vrow;

			TT = T1 + iDiff;
			for (long long j = i + 1; j < nRows; j++) {
				vrow2 = bufInput + j * LEN_D;
				iCB2 = (int)(*(vrow2 + 16)) >> 1;
				if (iCB1 == iCB2) continue;
				T2 = *(long long*)vrow2;
				if (T2 > TT) break;

				matCounter[128 * iCB1 + iCB2]++;
				matCounter[128 * iCB2 + iCB1]++;

				memcpy(bufOutput + indBuf * LEN_D, vrow, LEN_D * sizeof(char));
				indBuf++;
				memcpy(bufOutput + indBuf * LEN_D, vrow2, LEN_D * sizeof(char));
				indBuf++;
				nPairs++;

				if (indBuf == BUFF_LENGTH) {
					fwrite(bufOutput, sizeof(char), BUFF_LENGTH * LEN_D, fo);
					indBuf = 0;
				}
			}
		}

		nBefore = nRows - nDone;
		if (nBefore > 0) {
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
	if (processTwo() != 0) return -2;
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
	cbs->iDiffOther = atoi(argv[3]);
	
	iResult = cbs->startProcessing();
	delete cbs; cbs = nullptr;

	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsed = finish - start;
	printf("Elapsed Time: %10.2f s", 0.001*elapsed.count());

	return iResult;
}
