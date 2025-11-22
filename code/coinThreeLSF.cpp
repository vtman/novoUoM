//D:\NOVO\conData\out test 400

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <chrono>

const long long BUFF_LENGTH = 3*40000;
const long long LEN_D = 22 * 2;

class CoinThree {
public:
	CoinThree();
	~CoinThree();
	
	int processThree();
	int startProcessing();
	int allocateMemory();
	int printInfo();

	char inputFileName[1000], outputFileName[1000], ioFolder[1000], prefix[1000];

	FILE *fo, *fi, *fInfo;
	int iDiffOther;
	char* bufTemp, * bufInput, * bufOutput;
	long long* matCounter, nTriples;
};

CoinThree::CoinThree() {
	bufInput = nullptr;
	bufOutput = nullptr;
	bufTemp = nullptr;
	matCounter = nullptr;

	fo = nullptr;
	fi = nullptr;
	fInfo = nullptr;
}

CoinThree::~CoinThree() {
	if (bufInput != nullptr) { free(bufInput); bufInput = nullptr;}
	if (bufOutput != nullptr) { free(bufOutput); bufOutput = nullptr;}
	if (bufTemp != nullptr) { free(bufTemp); bufTemp = nullptr;}
	if (matCounter != nullptr) { free(matCounter); matCounter = nullptr; }
	if (fo != nullptr) { fclose(fo); fo = nullptr;}
	if (fi != nullptr) { fclose(fi); fi = nullptr;}
	if (fInfo != nullptr) { fclose(fInfo); fInfo = nullptr; }
}

int CoinThree::printInfo() {
	int iB1, iB2, iB3, iC1, iC2, iC3;

	for (int k = 0; k < 128; k++) {
		iB1 = k >> 3;
		iC1 = (k & 7) << 1;
		for (int i = k + 1; i < 128; i++) {
			iB2 = i >> 3;
			iC2 = (i & 7) << 1;
			for (int j = i + 1; j < 128; j++) {
				iB3 = j >> 3;
				iC3 = (j & 7) << 1;
				if (matCounter[128 * (128 * k + i) + j] == 0)continue;
				fprintf(fInfo, "dt%c_%i\tdt%c_%i\tdt%c_%i\t%lli\n", 'a' + iB1, iC1, 'a' + iB2, iC2, 'a' + iB3, iC3, matCounter[128 * (128 * k + i) + j]);
			}
		}
	}

	fprintf(fInfo, "\nNumber of triples: %lli\n", nTriples);
	return 0;
}

int CoinThree::allocateMemory() {
	bufTemp = (char*)malloc(sizeof(char) * LEN_D * BUFF_LENGTH);
	bufInput = (char*)malloc(sizeof(char) * LEN_D * BUFF_LENGTH);
	bufOutput = (char*)malloc(sizeof(char) * LEN_D * BUFF_LENGTH);
	matCounter = (long long*)malloc(sizeof(long long) * 128 * 128 *128);

	for (int i = 0; i < 128*128*128; i++) {
		matCounter[i] = 0;
	}

	if (bufInput == nullptr || bufOutput == nullptr || bufTemp == nullptr || matCounter == nullptr) {
		printf("Error: buffer\n");
		return -1;
	}

	sprintf(inputFileName, "%s/%s_cleanPairs.bin", ioFolder, prefix);
	sprintf(outputFileName, "%s/%s_coinThree.bin", ioFolder, prefix);

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

	sprintf(outputFileName, "%s/%s_coinThree.txt", ioFolder, prefix);
	fInfo = fopen(outputFileName, "w");
	if (fInfo == nullptr) {
		printf("Error: cannot open output file %s\n", outputFileName);
		return -4;
	}

	return 0;
}


int CoinThree::processThree() {
	char* vrow1, * vrow2, *vrow3, *vr[3], *vtemp;
	long long fileSize, startPosition, nDone, nBefore, nRead;
	long long nRows, nEntries, nLeft;
	long long Tlast, T, T1, T2, T3, TT;

	int iCB1, iCB2, iCB3, iCBval[3],  iCBtemp, iPerc, iPercOld, indBuf, iDiff;
	bool isLast, isAgain;

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

	nTriples = 0;

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

		for (long long i = 0; i < nDone; i++) {
			vrow1 = bufInput + i * LEN_D;
			iCB1 = (int)(*(vrow1 + 16)) >> 1;

			T1 = *(long long*)vrow1;

			TT = T1 + iDiff;
			for (long long j = i + 1; j < nRows; j++) {
				vrow2 = bufInput + j * LEN_D;
				iCB2 = (int)(*(vrow2 + 16)) >> 1;
				if (iCB1 == iCB2) continue;
				T2 = *(long long*)vrow2;
				if (T2 > TT) break;

				for (long long k = j + 1; k < nRows; k++) {
					vrow3 = bufInput + k * LEN_D;
					iCB3 = (int)(*(vrow3 + 16)) >> 1;
					if (iCB1 == iCB3) continue;
					if (iCB2 == iCB3) continue;
					T3 = *(long long*)vrow3;
					if (T3 > TT) break;

					vr[0] = vrow1;
					vr[1] = vrow2;
					vr[2] = vrow3;
					iCBval[0] = iCB1;
					iCBval[1] = iCB2;
					iCBval[2] = iCB3;

					do {
						isAgain = false;
						for (int m = 0; m < 2; m++) {
							if (iCBval[m] > iCBval[m + 1]) {
								isAgain = true;
								iCBtemp = iCBval[m];
								iCBval[m] = iCBval[m + 1];
								iCBval[m + 1] = iCBtemp;
								vtemp = vr[m];
								vr[m] = vr[m + 1];
								vr[m + 1] = vtemp;
							}
						}
					} while (isAgain);

					matCounter[128 * (128 * iCBval[0] + iCBval[1]) + iCBval[2]]++;
					memcpy(bufOutput + indBuf * LEN_D, vr[0], LEN_D * sizeof(char));
					memcpy(bufOutput + (indBuf + 1) * LEN_D, vr[1], LEN_D * sizeof(char));
					memcpy(bufOutput + (indBuf + 2) * LEN_D, vr[2], LEN_D * sizeof(char));
					nTriples++;
					indBuf += 3;
					if (indBuf == BUFF_LENGTH) {
						fwrite(bufOutput, sizeof(char), BUFF_LENGTH * LEN_D, fo);
						indBuf = 0;
					}
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

	printf("Number of triples: %lli\n", nTriples);

	return 0;
}

int CoinThree::startProcessing() {
	if (allocateMemory() != 0) return -1;
	if (processThree() != 0) return -2;
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
	CoinThree* cbs;
	
	if (argc != 4) {
		printUsage();
		return -99;
	}

	auto start = std::chrono::high_resolution_clock::now();

	cbs = new CoinThree();

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
