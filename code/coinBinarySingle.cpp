//D:\NOVO\conData\out\test_merged.bin D:\NOVO\conData\out\test_cbs.bin 2 200

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <chrono>

const int BUFF_LENGTH = 10000;

class CoinBinarySingle {
public:
	CoinBinarySingle();
	~CoinBinarySingle();
	
	int processNew();
	int startProcessing();
	int allocateMemory();

	char inputFileName[1000], outputFileName[1000];

	FILE *fo, *fi;
	int iDiffSame, iDiffOther, lenD;
	char* bufTemp, * bufInput, * bufOutput;
	int* vStatus;
};

CoinBinarySingle::CoinBinarySingle() {
	lenD = 18;

	bufInput = nullptr;
	bufOutput = nullptr;
	bufTemp = nullptr;
	vStatus = nullptr;

	fo = nullptr;
	fi = nullptr;
}

CoinBinarySingle::~CoinBinarySingle() {
	if (bufInput != nullptr) {
		free(bufInput); bufInput = nullptr;
	}
	if (bufOutput != nullptr) {
		free(bufOutput); bufOutput = nullptr;
	}
	if (bufTemp != nullptr) {
		free(bufTemp); bufTemp = nullptr;
	}
	if (vStatus != nullptr) {
		free(vStatus); vStatus = nullptr;
	}
	if (fo != nullptr) {
		fclose(fo); fo = nullptr;
	}
	if (fi != nullptr) {
		fclose(fi); fi = nullptr;
	}
}

int CoinBinarySingle::allocateMemory() {
	bufTemp = (char*)malloc(sizeof(char) * lenD * BUFF_LENGTH);
	bufInput = (char*)malloc(sizeof(char) * lenD * BUFF_LENGTH);
	bufOutput = (char*)malloc(sizeof(char) * lenD * BUFF_LENGTH);
	vStatus = (int*)malloc(sizeof(int) * BUFF_LENGTH);

	if (bufInput == nullptr || bufOutput == nullptr || bufTemp == nullptr || vStatus == nullptr) {
		printf("Error: buffer\n");
		return -1;
	}

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

	return 0;
}


int CoinBinarySingle::processNew() {
	char* vrow, * vrow2, * vrow3, * vrow4;
	long long fileSize, startPosition, nDone, nBefore, nRead;
	long long nRows, nEntries, nLeft;
	long long Tlast, T, T1, T2, T3, T4, TT;

	int iPerc, iPercOld, iChannel1, iChannel2, iBoard1, iBoard2;
	int ind1, ind2, ind3, ind4, indBuf;
	bool isFound, isLast;

        #ifdef _WIN32
	_fseeki64(fi, 0, SEEK_END);
	fileSize = _ftelli64(fi);
	nEntries = fileSize / lenD;
	_fseeki64(fi, 0, SEEK_SET);
	#else
	fseeko64(fi, 0, SEEK_END);
	fileSize = ftello64(fi);
	nEntries = fileSize / lenD;
	fseeko64(fi, 0, SEEK_SET);
	#endif

	printf("Number of entries (input): %lli\n", nEntries);

	startPosition = 0;
	nBefore = 0;
	nRows = 0;
	indBuf = 0;

	memset(vStatus, 0, sizeof(int) * BUFF_LENGTH);

	iPercOld = -1;

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
		fread(bufInput + nRows * lenD, sizeof(char), lenD * nRead, fi);
		memset(vStatus + nRows, 0, sizeof(int) * nRead);
		nRows += nRead;
		nDone = nRows - 1;
		if (!isLast) {
			Tlast = *(long long*)(bufInput + nDone * lenD) - iDiffOther;
			do {
				nDone--;
				T = *(long long*)(bufInput + nDone * lenD);
			} while (T > Tlast);
		}
		nDone++;

		for (long long i = 0; i < nDone; i++) {
			vrow = bufInput + i * lenD;
			iChannel1 = (int)(*(vrow + 16));
			if ((iChannel1 & 1) == 1) continue;
			if (iChannel1 == 14) continue;
			T1 = *(long long*)vrow;
			iBoard1 = (int)(*(vrow + 17));

			ind1 = i;
			TT = T1 + iDiffOther;
			isFound = false;
			for (long long j = i + 1; j < nRows; j++) {
				vrow2 = bufInput + j * lenD;
				iChannel2 = (int)(*(vrow2 + 16));
				if ((iChannel2 & 1) == 1) continue;
				if (iChannel2 == 14) continue;
				iBoard2 = (int)(*(vrow2 + 17));
				T2 = *(long long*)vrow2;
				if (T2 > TT) break;
				ind2 = j;
				isFound = true;
				break;
			}
			if (!isFound) continue;

			ind3 = ind1 - 1;

			TT = T1 - iDiffSame;
			isFound = false;
			while (ind3 >= 0) {
				vrow3 = bufInput + ind3 * lenD;
				T3 = *(long long*)vrow3;
				if (T3 < TT) break;
				if (iChannel1 + 1 == (int)(*(vrow3 + 16)) && iBoard1 == (int)(*(vrow3 + 17))) {
					isFound = true;
					break;
				}
				ind3--;
			}

			if (!isFound) {
				ind3 = ind1 + 1;

				TT = T1 + iDiffSame;
				while (ind3 < nRows) {
					vrow3 = bufInput + ind3 * lenD;
					T3 = *(long long*)vrow3;
					if (T3 > TT) break;
					if (iChannel1 + 1 == (int)(*(vrow3 + 16)) && iBoard1 == (int)(*(vrow3 + 17))) {
						isFound = true;
						break;
					}
					ind3++;
				}
			}

			if (!isFound) continue;

			ind4 = ind2 - 1;
			TT = T2 - iDiffSame;
			isFound = false;
			while (ind4 >= 0) {
				vrow4 = bufInput + ind4 * lenD;
				T4 = *(long long*)vrow4;
				if (T4 < TT) break;
				if (iChannel2 + 1 == (int)(*(vrow4 + 16)) && iBoard2 == (int)(*(vrow4 + 17))) {
					isFound = true;
					break;
				}
				ind4--;
			}


			if (!isFound) {
				ind4 = ind2 + 1;
				TT = T2 + iDiffSame;
				while (ind4 < nRows) {
					vrow4 = bufInput + ind4 * lenD;
					T4 = *(long long*)vrow4;
					if (T4 > TT) break;
					if (iChannel2 + 1 == (int)(*(vrow4 + 16)) && iBoard2 == (int)(*(vrow4 + 17))) {
						isFound = true;
						break;
					}
					ind4++;
				}
			}

			if (!isFound) continue;

			vStatus[ind1] = 1;
			vStatus[ind2] = 1;
			vStatus[ind3] = 1;
			vStatus[ind4] = 1;
		}

		for (long long i = 0; i < nDone; i++) {
			if (vStatus[i] == 0) continue;
			vrow = bufInput + i * lenD;
			memcpy(bufOutput + indBuf*lenD, vrow, lenD * sizeof(char));
			indBuf++;
			if (indBuf == BUFF_LENGTH) {
				fwrite(bufOutput, sizeof(char), BUFF_LENGTH* lenD, fo);
				indBuf = 0;
			}
		}

		nBefore = nRows - nDone;
		if (nBefore > 0) {
			for (int i = 0; i < nBefore; i++) {
				vStatus[i] = vStatus[i + nDone];
			}
			memcpy(bufTemp, bufInput + lenD * nDone, sizeof(char) * lenD * nBefore);
			memcpy(bufInput, bufTemp, sizeof(char) * lenD * nBefore);
		}

		startPosition += nDone;

	} while (true);

	if (indBuf > 0) {
		fwrite(bufOutput, sizeof(char), indBuf * lenD, fo);
	}

	return 0;
}

int CoinBinarySingle::startProcessing() {
	if (allocateMemory() != 0) return -5;
	if (processNew() != 0) return -6;
	
	return 0;
}

int printUsage() {
	printf("Input parameters:\n");
	printf("\t1) Input binary file (string)\n");
	printf("\t2) Output binary file (string)\n");
	printf("\t3) timeStamp for the same bar (integer, positive)\n");
	printf("\t4) timeStamp for the different bars (integer, positive)\n\n");
	return 0;
}

int main(int argc, char** argv) {
	int iResult;
	CoinBinarySingle* cbs;
	
	if (argc != 5) {
		printUsage();
		return -99;
	}

	auto start = std::chrono::high_resolution_clock::now();

	cbs = new CoinBinarySingle();

	sprintf(cbs->inputFileName, "%s", argv[1]);
	sprintf(cbs->outputFileName, "%s", argv[2]);
	cbs->iDiffSame = atoi(argv[3]);
	cbs->iDiffOther = atoi(argv[4]);
	
	iResult = cbs->startProcessing();
	delete cbs; cbs = nullptr;

	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsed = finish - start;
	printf("Elapsed Time: %10.2f s", 0.001*elapsed.count());

	return iResult;
}
