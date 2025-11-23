//D:\NOVO\2025_02_Cf252\test\temp\153_cleanPairs.bin D:\NOVO\2025_02_Cf252\test\temp\153_cleanPairs.root 2 1

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <chrono>

const long long BUFF_LENGTH = 32000;
const long long LEN_RECORD = 22;
const long long LEN_BLOCK = 2 * LEN_RECORD;
const float RAT_VALUE = 1.0 / 1024.0;
const long long NBITS = 10ULL;

class BlockData {
public:
	BlockData();
	~BlockData();

	int allocateMemory();
	int convertInput2Output();

	int nBoards;
	int nChannels;
	int nHits;
	int nBlocks;
	int nArrays;
	int nRows;

	char* vInput, *vOutput;
	char** pvEshort, ** pvElong, ** pvTime, ** pvTimestamp, ** pvFlag;
	char** pvArrays;
	int* vStartPoistion;
};

class Bin2Root {
public:
	Bin2Root();
	~Bin2Root();
	
	int openInputFile();
	int openOutputFile();
	int startProcessing();
	int createBlocks();
	int allocateRootArrays();
	int processBlocks();
	int createBranches();

	int printBuffer();

	char inputFileName[1000], outputFileName[1000];

	long long inputBlockSize, nBlocks;
	long long nHits, nBoards;
	long long fileSize, chunkSize, bufSize;

	int nThreads, nChannels;
	char* vData, *vBuffer;

	TFile* outputTFile;
	TTree* tree;

	FILE *fi;

	BlockData* BD;
};

Bin2Root::Bin2Root() {
	fi = nullptr;
	vBuffer = nullptr;
	vData = nullptr;

	outputTFile = nullptr;
	BD = nullptr;
	nThreads = 0;
	nChannels = 0;
}

Bin2Root::~Bin2Root() {
	if (BD != nullptr) { delete BD; BD = nullptr; }
	if (vData != nullptr) { free(vData); vData = nullptr; }
	if (vBuffer != nullptr) { free(vBuffer); vBuffer = nullptr; }
	if (fi != nullptr) { fclose(fi); fi = nullptr;}
}

int Bin2Root::openInputFile() {
	printf("Input file: %s\n", inputFileName);
	fi = fopen(inputFileName, "rb");
	
	if (fi == nullptr) {
		printf("Error: cannot open the file\n");
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

	inputBlockSize = nHits * LEN_BLOCK;
	printf("Input block size: %lli\n", inputBlockSize);

	nBlocks = fileSize / inputBlockSize;
	printf("Number of blocks: %lli\n", nBlocks);

	if (nBlocks * inputBlockSize != fileSize) {
		printf("Warning: the file size is not a multiple of block size");
	}

	return 0;
}

int Bin2Root::openOutputFile() {
	printf("Output file: %s\n", outputFileName);
	outputTFile = new TFile(outputFileName, "RECREATE", "", 5);
	if (outputTFile == nullptr) {
		printf("Error: null pointer for the output file\n");
		return -1;
	}
	if (outputTFile->IsZombie()) {
		printf("Error: zombie\n");
		return -2;
	}
	if (!outputTFile->IsOpen()) {
		printf("Error: file is not open\n");
		return -3;
	}
	if (!outputTFile->IsWritable()) {
		printf("Error: file is not writable\n");
		return -4;
	}
	return 0;
}

BlockData::BlockData() {
	pvElong = nullptr;
	pvEshort = nullptr;
	pvFlag = nullptr;
	pvTime = nullptr;
	pvTimestamp = nullptr;
	vInput = nullptr;
	vOutput = nullptr;

	pvArrays = nullptr;
	vStartPoistion = nullptr;

	nArrays = 0;
	nChannels = 0;
}

int BlockData::allocateMemory() {
	nChannels = 14 * nBoards;
	nArrays = nChannels * 5;

	pvArrays = (char**)malloc(sizeof(char*) * nArrays);
	vStartPoistion = (int*)malloc(sizeof(int) * (nArrays + 1));

	pvElong = (char**)malloc(sizeof(char*) * nChannels);
	pvEshort = (char**)malloc(sizeof(char*) * nChannels);
	pvFlag = (char**)malloc(sizeof(char*) * nChannels);
	pvTime = (char**)malloc(sizeof(char*) * nChannels);
	pvTimestamp = (char**)malloc(sizeof(char*) * nChannels);

	vStartPoistion[0] = 0;
	nArrays = 0;

	for (int i = 0; i < nChannels; i++) {
		pvElong[i] = (char*)malloc(sizeof(int16_t) * BUFF_LENGTH);
		pvArrays[nArrays] = pvElong[i];
		pvEshort[i] = (char*)malloc(sizeof(int16_t) * BUFF_LENGTH);
		pvArrays[nArrays + 1] = pvEshort[i];
		pvFlag[i] = (char*)malloc(sizeof(char) * BUFF_LENGTH);
		pvArrays[nArrays + 2] = pvFlag[i];
		pvTime[i] = (char*)malloc(sizeof(float) * BUFF_LENGTH);
		pvArrays[nArrays + 3] = pvTime[i];
		pvTimestamp[i] = (char*)malloc(sizeof(int64_t) * BUFF_LENGTH);
		pvArrays[nArrays + 4] = pvTimestamp[i];

		vStartPoistion[nArrays + 1] = vStartPoistion[nArrays] + 2; nArrays++;//Elong
		vStartPoistion[nArrays + 1] = vStartPoistion[nArrays] + 2; nArrays++;//Eshort
		vStartPoistion[nArrays + 1] = vStartPoistion[nArrays] + 1; nArrays++;//Flag
		vStartPoistion[nArrays + 1] = vStartPoistion[nArrays] + 4; nArrays++;//Time
		vStartPoistion[nArrays + 1] = vStartPoistion[nArrays] + 8; nArrays++;//Timestamp

		if (pvElong[i] == nullptr || pvEshort[i] == nullptr || pvFlag[i] == nullptr || pvTime[i] == nullptr || pvTimestamp[i] == nullptr) {
			printf("Error: memory allocation (output blocks)\n");
			return -1;
		}
	}

	vInput = (char*)malloc(sizeof(char) * BUFF_LENGTH * nHits * LEN_BLOCK);
	if (vInput == nullptr) {
		printf("Error: memory allocation (input block)\n");
		return -2;
	}

	vOutput = (char*)malloc(sizeof(char) * BUFF_LENGTH * nChannels * 17);
	if (vOutput == nullptr) {
		printf("Error: memory allocation (output block)\n");
		return -3;
	}

	return 0;
}

BlockData::~BlockData() {
	if (vInput != nullptr) {free(vInput); vInput = nullptr;}
	if (vOutput != nullptr) {free(vOutput); vOutput = nullptr;}
	if (vStartPoistion != nullptr) {free(vStartPoistion); vStartPoistion = nullptr;}
	if (pvArrays != nullptr) {free(pvArrays); pvArrays = nullptr;}
	if (pvElong != nullptr) {
		for (int i = 0; i < nChannels; i++) {
			if (pvElong[i] != nullptr) { free(pvElong[i]); pvElong[i] = nullptr;}
		}
		free(pvElong); pvElong = nullptr;
	}
	if (pvEshort != nullptr) {
		for (int i = 0; i < nChannels; i++) {
			if (pvEshort[i] != nullptr) { free(pvEshort[i]); pvEshort[i] = nullptr;}
		}
		free(pvEshort); pvEshort = nullptr;
	}
	if (pvFlag != nullptr) {
		for (int i = 0; i < nChannels; i++) {
			if (pvFlag[i] != nullptr) { free(pvFlag[i]); pvFlag[i] = nullptr;}
		}
		free(pvFlag); pvFlag = nullptr;
	}
	if (pvTime != nullptr) {
		for (int i = 0; i < nChannels; i++) {
			if (pvTime[i] != nullptr) { free(pvTime[i]); pvTime[i] = nullptr;}
		}
		free(pvTime); pvTime = nullptr;
	}
	if (pvTimestamp != nullptr) {
		for (int i = 0; i < nChannels; i++) {
			if (pvTimestamp[i] != nullptr) {free(pvTimestamp[i]); pvTimestamp[i] = nullptr;}
		}
		free(pvTimestamp); pvTimestamp = nullptr;
	}
}

int BlockData::convertInput2Output(){
	int lenO, lenI, m, iCB, iCBnew;
	int iChannel, iBoard;
	char* vLocI, * vLocO, *vi, *vo;

	m = 2 * nHits;
	lenI = LEN_BLOCK * nHits;
	lenO = vStartPoistion[nArrays];

	memset(vOutput, 0, lenO * nRows);

	for (int i = 0; i < nRows; i++) {
		vLocI = vInput + lenI * i;
		vLocO = vOutput + lenO * i;
		for (int ih = 0; ih < m; ih++) {
			vi = vLocI + LEN_RECORD * ih;
			iCB = vi[16];
			iChannel = iCB & 15;
			iBoard = iCB >> 4;
			if (iBoard >= nBoards) continue;
			iCBnew = iBoard * 14 + iChannel;
			vo = vLocO + iCBnew * 17;

			//Eshort
			*(uint16_t*)(vo + 0) = *(uint16_t*)(vi + 19);
			//Elong
			*(uint16_t*)(vo + 2) = *(uint16_t*)(vi + 17);
			//time
			*(float*)(vo + 4) = (float)(1023 & *(uint32_t*)(vi + 0)) * RAT_VALUE;
			//timestamp
			*(long long*)(vo + 8) = *(uint64_t*)(vi + 0) >> NBITS;
			//flags
			vo[16] = vi[21];
		}
	}
	return 0;
}

int Bin2Root::allocateRootArrays() {
	chunkSize = LEN_BLOCK * nHits * BUFF_LENGTH;
	printf("Chunk size: %lli\n", chunkSize);

	vData = (char*)malloc(sizeof(char) * chunkSize);
	if (vData == nullptr) {
		printf("Error: chunk data\n");
		return -1;
	}

	bufSize = nChannels * 17;

	vBuffer = (char*)malloc(sizeof(char) * bufSize);

	return 0;
}

int Bin2Root::createBlocks() {
	BD = new BlockData();
	BD->nBoards = nBoards;
	BD->nHits = nHits;
	if (BD->allocateMemory() != 0) {
		printf("Error: block data\n");
		return -1;
	}

	nChannels = 14 * nBoards;
	return 0;
}

int Bin2Root::processBlocks() {
	long long chunkLeft, nChunks, nRows;
	
	nRows = BUFF_LENGTH;

	nChunks = fileSize / chunkSize;
	if (nChunks * chunkSize < fileSize) nChunks++;

	printf("Number of chunks: %lli\n", nChunks);

	chunkLeft = chunkSize;
	for (long long ic = 0; ic < nChunks; ic++) {
		printf("IC: %lli\n", ic);
		if (ic == nChunks - 1) {
			chunkLeft = fileSize - (nChunks - 1) * chunkSize;
			nRows = chunkLeft / (LEN_BLOCK * nHits);
		}

#ifdef _WIN32
		_fseeki64(fi, ic * chunkSize, SEEK_SET);
#else
		fseeko64(fi, ic * chunkSize, SEEK_SET);
#endif
		fread(vData, sizeof(char), chunkLeft, fi);

		memcpy(BD->vInput, vData, chunkLeft);
		BD->nRows = nRows;
		BD->convertInput2Output();

		for (long long i = 0; i < nRows; i++) {
			memcpy(vBuffer, BD->vOutput + bufSize * i, bufSize * sizeof(char));
			tree->Fill();
		}

	}
	tree->Write();
	outputTFile->Close();
	delete outputTFile; outputTFile = nullptr;

	return 0;
}

int Bin2Root::printBuffer() {
	int ipos;
	ipos = 0;
	for (int k = 0; k < nChannels; k++) {
		printf("%i\t", k);
		printf("%i\t", *(int16_t*)(vBuffer + ipos)); ipos += 2;
		printf("%i\t", *(int16_t*)(vBuffer + ipos)); ipos += 2;
		printf("%f\t", *(float*)(vBuffer + ipos)); ipos += 4;
		printf("%lli\t", *(long long*)(vBuffer + ipos)); ipos += 8;
		printf("%i", *(vBuffer + ipos)); ipos++;
		printf("\n");
	}
	return 0;
}

int Bin2Root::createBranches() {
	char bName[100], bPrefix[100];
	int ipos;

	tree = new TTree("data", "coincidence data");

	ipos = 0;

	for (int ib = 0; ib < nBoards; ib++) {
		for (int ic = 0; ic < 14; ic++) {
			sprintf(bPrefix, "dt%c_ch%i_", 'a' + ib, ic);
			
			sprintf(bName, "%sEshort", bPrefix);
			tree->Branch(bName, vBuffer + ipos, "Eshort/s");
			ipos += 2;

			sprintf(bName, "%sElong", bPrefix);
			tree->Branch(bName, vBuffer + ipos, "Elong/s");
			ipos += 2;

			sprintf(bName, "%stime", bPrefix);
			tree->Branch(bName, vBuffer + ipos, "time/F");
			ipos += 4;

			sprintf(bName, "%stimestamp", bPrefix);
			tree->Branch(bName, vBuffer + ipos, "timestamp/l");
			ipos += 8;

			sprintf(bName, "%sisClipped", bPrefix);
			tree->Branch(bName, vBuffer + ipos, "isClipped/O");
			ipos += 1;
		}
	}
	return 0;
}

int Bin2Root::startProcessing() {
	if (openInputFile() != 0) return -1;
	if (openOutputFile() != 0) return -2;
	if (createBlocks() != 0) return -3;
	if (allocateRootArrays() != 0) return -4;
	if (createBranches() != 0) return -5;
	if (processBlocks() != 0) return -6;
	return 0;
}

int printUsage() {
	printf("Input parameters:\n");
	printf("\t1) Input binary file (string)\n");
	printf("\t2) Output ROOT file (string)\n");
	printf("\t3) Number of boards (integer, positive)\n");
	printf("\t4) Number of hits (integer, positive)\n\n");
	return 0;
}

int main(int argc, char** argv) {
	int iResult;
	Bin2Root* b2r;
	
	if (argc != 5) {
		printUsage();
		return -99;
	}

	/*
	TFile* tf = TFile::Open("D:\\NOVO\\2025_02_Cf252\\test\\temp\\153_coinTwo.root", "READ");
	TTree* tree = (TTree*)tf->Get("data");
	tree->Scan("*", "", "", 100);
	return 0;
	*/

	auto start = std::chrono::high_resolution_clock::now();

	b2r = new Bin2Root();

	sprintf(b2r->inputFileName, "%s", argv[1]);
	sprintf(b2r->outputFileName, "%s", argv[2]);
	b2r->nBoards = atoi(argv[3]);
	b2r->nHits = atoi(argv[4]);
	
	iResult = b2r->startProcessing();
	delete b2r; b2r = nullptr;

	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsed = finish - start;
	printf("Elapsed Time: %10.2f s", 0.001*elapsed.count());

	return iResult;
}
