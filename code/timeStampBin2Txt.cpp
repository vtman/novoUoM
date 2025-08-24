//D:\NOVO\conData\out\test_merged.bin D:\NOVO\conData\out\test_merged.txt 123456 1000

#include <iostream>
#include <fstream>
#include <chrono>

const int BUFF_LENGTH = 100000;

int printUsage() {
	printf("Input parameters:\n");
	printf("\t1) Input binary file (string)\n");
	printf("\t2) Output text file (string)\n");
	printf("\t3) The first entry (start with 0, int)\n");
	printf("\t4) The number of entries (if < 0, then all remaining, int)\n\n");
	return 0;
}

int main(int argc, char** argv) {
	char inputFile[1000], outputFile[1000], *vloc, *vData;
	long long indFirst, nTotal, lenD, indLast, fileSize, nEntries;
	long long ind, memRead, bufSize, nLeft;
	FILE* fi, * fo;
	
	if (argc != 5) {
		printUsage();
		return -99;
	}

	auto start = std::chrono::high_resolution_clock::now();

	lenD = 18;

	sprintf(inputFile, "%s", argv[1]);
	sprintf(outputFile, "%s", argv[2]);
	#ifdef _WIN32
	indFirst = _atoi64(argv[3]);
	nTotal = _atoi64(argv[4]);
	#else
	indFirst = atoll(argv[3]);
	nTotal = atoll(argv[4]);
	#endif
	
	printf("Input file: %s\n", inputFile);

	fi = fopen(inputFile, "rb");
	if (fi == nullptr) {
		printf("Error: cannot open intput file %s\n", inputFile);
		return -2;
	}

#ifdef _WIN32
	_fseeki64(fi, 0, SEEK_END);
	fileSize = _ftelli64(fi);
	#else
	fseeko64(fi, 0, SEEK_END);
	fileSize = ftello64(fi);
	#endif
	printf("File size: %lli\n", fileSize);

	nEntries = fileSize / lenD;

	printf("Number of entries: %lli\n", nEntries);

	printf("The first entry (input): %lli\n", indFirst);
	if (indFirst < 0) {
		indFirst = 0;
	}

	if (indFirst >= nEntries) {
		printf("Error: the first entry is out of the range\n");
		return -3;
	}

	printf("The first entry (corrected): %lli\n", indFirst);

	printf("The number of output entries (input): %lli\n", nTotal);

	if (nTotal == 0) {
		printf("Warning: nothing to do\n");
		return -4;
	}

	if (nTotal < 0) {
		nTotal = nEntries - indFirst;
	}

	if (nTotal + indFirst > nEntries) {
		nTotal = nEntries - indFirst;
	}

	printf("The number of output entries (corrected): %lli\n", nTotal);

	fo = fopen(outputFile, "w");
	if (fo == nullptr) {
		printf("Error: cannot open output file %s\n", outputFile);
		fclose(fi); fi = nullptr;
		return -4;
	}

	bufSize = lenD * BUFF_LENGTH;
	vData = (char*)malloc(sizeof(char) * bufSize);

	if (vData == nullptr) {
		printf("Error: buffer memory\n");
		fclose(fi); fi = nullptr;
		fclose(fo); fo = nullptr;
		return -5;
	}

	ind = indFirst;
	indLast = indFirst + nTotal;

#ifdef _WIN32
	_fseeki64(fi, lenD * indFirst, SEEK_SET);
	#else
	fseeko64(fi, lenD * indFirst, SEEK_SET);
	#endif

	do {
		nLeft = indLast - ind;
		if (nLeft <= 0) break;
		if (nLeft > BUFF_LENGTH) nLeft = BUFF_LENGTH;
		memRead = nLeft * lenD;
		fread(vData, sizeof(char), memRead, fi);
		for (long long i = 0; i < nLeft; i++) {
			vloc = vData + i * lenD;
			fprintf(fo, "%lli\t%lli\t%lli\t%02i\t%02i\n", ind + i, *(long long*)(vloc), *(long long*)(vloc + 8), (int)(*(vloc + 16)), (int)(*(vloc + 17)));
		}
		ind += nLeft;
	} while (true);

	free(vData); vData = nullptr;
	fclose(fi); fi = nullptr;
	fclose(fo); fo = nullptr;

	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsed = finish - start;
	printf("Elapsed Time: %10.2f s", 0.001 * elapsed.count());

	return 0;
}
