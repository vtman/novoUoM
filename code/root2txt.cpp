//D:\NOVO\conData\out\ts_000206_out.root D:\NOVO\conData\out cc all all 1234 1000

//D:\NOVO\conData\out\ts_000206_out.root D:\NOVO\conData\out cc all timestamp,timestampExtended,channel,board 1234 1000

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TKey.h>
#include <TLeaf.h>
#include <stdio.h>
#include <fstream>
#include <chrono>

const int ROOT_TYPE_BOOL = 0;
const int ROOT_TYPE_CHAR = 1;
const int ROOT_TYPE_UCHAR = 2;
const int ROOT_TYPE_SHORT = 3;
const int ROOT_TYPE_USHORT = 4;
const int ROOT_TYPE_INT = 5;
const int ROOT_TYPE_UINT = 6;
const int ROOT_TYPE_LONG = 7;
const int ROOT_TYPE_ULONG = 8;
const int ROOT_TYPE_LONG64 = 9;
const int ROOT_TYPE_ULONG64 = 10;
const int ROOT_TYPE_FLOAT = 11;
const int ROOT_TYPE_DOUBLE = 12;

class Root2Text {
public:
    Root2Text();
    ~Root2Text();

    int printUsage();
    int startProcessing();
    int processInputParameters();
    int openRootFile();
    int findTrees();
    int processTree();
    int printValue(char* vc, int iType);

    char inputFileName[1000], outputFileName[1000], outputPrefix[1000], outputFolder[1000];
    char* locVariable;
    char inputTreeNames[1000], inputBranchNames[1000];
    char* vTreeNames, * vBranchNames;
    int nSelTrees, nSelBranches;
    int nInputTrees;
    int lenVariable;
    int* vIndTrees, * vIndBranches;
    int nUniqueTrees;

    long long indFirstEntry, nOutEntries;

    bool isAllTrees, isAllBranches;

    TFile* inFile;
    TTree* tree;
    TList* listK;
    FILE* fo, * flog;
    TLeaf* leaf;
};

Root2Text::Root2Text() {
    fo = nullptr;
    flog = nullptr;
    inFile = nullptr;

    vTreeNames = nullptr;
    vBranchNames = nullptr;
}

Root2Text::~Root2Text() {
    if (fo != nullptr) {
        fclose(fo); fo = nullptr;
    }
    if (flog != nullptr) {
        fclose(flog); flog = nullptr;
    }
    if (inFile != nullptr) {
        inFile->Close(); inFile = nullptr;
    }
    if (vTreeNames != nullptr) {
        free(vTreeNames); vTreeNames = nullptr;
    }
    if (vBranchNames != nullptr) {
        free(vBranchNames); vBranchNames = nullptr;
    }
}

int Root2Text::printValue(char* vc, int iType) {
    switch (iType) {
    case ROOT_TYPE_BOOL: fprintf(fo, "%d", *(bool*)vc); break;
    case ROOT_TYPE_CHAR: fprintf(fo, "%c", *(char*)vc); break;
    case ROOT_TYPE_UCHAR: fprintf(fo, "%u", *(unsigned char*)vc); break;
    case ROOT_TYPE_SHORT: fprintf(fo, "%hd", *(short*)vc); break;
    case ROOT_TYPE_USHORT: fprintf(fo, "%hu", *(unsigned short*)vc); break;
    case ROOT_TYPE_INT: fprintf(fo, "%d", *(int*)vc); break;
    case ROOT_TYPE_UINT: fprintf(fo, "%u", *(unsigned int*)vc); break;
    case ROOT_TYPE_LONG: fprintf(fo, "%ld", *(long*)vc); break;
    case ROOT_TYPE_ULONG: fprintf(fo, "%lu", *(unsigned long*)vc); break;
    case ROOT_TYPE_LONG64: fprintf(fo, "%lld", *(long long*)vc); break;
    case ROOT_TYPE_ULONG64: fprintf(fo, "%llu", *(unsigned long long*)vc); break;
    case ROOT_TYPE_FLOAT: fprintf(fo, "%f", *(float*)vc); break;
    case ROOT_TYPE_DOUBLE: fprintf(fo, "%lf", *(double*)vc); break;
    default: fprintf(fo, "_"); break;
    }
    return 0;
}

int Root2Text::printUsage() {
    printf("Usage info:\n");
    printf("\t1) Input ROOT file (string)\n");
    printf("\t2) Output folder (string)\n");
    printf("\t3) Output prefix (string)\n");
    printf("\t4) List of trees to be processed (comma-separated or all for all trees; string)\n");
    printf("\t5) List of branches to be processed (comma-separated or all for all trees; string)\n");
    printf("\t6) Index of the first entry (0 for the first; integer)\n");
    printf("\t7) Number of entries to be processed ( if < 0, then all; integer)\n\n");
    return 0;
}

int Root2Text::processInputParameters() {
    int n, ind1, ind2;
    printf("1) Input ROOT file: %s\n", inputFileName);
    printf("2) Output folder: %s\n", outputFolder);
    printf("3) Output prefix: %s\n", outputPrefix);
    printf("4) List of trees to be processed: %s\n", inputTreeNames);
    printf("5) List of branches to be processed: %s\n", inputBranchNames);
    printf("6) Index of the first entry: %lli\n", indFirstEntry);
    printf("7) Number of entries to be processed: %lli\n\n", nOutEntries);

    sprintf(outputFileName, "%s/%s.log", outputFolder, outputPrefix);
    printf("Log file: %s\n", outputFileName);

    flog = fopen(outputFileName, "w");
    if (flog == nullptr) {
        printf("Error: cannot open file %s\n", outputFileName);
        return -1;
    }

    fprintf(flog, "Input parameters:\n");
    fprintf(flog, "\t1) Input ROOT file: %s\n", inputFileName);
    fprintf(flog, "\t2) Output folder: %s\n", outputFolder);
    fprintf(flog, "\t3) Output prefix: %s\n", outputPrefix);
    fprintf(flog, "\t4) List of trees to be processed: %s\n", inputTreeNames);
    fprintf(flog, "\t5) List of branches to be processed: %s\n", inputBranchNames);
    fprintf(flog, "\t6) Index of the first entry: %lli\n", indFirstEntry);
    fprintf(flog, "\t7) Number of entries to be processed: %lli\n\n", nOutEntries);

    isAllTrees = false;
    if (strcmp(inputTreeNames, "all") == 0) {
        isAllTrees = true;
    }
    else {
        nSelTrees = 1;
        n = strlen(inputTreeNames);
        for (int i = 0; i < n; i++) {
            if (inputTreeNames[i] == ',' || inputTreeNames[i] == ';') nSelTrees++;
        }
        vTreeNames = (char*)malloc(sizeof(char) * 100 * nSelTrees);
        ind1 = 0;
        nSelTrees = 0;
        for (int i = 0; i < n; i++) {
            if (!(inputTreeNames[i] == ',' || inputTreeNames[i] == ';')) continue;
            ind2 = i;
            strncpy(vTreeNames + 100 * nSelTrees, inputTreeNames + ind1, ind2 - ind1);
            vTreeNames[100 * nSelTrees + ind2 - ind1] = '\0';
            nSelTrees++;
            ind1 = i + 1;
        }
        ind2 = n;
        strncpy(vTreeNames + 100 * nSelTrees, inputTreeNames + ind1, ind2 - ind1);
        vTreeNames[100 * nSelTrees + ind2 - ind1] = '\0';
        nSelTrees++;

        fprintf(flog, "Number of selected trees: %i\n", nSelTrees);

        for (int i = 0; i < nSelTrees; i++) {
            fprintf(flog, "\t%s\n", vTreeNames + 100 * i);
        }
        fprintf(flog, "\n");
    }

    isAllBranches = false;
    if (strcmp(inputBranchNames, "all") == 0) {
        isAllBranches = true;
    }
    else {
        nSelBranches = 1;
        n = strlen(inputBranchNames);
        for (int i = 0; i < n; i++) {
            if (inputBranchNames[i] == ',' || inputBranchNames[i] == ';') nSelBranches++;
        }
        vBranchNames = (char*)malloc(sizeof(char) * 100 * nSelBranches);
        ind1 = 0;
        nSelBranches = 0;
        for (int i = 0; i < n; i++) {
            if (!(inputBranchNames[i] == ',' || inputBranchNames[i] == ';')) continue;
            ind2 = i;
            strncpy(vBranchNames + 100 * nSelBranches, inputBranchNames + ind1, ind2 - ind1);
            vBranchNames[100 * nSelBranches + ind2 - ind1] = '\0';
            nSelBranches++;
            ind1 = i + 1;
        }
        ind2 = n;
        strncpy(vBranchNames + 100 * nSelBranches, inputBranchNames + ind1, ind2 - ind1);
        vBranchNames[100 * nSelBranches + ind2 - ind1] = '\0';
        nSelBranches++;

        fprintf(flog, "Number of selected branches: %i\n", nSelBranches);

        for (int i = 0; i < nSelBranches; i++) {
            fprintf(flog, "\t%s\n", vBranchNames + 100 * i);
        }
        fprintf(flog, "\n");
    }

    return 0;
}

int Root2Text::openRootFile() {
    printf("Input File: %s\n", inputFileName);
    fprintf(flog, "Input File: %s\n", inputFileName);
    inFile = TFile::Open(inputFileName);
    if (!inFile || inFile->IsZombie()) {
        printf("Error: cannot open file %s\n", inputFileName);
        fprintf(flog, "Error: cannot open file %s\n", inputFileName);
        return -1;
    }
    printf("The file is openned\n\n");
    fprintf(flog, "The file is openned\n\n");
    return 0;
    return 0;
}

int Root2Text::findTrees() {
    Int_t nKeys, nTreeKeys;
    TKey* key;
    int* vStatus, * vSame;
    int nSame, indCycle, nBest;
    bool isFound;
    
    listK = inFile->GetListOfKeys();
    nKeys = listK->GetEntries();

    printf("Number of keys: %i\n", nKeys);
    fprintf(flog, "Number of keys: %i\n", nKeys);

    if (nKeys == 0) return -1;

    vStatus = (int*)malloc(sizeof(int) * nKeys);
    vSame = (int*)malloc(sizeof(int) * nKeys);

    nTreeKeys = 0;

    for (int i = 0; i < nKeys; i++) {
        vStatus[i] = 0;
        key = (TKey*)listK->At(i);
        if (strcmp(key->GetClassName(), "TTree") != 0) continue;
        vStatus[i] = 1;
        nTreeKeys++;
    }

    printf("Number of keys (TTree): %i\n", nTreeKeys);
    fprintf(flog, "Number of keys (TTree): %i\n", nTreeKeys);

    nInputTrees = 0;

    for (int i = 0; i < nKeys; i++) {
        if (vStatus[i] == 0) continue;
        nSame = 0;
        vSame[nSame] = i;
        nSame++;
        for (int j = i + 1; j < nKeys; j++) {
            if (vStatus[j] == 0) continue;
            if (strcmp(listK->At(j)->GetName(), listK->At(i)->GetName()) != 0) continue;
            vSame[nSame] = j;
            nSame++;
        }
        if (nSame == 1) {
            nInputTrees++;
            continue;
        }
        indCycle = -1;
        nBest = -1;
        for (int j = 0; j < nSame; j++) {
            key = (TKey*)listK->At(vSame[j]);
            if (key->GetCycle() > indCycle) {
                indCycle = key->GetCycle();
                nBest = j;
            }
        }

        for (int j = 0; j < nSame; j++) {
            if (j != nBest) vStatus[vSame[j]] = 0;
        }
        nInputTrees++;
    }

    free(vSame); vSame = nullptr;

    printf("Number of unique trees: %i\n", nInputTrees);
    fprintf(flog, "Number of unique trees: %i\n", nInputTrees);

    if (nInputTrees == 0) {
        free(vStatus); vStatus = nullptr;
        return -2;
    }

    if (!isAllTrees) {
        nUniqueTrees = 0;
        for (int i = 0; i < nKeys; i++) {
            if (vStatus[i] == 0)continue;
            isFound = false;
            for (int j = 0; j < nSelTrees; j++) {
                if (strcmp(vTreeNames + 100 * j, listK->At(i)->GetName()) == 0) {
                    isFound = true;
                    break;
                }
            }
            if (!isFound) {
                vStatus[i] = 0;
                continue;
            }
            nUniqueTrees++;
        }
    }
    else {
        nUniqueTrees = nInputTrees;
    }

    printf("Number of trees (final): %i\n", nUniqueTrees);
    fprintf(flog, "Number of trees (final): %i\n", nUniqueTrees);

    if (nUniqueTrees == 0) {
        free(vStatus); vStatus = nullptr;
        return -3;
    }

    vIndTrees = (int*)malloc(sizeof(int) * nUniqueTrees);

    nUniqueTrees = 0;
    for (int i = 0; i < nKeys; i++) {
        if (vStatus[i] == 0) continue;
        vIndTrees[nUniqueTrees] = i;
        nUniqueTrees++;
    }

    free(vStatus); vStatus = nullptr;
    
    for (int i = 0; i < nUniqueTrees; i++) {
        printf("%i) %s\n", i + 1, listK->At(vIndTrees[i])->GetName());
        fprintf(flog, "%i) %s\n", i + 1, listK->At(vIndTrees[i])->GetName());
    }
    printf("\n");
    fprintf(flog, "\n");

    return 0;
}

int Root2Text::processTree() {
    TObjArray* allBranches;
    TBranch* locBranch;
    int nAllBranches, nB;
    int nCount, iStatus, iType;
    int* vStatus, *vOffset, * vType;
    long long ind;
    bool isTrue;
    Long64_t nEntriesAvail, indFirst, nTotal;

    allBranches = tree->GetListOfBranches();

    nAllBranches = allBranches->GetEntries();

    printf("Tree: %s\n", tree->GetName());
    fprintf(flog, "Tree: %s\n", tree->GetName());
    printf("Number of branches available: %i\n", nAllBranches);
    fprintf(flog, "Number of branches available: %i\n", nAllBranches);

    if (nAllBranches == 0) return 0;

    for (int i = 0; i < nAllBranches; i++) {
        locBranch = (TBranch*)allBranches->At(i);
        fprintf(flog, "%i) %s\n", i + 1, locBranch->GetName());
    }

    printf("\n");
    fprintf(flog, "\n");

    printf("Branches from the input list:\n");
    fprintf(flog, "Branches from the input list:\n");

    if (isAllBranches) {
        nB = nAllBranches;
    }
    else {
        nB = nSelBranches;
    }
    
    vType = (int*)malloc(sizeof(int) * nAllBranches);
    vOffset = (int*)malloc(sizeof(int) * (nAllBranches + 1));
    vStatus = (int*)malloc(sizeof(int) * nB);

    vOffset[0] = 0;

    for (int j = 0; j < nAllBranches; j++) {
        locBranch = (TBranch*)allBranches->At(j);
        leaf = locBranch->GetLeaf(locBranch->GetName());
        vOffset[j + 1] = vOffset[j] + leaf->GetLenType() * leaf->GetNdata();
        vType[j] = -1;
        if (strcmp(leaf->GetTypeName(), "Bool_t") == 0) vType[j] = ROOT_TYPE_BOOL;
        if (strcmp(leaf->GetTypeName(), "Char_t") == 0) vType[j] = ROOT_TYPE_CHAR;
        if (strcmp(leaf->GetTypeName(), "UChar_t") == 0) vType[j] = ROOT_TYPE_UCHAR;
        if (strcmp(leaf->GetTypeName(), "Short_t") == 0) vType[j] = ROOT_TYPE_SHORT;
        if (strcmp(leaf->GetTypeName(), "UShort_t") == 0) vType[j] = ROOT_TYPE_USHORT;
        if (strcmp(leaf->GetTypeName(), "Int_t") == 0) vType[j] = ROOT_TYPE_INT;
        if (strcmp(leaf->GetTypeName(), "UInt_t") == 0) vType[j] = ROOT_TYPE_UINT;
        if (strcmp(leaf->GetTypeName(), "Long_t") == 0) vType[j] = ROOT_TYPE_LONG;
        if (strcmp(leaf->GetTypeName(), "ULong_t") == 0) vType[j] = ROOT_TYPE_ULONG;
        if (strcmp(leaf->GetTypeName(), "Long64_t") == 0) vType[j] = ROOT_TYPE_LONG64;
        if (strcmp(leaf->GetTypeName(), "ULong64_t") == 0) vType[j] = ROOT_TYPE_ULONG64;
        if (strcmp(leaf->GetTypeName(), "Float_t") == 0) vType[j] = ROOT_TYPE_FLOAT;
        if (strcmp(leaf->GetTypeName(), "Double_t") == 0) vType[j] = ROOT_TYPE_DOUBLE;
    }

    nCount = 0;

    for (int i = 0; i < nB; i++) {
        vStatus[i] = -1;
        for (int j = 0; j < nAllBranches; j++) {
            if (isAllBranches && i != j)continue;
            isTrue = false;
            if (isAllBranches) isTrue = true;
            if (!isAllBranches) {
                if (strcmp(allBranches->At(j)->GetName(), vBranchNames + 100 * i) == 0) isTrue = true;
            }
            if (isTrue) {
                locBranch = (TBranch*)allBranches->At(j);
                leaf = locBranch->GetLeaf(locBranch->GetName());
                if (leaf->GetNdata() != 1) continue;
                vStatus[i] = j;
                nCount++;
                fprintf(flog, "%s\t%i\t%i\t%s\t%i\t%i\t%i\n", locBranch->GetName(), j, leaf->GetNdata(), leaf->GetTypeName(), leaf->GetLenType(), vOffset[j], vOffset[j+1]);
                break;
            }
        }
        
    }

    fprintf(flog, "Total number of variables: %i\n", nCount);

    lenVariable = vOffset[nAllBranches];
    fprintf(flog, "Length of storage: %i\n", lenVariable);

    if (lenVariable == 0) {
        free(vStatus); vStatus = nullptr;
        free(vOffset); vOffset = nullptr;
        return 0;
    }

    locVariable = (char*)malloc(sizeof(char) * lenVariable);

    tree->SetMakeClass(1);

    for (int i = 0; i < nB; i++) {
        if (vStatus[i] < 0) continue;
        printf("%s\n", allBranches->At(vStatus[i])->GetName());
        tree->SetBranchAddress(allBranches->At(vStatus[i])->GetName(), locVariable + vOffset[vStatus[i]]);
    }

    sprintf(outputFileName, "%s/%s_%s.txt", outputFolder, outputPrefix, tree->GetName());
    fo = fopen(outputFileName, "wb");
    if (fo == nullptr) {
        printf("Error: cannot open file %s\n", outputFileName);
        free(locVariable); locVariable = nullptr;
        free(vStatus); vStatus = nullptr;
        free(vOffset); vOffset = nullptr;
        return -1;
    }

    nEntriesAvail = tree->GetEntries();
    fprintf(flog, "Number of entries available: %lli\n", nEntriesAvail);

    if (indFirstEntry < 0) {
        indFirst = 0;
    }
    else {
        indFirst = indFirstEntry;
    }

    if (indFirst >= nEntriesAvail) {
        fclose(fo); fo = nullptr;
        free(locVariable); locVariable = nullptr;
        free(vStatus); vStatus = nullptr;
        free(vOffset); vOffset = nullptr;
    }

    nTotal = nEntriesAvail - indFirst;
    if (nOutEntries >= 0) {
        if (nTotal > nOutEntries) nTotal = nOutEntries;
    }
    fprintf(flog, "Number of entries to write: %lli\n", nTotal);
    
    for (long long i = 0; i < nTotal; i++) {
        ind = indFirst + i;
        tree->GetEntry(ind);
        fprintf(fo, "%lli", ind);
        for (int j = 0; j < nB; j++) {
            fprintf(fo, "\t");
            iStatus = vStatus[j];
            if (iStatus < 0)continue;
            iType = vType[iStatus];
            printValue(locVariable + vOffset[iStatus], iType);
        }
        fprintf(fo, "\n");
    }

    fprintf(flog, "\n");
    printf("\n");

    fclose(fo); fo = nullptr;
    free(locVariable); locVariable = nullptr;
    free(vStatus); vStatus = nullptr;
    free(vOffset); vOffset = nullptr;
    return 0;
}

int Root2Text::startProcessing() {
    if (processInputParameters() != 0) return -1;
    if (openRootFile() != 0) return -2;
    if (findTrees() != 0) return -3;

    for (int i = 0; i < nUniqueTrees; i++) {
        tree = (TTree*)inFile->Get(listK->At(vIndTrees[i])->GetName());
        processTree();
    }
    return 0;
}

int main(int argc, char** argv) {
    int iResult;
    Root2Text* r2t;

    r2t = new Root2Text();
    if (argc != 8) {
        r2t->printUsage();
        delete r2t; r2t = nullptr;
        return -99;
    }

    auto start = std::chrono::high_resolution_clock::now();

    sprintf(r2t->inputFileName, "%s", argv[1]);
    sprintf(r2t->outputFolder, "%s", argv[2]);
    sprintf(r2t->outputPrefix, "%s", argv[3]);
    sprintf(r2t->inputTreeNames, "%s", argv[4]);
    sprintf(r2t->inputBranchNames, "%s", argv[5]);
    #ifdef _WIN32
    r2t->indFirstEntry = _atoi64(argv[6]);
    r2t->nOutEntries = _atoi64(argv[7]);
    #else
    r2t->indFirstEntry = atoll(argv[6]);
    r2t->nOutEntries = atoll(argv[7]);
    #endif

    iResult = r2t->startProcessing();

    delete r2t; r2t = nullptr;

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = finish - start;
    printf("Elapsed Time: %10.2f s", 0.001 * elapsed.count());

    return iResult;
}
