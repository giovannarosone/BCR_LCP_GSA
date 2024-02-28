#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <assert.h>
#include <sys/stat.h>
#include <math.h>
#include <string.h>   

#define SIZEBUFFER 1024

int main (int argn, char** argv) {

    char* fileOutput=argv[1];
    char *bufferBWT = new char[SIZEBUFFER];
    char *bufferSAP = new char[SIZEBUFFER];
    std::string inserted_symbols("");

    bool found_sap=0;
    ulong number_sap=0;
	ulong sap_length=0;
    ulong total_sap_length=0;
    ulong total_diff_char=0;
	char symbol;
    char lastSymbol;
    ulong total_run_sap=0;

    ulong lung = strlen(fileOutput);
	char *fileOutRes = new char[lung+100];
	char *fnBWT = new char[lung+100];
    sprintf (fileOutRes,"%s%s",fileOutput,".txt");
    sprintf (fnBWT,"%s%s",fileOutput,".ebwt");
    FILE *InFileBWT = fopen(fnBWT, "rb");
    if (InFileBWT==NULL) {
        std::cerr << "printOutput: Entire BWT file: Error opening "  << fnBWT <<  "." << std::endl;
        exit (EXIT_FAILURE);
    }

    char *fnSAP = new char[lung+100];
    sprintf (fnSAP,"%s%s",fileOutput,".bwt.red_sap");
    FILE *InFileSAP = fopen(fnSAP, "rb");
    if (InFileSAP==NULL) {
        std::cerr << "printOutput: Entire SAP file: Error opening " << fnSAP <<  "." << std::endl;
        exit (EXIT_FAILURE);
    }

    ulong numcharBWT;
    ulong numcharSAP;

    while ( (numcharBWT = fread(bufferBWT,sizeof(char),SIZEBUFFER,InFileBWT)) && (numcharBWT>0) ) {

        numcharSAP = fread(bufferSAP,sizeof(char),SIZEBUFFER,InFileSAP);
		assert(numcharSAP == numcharBWT);

        for (unsigned int i=0; i < numcharBWT; i++) {
            if (bufferSAP[i] == '1') {
                if(found_sap == 0) {
                    number_sap++;
                    found_sap=1;
                    inserted_symbols+=symbol;
                    lastSymbol=symbol;
                    total_run_sap++;
                }
                sap_length++;

                if(bufferBWT[i]!=lastSymbol) {
                    total_run_sap++;
                    lastSymbol=bufferBWT[i];
                }

                if(inserted_symbols.find(bufferBWT[i]) == std::string::npos) {
                    inserted_symbols+=bufferBWT[i];
                }
            }
             else {
                if(found_sap == 1) {
                    found_sap=0;
                    total_sap_length=total_sap_length+sap_length+1;
                    sap_length=0;
                    total_diff_char=total_diff_char+inserted_symbols.length();
                    inserted_symbols="";
                }
                symbol=bufferBWT[i];
                lastSymbol=bufferBWT[i];
            }
        }
    }

    if(found_sap == 1) {
        total_sap_length=total_sap_length+sap_length+1;
        total_diff_char=total_diff_char+inserted_symbols.length();
    }

    fclose(InFileBWT);
    fclose(InFileSAP);

	std::cerr << "total sap interval of type two: " << number_sap << "\n";
	std::cerr << "total length of sap intervals of type two: " << total_sap_length << "\n";
	std::cerr << "total distinctive chracters in sap intervals of type two: " << total_diff_char << "\n";
    std::cerr << "total runs in sap intervals of type two: " << total_run_sap << "\n";
    return 0;
}
