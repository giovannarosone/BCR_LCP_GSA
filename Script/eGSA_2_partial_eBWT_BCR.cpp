//
//  main.cpp
//  ScriptBCR
//
//  Created by Giovanna on 23/12/18.
//  Copyright © 2018 Giovanna. All rights reserved.
//

//Computes the BCR partial files starting to another EGSA
//See output EGSA of Felipe Louze


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <assert.h>
#include <sys/stat.h>
#include <math.h>
#include "Parameters.h"
#include "Tools.h"

using namespace std;


//See Felipe Louza tool: eGSA

typedef unsigned int int_text;
typedef unsigned int int_suff;     //-2^31 to 2^31
typedef unsigned int int_lcp;
typedef unsigned char int8; //0 to 2^8

#define PREFIX_SIZE         10

typedef struct{
    dataTypeNChar    text;
    dataTypelenSeq    suff;
    dataTypelenSeq     lcp;
    
    uchar        bwt;
} t_GSA;

bool findLengthNseq( const string& input, const string& fileOutput);

int main(int argc, char **argv) {
    
    if( argc != 3 )
    {
        std::cerr << "usage: " << argv[0] << " input" << std::endl;
        std::cerr << "where:" << std::endl;
        std::cerr << "  input is the filename without 0.gesa " << std::endl;
        std::cerr << "  output is the filename for BCR partial files " << std::endl;
        exit(1);
    }
    
    string fileInput = argv[1];
    string fileOutput = argv[2];
    
    
    std::cerr << "ReadEGSA: File: " << fileInput << std::endl;
    findLengthNseq( fileInput, fileOutput);
    
    //Open EGSA file
    FILE         *f_ESA;            // pointer to the ESA input file
    string fnEGSA;
    fnEGSA = fileInput + "." + "0" + ".gesa\0";
    std::cerr << "ReadEGSA: File: " << fnEGSA << std::endl;
    
    f_ESA = fopen(fnEGSA.c_str(), "rb");
    if (f_ESA == NULL) {
        std::cerr << "readEGSA: Error opening: " << fnEGSA << std::endl;
        exit (EXIT_FAILURE);
    }
    fseek(f_ESA , 0, SEEK_SET);
    t_GSA GSA;
    
    const char *ext = ".aux";
    
    dataTypeNChar numcharBWT, numcharPairSA, numcharLCP, numEle=0, numRead=0, numWrite=0;
    
    uchar c;
    
    dataTypeNChar freq[256];  //contains the distribution of the symbols. It is useful only for testing. It depends on the #characters
    dataTypedimAlpha alpha[SIZE_ALPHA]; //Corresponding between the alphabet, the piles and tableOcc
    dataTypedimAlpha sizeAlpha;  //number of the different symbols in the input texts
    vector<dataTypedimAlpha> alphaInverse;  //Corresponding between alpha[i] and the symbol as char
    dataTypeNChar** tableOcc; //contains the number of occurrences of each symbol
    
    for (dataTypedimAlpha z = 0 ; z < SIZE_ALPHA-1; z++)
        freq[z]=0;
    freq[SIZE_ALPHA-1]=0;
    
    
    //First reading in order to find the alphabet
    numRead=fread(&GSA.text, sizeof(int_text), 1, f_ESA);
    if (numRead > 0) {
        fread(&GSA.suff, sizeof(int_suff), 1, f_ESA);
        fread(&GSA.lcp, sizeof(int_lcp), 1, f_ESA);
        fread(&GSA.bwt, sizeof(int8), 1, f_ESA);
    }
    
    dataTypelenSeq minLCP=GSA.lcp, maxLCP=GSA.lcp;
    if (GSA.bwt == '\0')
        c=TERMINATE_CHAR;
    else
        c=GSA.bwt;
    freq[(unsigned int)(c)]=1;
    
    while ( fread(&GSA.text, sizeof(int_text), 1, f_ESA) > 0 )   {
        //std::cerr << "bwt\tlcp\tpos\tnumSeq\tSA\tQS\n";
        fread(&GSA.suff, sizeof(int_suff), 1, f_ESA);
        fread(&GSA.lcp, sizeof(int_lcp), 1, f_ESA);
        fread(&GSA.bwt, sizeof(int8), 1, f_ESA);
        //std::cerr << "-->" << numEle << "\t" << GSA.text << "\t" << GSA.suff <<  "\t(" << GSA.lcp <<  "," << minLCP << ")\t" << c << "\n";
        if (GSA.bwt == '\0')
            c=TERMINATE_CHAR;
        else
            c=GSA.bwt;

        freq[(unsigned int)(c)]++;
            
        if (minLCP > GSA.lcp)
                minLCP=GSA.lcp;
        if (maxLCP < GSA.lcp)
                maxLCP=GSA.lcp;
            //std::cerr << "-->" << numEle << "\t" << GSA.text << "\t" << GSA.suff <<  "\t(" << GSA.lcp <<  "," << minLCP << ")\t" << c << "\n";
            
        numEle++;
    }//end-while
    
    sizeAlpha=0;
    for (dataTypedimAlpha i = 0; i < SIZE_ALPHA-1; ++i)
        if (freq[i] > 0) {
            sizeAlpha++;
        }
    if (freq[SIZE_ALPHA-1] > 0) {
        sizeAlpha++;
    }
    
    /*
    alphaInverse.resize(sizeAlpha);
    sizeAlpha=0;
    for (dataTypedimAlpha i = 0; i < SIZE_ALPHA-1; ++i)
        if (freq[i] > 0) {
            alpha[i] = sizeAlpha;
            alphaInverse[sizeAlpha]=i;
            sizeAlpha++;
        }
    if (freq[SIZE_ALPHA-1] > 0) {
        alpha[SIZE_ALPHA-1] = sizeAlpha;
        alphaInverse[sizeAlpha]=SIZE_ALPHA-1;
        sizeAlpha++;
    }
    */
    
    std::cerr << "We supposed that the " <<numEle << " symbols in the input file are (ASCII, char, freq, code):\n";
    for (dataTypedimAlpha i = 0; i < SIZE_ALPHA-1; ++i)
        if (freq[i] > 0)
            std::cerr << (unsigned int)i << " " << i << " " << freq[i] << "\n";
    if (freq[SIZE_ALPHA-1] > 0)
        std::cerr << (unsigned int)SIZE_ALPHA-1 << " " << SIZE_ALPHA-1 << " " << freq[SIZE_ALPHA-1] << " " << "\n";
    
    

    fseek(f_ESA , 0, SEEK_SET);
    
    //Open BCR partial files
    
    
    std::cerr << "Build the partial BCR files for eBWT/LCP/DA/SA.\n";
    char *fnBWT = new char[fileOutput.size()+100];
    static FILE *OutFileBWT;
    
    #if BUILD_LCP == 1
        char *fnLCP = new char[fileOutput.size()+100];
        static FILE *OutFileLCP;
    #endif
    #if BUILD_DA == 1
        char *fnDA = new char[fileOutput.size()+100];
        static FILE *OutFileDA;
    #endif
    #if BUILD_SA == 1
        char *fnSA = new char[fileOutput.size()+100];
        static FILE *OutFileSA;
    #endif

    numEle=0;
    numRead=0;
    numWrite=0;
    minLCP=0;
    //std::cout << "nSeq" << "\t" << "Pos" << "\t" <<  "LCP" << "\t" <<  "minLCP" << "\t" <<  "BWT" << "\t" << "\n";
    
    for (dataTypedimAlpha currentPile = 0 ; currentPile < SIZE_ALPHA-1; currentPile++) {
        
        if (freq[currentPile] > 0) {
            sprintf (fnBWT,"%s%s%d",fileOutput.c_str(),".ebwt_",currentPile);
            OutFileBWT = fopen(fnBWT, "wb");
            if (OutFileBWT==NULL) {
                std::cerr << "Error opening: " << fnBWT << std::endl;
                exit (EXIT_FAILURE);
            }
            #if BUILD_LCP == 1
                sprintf (fnLCP,"%s%s%d",fileOutput.c_str(),".lcp_",currentPile);
                OutFileLCP = fopen(fnLCP, "wb");
                if (OutFileLCP==NULL) {
                    std::cerr << "Error opening: " << fnLCP << std::endl;
                    exit (EXIT_FAILURE);
                }
            #endif
            #if BUILD_DA == 1
                sprintf (fnDA,"%s%s%d",fileOutput.c_str(),".da_",currentPile);
                OutFileDA = fopen(fnDA, "wb");
                if (OutFileDA==NULL) {
                    std::cerr << "Error opening: " << fnDA << std::endl;
                    exit (EXIT_FAILURE);
                }
            #endif
            #if BUILD_SA == 1
                sprintf (fnSA,"%s%s%d",fileOutput.c_str(),".sa_",currentPile);
                OutFileSA = fopen(fnSA, "wb");
                if (OutFileSA==NULL) {
                    std::cerr << "Error opening: " << fnSA << std::endl;
                    exit (EXIT_FAILURE);
                }
            #endif
            //std::cerr << "bwt\tlcp\tpos\tnumSeq\tSA\tQS\n";
            
            for (dataTypeNChar j = 0 ; j < freq[currentPile]; j++) {

                numRead=fread(&GSA.text, sizeof(int_text), 1, f_ESA);
                fread(&GSA.suff, sizeof(int_suff), 1, f_ESA);
                fread(&GSA.lcp, sizeof(int_lcp), 1, f_ESA);
                fread(&GSA.bwt, sizeof(int8), 1, f_ESA);
                
                if (GSA.bwt == '\0')
                    c=TERMINATE_CHAR;
                else
                    c=GSA.bwt;
                
               fwrite(&c, sizeof(uchar), 1, OutFileBWT);
                #if BUILD_LCP == 1
                    fwrite (&GSA.lcp, sizeof(dataTypelenSeq), 1 , OutFileLCP);
                #endif
                #if BUILD_DA == 1
                    fwrite (&GSA.text, sizeof(dataTypeNSeq), 1 , OutFileDA);
                #endif
                #if BUILD_SA == 1
                    fwrite (&GSA.suff, sizeof(dataTypelenSeq), 1 , OutFileSA);
                #endif
            
                
                //std::cerr << "-->" << numEle << "\t" << GSA.text << "\t" << GSA.suff <<  "\t(" << GSA.lcp <<  "," << minLCP << ")\t" << c << "\n";
                
                numEle++;
            }  //end-for
        } //end-if
    }  //end-for
    
    if (freq[SIZE_ALPHA-1] > 0) {
        sprintf (fnBWT,"%s%s%d",fileOutput.c_str(),".ebwt_",SIZE_ALPHA-1);
        OutFileBWT = fopen(fnBWT, "wb");
        if (OutFileBWT==NULL) {
            std::cerr << "Error opening: " << fnBWT << std::endl;
            exit (EXIT_FAILURE);
        }
        #if BUILD_LCP == 1
            sprintf (fnLCP,"%s%s%d",fileOutput.c_str(),".lcp_",SIZE_ALPHA-1);
            OutFileLCP = fopen(fnLCP, "wb");
            if (OutFileLCP==NULL) {
                std::cerr << "Error opening: " << fnLCP << std::endl;
                exit (EXIT_FAILURE);
            }
        #endif
        #if BUILD_DA == 1
            sprintf (fnDA,"%s%s%d",fileOutput.c_str(),".da_",SIZE_ALPHA-1);
            OutFileDA = fopen(fnDA, "wb");
            if (OutFileDA==NULL) {
                std::cerr << "Error opening: " << fnDA << std::endl;
                exit (EXIT_FAILURE);
            }
        #endif
        #if BUILD_SA == 1
            sprintf (fnSA,"%s%s%d",fileOutput.c_str(),".sa_",SIZE_ALPHA-1);
            OutFileSA = fopen(fnSA, "wb");
            if (OutFileSA==NULL) {
                std::cerr << "Error opening: " << fnSA << std::endl;
                exit (EXIT_FAILURE);
            }
        #endif
        //std::cerr << "bwt\tlcp\tpos\tnumSeq\tSA\tQS\n";
        
        for (dataTypeNChar j = 0 ; j < freq[SIZE_ALPHA-1]; j++) {
            
            numRead=fread(&GSA.text, sizeof(int_text), 1, f_ESA);
            fread(&GSA.suff, sizeof(int_suff), 1, f_ESA);
            fread(&GSA.lcp, sizeof(int_lcp), 1, f_ESA);
            fread(&GSA.bwt, sizeof(int8), 1, f_ESA);
            
            if (GSA.bwt == '\0')
                c=TERMINATE_CHAR;
            else
                c=GSA.bwt;
            
            fwrite(&c, sizeof(uchar), 1, OutFileBWT);
            #if BUILD_LCP == 1
                fwrite (&GSA.lcp, sizeof(dataTypelenSeq), 1 , OutFileLCP);
            #endif
            #if BUILD_DA == 1
                fwrite (&GSA.text, sizeof(dataTypeNSeq), 1 , OutFileDA);
            #endif
            #if BUILD_SA == 1
                fwrite (&GSA.suff, sizeof(dataTypelenSeq), 1 , OutFileSA);
            #endif
            
            
            //std::cerr << "-->" << numEle << "\t" << GSA.text << "\t" << GSA.suff <<  "\t(" << GSA.lcp <<  "," << minLCP << ")\t" << c << "\n";
            
            numEle++;
        }  //end-for
    } //end-if
    
    
    std::cerr <<  " The total number of elements is " << numEle << "\n";
    
    fclose(f_ESA);
    
    return 0;
}

bool findLengthNseq( const string& input, const string& fileOutput)
{
    std::ifstream infile(input.c_str());
    
    if (infile.is_open() == false) {
        std::cerr << "findLengthNseq: could not open file \"" << input.c_str() << "\"!"<< std::endl;
        exit (EXIT_FAILURE);
    }
    
    char *fileLen = new char[fileOutput.length()+10];
    sprintf(fileLen, "%s.len", fileOutput.c_str());
    
    static FILE *OutFileLen;                  // output file of the end positions;
    OutFileLen = fopen(fileLen, "wb");
    if (OutFileLen==NULL) {
        std::cerr << "findLengthNseq: could not open file \"" << fileLen << "\"!"<< std::endl;
        exit (EXIT_FAILURE);
    }
    
    string bufChar;
    
    
    bool lenSeq = false;
    vector <dataTypelenSeq> lengthSeqVector;
    
    //Find max length and number of reads
    dataTypelenSeq charsNumber=0;
    dataTypelenSeq lengthRead=0;
    dataTypeNSeq nSeq = 0;
    dataTypeNChar lengthTexts = 0;       //Total length of all texts without $-symbols
    
    while (getline(infile, bufChar))  {
        
        if ( (bufChar.length() > 0) && ((bufChar[bufChar.length()-2] == '\r') || (bufChar[bufChar.length()-2] == '\n')) )
            bufChar[bufChar.length()-2] = '\0';
        else if ( (bufChar.length() > 0) && ((bufChar[bufChar.length()-1] == '\r') || (bufChar[bufChar.length()-1] == '\n')) )
            bufChar[bufChar.length()-1] = '\0';
        
        //std::cerr << "Nuova lettura: seq:" << bufChar << ". Length " << bufChar.length() << " the last char is "<< (char)bufChar[bufChar.length()-1] << std::endl;
        
        if( bufChar[0] == '>' ) {     //it is the title of fasta file
            nSeq++;
            //std::cerr << "--Length Sequence N. " << (int)(nSeq) << " is " << (int)charsNumber << std::endl;
            if (charsNumber != 0) {
                //charsNumber--;
                if (lengthRead != charsNumber) {
                    if (lengthRead < charsNumber)
                        lengthRead = charsNumber;
                    if (nSeq > 2)
                        lenSeq = true;
                }
                
                //dataTypeNSeq numchar = fwrite (&charsNumber, sizeof(dataTypelenSeq), 1 , OutFileLen);
                //assert( numchar == 1); // we should always read the same number of characters
                lengthSeqVector.push_back (charsNumber);
                lengthTexts += charsNumber;
                //if ((verboseEncode==1) || (verboseDistance ==1))
                //std::cerr << "Length Sequence N. " << (int)(nSeq-1) << " is " << (int)charsNumber << " nSeq= "<< nSeq << std::endl;
            }
            
            charsNumber = 0;
        }
        else   //no title
        {
            // increase the counter of chars
            charsNumber = charsNumber + bufChar.length();
            //std::cerr << "Stessa sequenza --Length Sequence N. " << (int)(nSeq) << " is " << (int)charsNumber << std::endl;
        }
        
    }   //end-while
    //charsNumber--;
    
    //Update the max length of the reads
    if (lengthRead != charsNumber) {
        if (lengthRead < charsNumber)
            lengthRead = charsNumber;
        lenSeq = true;
    }
    lengthSeqVector.push_back (charsNumber);
    lengthTexts += charsNumber;
    
    if (lenSeq == true)
        cerr << "The  (new and-or old) reads have a different length." << endl;
    else
        cerr << "The (new and-or old) reads have the same length." << endl;
    
    //Store vector in file
    dataTypeNSeq numchar = fwrite (&lengthSeqVector[0], sizeof(dataTypelenSeq), lengthSeqVector.size(), OutFileLen);
    
    std::cerr << "lengthTexts= " << (lengthTexts) << " lengthRead= " << (int)lengthRead << std::endl;
    std::cerr << "lengthSeqVector.size= " << (lengthSeqVector.size()) << " nSeq= " << nSeq << std::endl;
    assert (lengthSeqVector.size() == numchar);
    assert (lengthSeqVector.size() == nSeq);
    
    lengthSeqVector.clear();
    //lengthSeqVector.shrink_to_fit();
    
    delete [] fileLen;
    fclose(OutFileLen);
    infile.close();
    
    return true;
}
