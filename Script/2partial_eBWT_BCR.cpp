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
#include "../Parameters.h"
#include "../Tools.h"

using namespace std;

#define BUFFERLCPSIZE 1024


bool findLengthNseq( const string& input, const string& fileOutput, dataTypeNSeq, dataTypeNChar);

int main(int argc, char **argv) {
    
    if( argc != 3 )
    {
        std::cerr << "usage: " << argv[0] << " input" << std::endl;
        std::cerr << "where:" << std::endl;
        std::cerr << "  input is the filename without extension " << std::endl;
        std::cerr << "  output is the filename for BCR partial files " << std::endl;
        exit(1);
    }
    
    string fileFasta = argv[1];
    string fileOutput = argv[2];
    
	dataTypeNSeq numTotSeq=0;
	dataTypeNChar lenghtBWT=0;
    
    string fnLCP, fnDA, fnBWT;
    fnLCP = fileFasta+".4.lcp\0";
    fnDA =fileFasta+".da\0";
    fnBWT = fileFasta+ ".bwt";
    
    //Open LCP and DA and BWT files
    FILE *InBWT = fopen(fnBWT.c_str(), "rb");
    if (InBWT==NULL) {
        std::cerr << "Error opening " << fnBWT << "." << std::endl;
        //printf("fopen failed, errno = %d\n", errno);
        exit (EXIT_FAILURE);
    }
    fseek(InBWT, 0, SEEK_SET);
    #if BUILD_LCP == 1
        FILE *InLCP = fopen(fnLCP.c_str(), "rb");
        if (InLCP==NULL) {
            std::cerr << "Error opening " << fnLCP << "." << std::endl;
            //printf("fopen failed, errno = %d\n", errno);
            exit (EXIT_FAILURE);
        }
        fseek(InLCP, 0, SEEK_SET);
    #endif
    #if BUILD_DA == 1
        FILE *InDA = fopen(fnDA.c_str(), "rb");
        if (InDA==NULL) {
            std::cerr << "Error opening " << fnDA << "." << std::endl;
            //printf("fopen failed, errno = %d\n", errno);
            exit (EXIT_FAILURE);
        }
        fseek(InDA,sizeof(dataTypeNSeq), SEEK_SET);
    #endif
    
    dataTypeNChar numcharLCP, numcharDA, numcharBWT;
    dataTypeNChar numEle=0, numRead=0, numWrite=0;
    
    /*
    dataTypelenSeq *bufferLCP = new dataTypelenSeq[BUFFERLCPSIZE];
    dataTypeNSeq *bufferEle = new dataTypeNSeq[BUFFERLCPSIZE];
    dataTypedimAlpha *BWTbuffer= new dataTypedimAlpha[BUFFERLCPSIZE];
    */
    
    dataTypeNChar freq[256];  //contains the distribution of the symbols.
    
    for (dataTypedimAlpha z = 0 ; z < SIZE_ALPHA-1; z++)
        freq[z]=0;
    freq[SIZE_ALPHA-1]=0;
    
    
    //First reading in order to find the alphabet
    std::cerr << "find the alphabet" << std::endl;
    uchar bwt;
    uchar c;
    dataTypelenSeq lcp;
    dataTypeNSeq idSeq;
    numcharBWT = fread(&bwt,sizeof(dataTypedimAlpha),1, InBWT);
    
	#if BUILD_LCP == 1
		numcharLCP = fread(&lcp,sizeof(dataTypelenSeq),1,InLCP);
		assert(numcharLCP==numcharBWT);
		dataTypelenSeq minLCP=lcp, maxLCP=lcp;
	#endif
	#if BUILD_DA == 1
		numcharDA = fread(&idSeq,sizeof(dataTypeNSeq),1,InDA);
		assert(numcharBWT==numcharDA);
    #endif
    if ((unsigned int)bwt == 0)
        c=TERMINATE_CHAR;
    else
        c=bwt;
    freq[(unsigned int)(c)]=1;
	numEle++;
    
	dataTypeNSeq numSeq=0;
	
    while ( fread(&bwt,sizeof(dataTypedimAlpha),1, InBWT) > 0 )   {
        //std::cerr << "bwt\tlcp\tpos\tnumSeq\tSA\tQS\n";
		#if BUILD_LCP == 1
			numcharLCP = fread(&lcp,sizeof(dataTypelenSeq),1,InLCP);
		#endif
		#if BUILD_DA == 1
			numcharDA = fread(&idSeq,sizeof(dataTypeNSeq),1,InDA);
		#endif
        //std::cerr << "-->" << numEle << "\t" << GSA.text << "\t" << GSA.suff <<  "\t(" << GSA.lcp <<  "," << minLCP << ")\t" << c << "\n";
        if ((unsigned int)bwt == 0) {
            c=TERMINATE_CHAR;
			numSeq++;
		}
        else
            c=bwt;

        freq[(unsigned int)(c)]++;
        #if BUILD_LCP == 1
			if (minLCP > lcp)
					minLCP=lcp;
			if (maxLCP < lcp)
					maxLCP = lcp;
		#endif
            //std::cerr << "-->" << numEle << "\t" << GSA.text << "\t" << GSA.suff <<  "\t(" << GSA.lcp <<  "," << minLCP << ")\t" << c << "\n";
            
        numEle++;
    }//end-while
    
	
    //Check
    assert(feof(InBWT));
    #if BUILD_LCP == 1
        int num = fread(&lcp,sizeof(dataTypelenSeq),1,InLCP);
		if (!feof(InLCP)) {
			std::cerr << "lcp: " << lcp << " num " << num << "\n";
		}
        assert(feof(InLCP));
    #endif
    #if BUILD_DA == 1
        int num1 = fread(&idSeq,sizeof(dataTypeNSeq),1,InDA);
		if (!feof(InLCP)) {
			std::cerr << "lcp: " << lcp << " num " << num1 << "\n";
		}
        assert(feof(InDA));
    #endif
    
    dataTypedimAlpha *alpha = new dataTypedimAlpha[SIZE_ALPHA];;  //Corresponding between alpha[i] and the symbol as char
    dataTypedimAlpha sizeAlpha=0;
    for (dataTypedimAlpha i = 0; i < SIZE_ALPHA-1; i++)
        if (freq[i] > 0) {
			alpha[i] = sizeAlpha;
            sizeAlpha++;			
        }
    if (freq[SIZE_ALPHA-1] > 0) {
		alpha[SIZE_ALPHA-1] = sizeAlpha;
        sizeAlpha++;
    }
    
	std::cerr << "Compute tableocc with where the size of the alphabet is " << (int)sizeAlpha << "." << std::endl;
	dataTypeNChar** tableOcc; //contains the number of occurrences of each symbol
	tableOcc = new dataTypeNChar*[sizeAlpha];
	//Counting for each pile, es. $-pile, A-pile, C-pile, G-pile, N-pile, T-pile
	for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++) {
			tableOcc[j] = new dataTypeNChar[sizeAlpha];
	  }
	for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++)
		for (dataTypedimAlpha h = 0 ; h < sizeAlpha; h++)
			tableOcc[j][h]=0;
	
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
    
    std::cerr << "There are " <<numEle << " symbols in the input file (ASCII, char, freq, code):\n";
    for (dataTypedimAlpha i = 0; i < SIZE_ALPHA-1; i++)
        if (freq[i] > 0)
            std::cerr << (unsigned int)i << " " << i << " " << freq[i] << "\n";
    if (freq[SIZE_ALPHA-1] > 0)
        std::cerr << (unsigned int)SIZE_ALPHA-1 << " " << SIZE_ALPHA-1 << " " << freq[SIZE_ALPHA-1] << " " << "\n";
    
    
    
    
    
    //Open BCR partial files
    
    
    std::cerr << "Build the partial BCR files for eBWT/LCP/DA/SA.\n";
    char *fnOutBWT = new char[fileOutput.size()+100];
    static FILE *OutFileBWT;
    fseek(InBWT, 0, SEEK_SET);
    #if BUILD_LCP == 1
        char *fnOutLCP = new char[fileOutput.size()+100];
        static FILE *OutFileLCP;
        fseek(InLCP, 0, SEEK_SET);
    #endif
    #if BUILD_DA == 1
        char *fnOutDA = new char[fileOutput.size()+100];
        static FILE *OutFileDA;
        fseek(InDA,0, SEEK_SET);
    #endif
    /*
    #if BUILD_SA == 1
        char *fnOutSA = new char[fileOutput.size()+100];
        static FILE *OutFileSA;
    #endif
     */
    numEle=0;
    numRead=0;
    numWrite=0;
	#if BUILD_LCP == 1
		minLCP=0;
	#endif
    //std::cout << "nSeq" << "\t" << "Pos" << "\t" <<  "LCP" << "\t" <<  "minLCP" << "\t" <<  "BWT" << "\t" << "\n";
    
    for (dataTypedimAlpha currentPile = 0 ; currentPile < SIZE_ALPHA-1; currentPile++) {
        
		if (freq[currentPile] > 0) {
            sprintf (fnOutBWT,"%s%s%d",fileOutput.c_str(),".ebwt_",currentPile);
            OutFileBWT = fopen(fnOutBWT, "wb");
            if (OutFileBWT==NULL) {
                std::cerr << "Error opening: " << fnBWT << std::endl;
                exit (EXIT_FAILURE);
            }
            #if BUILD_LCP == 1
                sprintf (fnOutLCP,"%s%s%d",fileOutput.c_str(),".lcp_",currentPile);
                OutFileLCP = fopen(fnOutLCP, "wb");
                if (OutFileLCP==NULL) {
                    std::cerr << "Error opening: " << fnOutLCP << std::endl;
                    exit (EXIT_FAILURE);
                }
            #endif
            #if BUILD_DA == 1
                sprintf (fnOutDA,"%s%s%d",fileOutput.c_str(),".da_",currentPile);
                OutFileDA = fopen(fnOutDA, "wb");
                if (OutFileDA==NULL) {
                    std::cerr << "Error opening: " << fnOutDA << std::endl;
                    exit (EXIT_FAILURE);
                }
            #endif
            /*
            #if BUILD_SA == 1
                sprintf (fnOutSA,"%s%s%d",fileOutput.c_str(),".sa_",currentPile);
                OutFileSA = fopen(fnOutSA, "wb");
                if (OutFileSA==NULL) {
                    std::cerr << "Error opening: " << fnOutSA << std::endl;
                    exit (EXIT_FAILURE);
                }
            #endif
             */
            //std::cerr << "bwt\tlcp\tpos\tnumSeq\tSA\tQS\n";
            
			for (dataTypeNChar j = 0 ; j < freq[currentPile]; j++) {
                numcharBWT = fread(&bwt,sizeof(dataTypedimAlpha),1, InBWT);
                if ((unsigned int)bwt == 0) {
                    c=TERMINATE_CHAR;
					numTotSeq++;	
				}
                else
                    c=bwt;
				lenghtBWT++;
                numWrite = fwrite(&c, sizeof(uchar), 1, OutFileBWT);				
				assert(numcharBWT == numWrite);
								
				tableOcc[alpha[(unsigned int)currentPile]][alpha[(unsigned int)c]]++;      
				//counting the number of occurrences in BWT of the $-pile

				
				
                #if BUILD_LCP == 1
                    numcharLCP = fread(&lcp,sizeof(dataTypelenSeq),1,InLCP);
                    assert(numcharLCP==numcharBWT);
                    numWrite = fwrite (&lcp, sizeof(dataTypelenSeq), 1 , OutFileLCP);
					assert(numcharLCP == numWrite);
                #endif
                #if BUILD_DA == 1
                    numcharDA = fread(&idSeq,sizeof(dataTypeNSeq),1,InDA);
					
					if (numcharDA!=numcharBWT)
						std::cerr <<  "numcharDA "<< numcharDA << " and numcharBWT" << numcharBWT << "\n";	
					
                    assert(numcharDA==numcharBWT);
                    numWrite = fwrite (&idSeq, sizeof(dataTypeNSeq), 1 , OutFileDA);
					assert(numcharDA == numWrite);
                #endif
                
                /*
                #if BUILD_SA == 1
                    fwrite (&suff, sizeof(dataTypelenSeq), 1 , OutFileSA);
                #endif
                 */
                
                //std::cerr << "-->" << numEle << "\t" << GSA.text << "\t" << GSA.suff <<  "\t(" << GSA.lcp <<  "," << minLCP << ")\t" << c << "\n";
                
                numEle++;
            }  //end-for
        } //end-if
		//else

    }  //end-for
    
    
    
    //std::cout << "currentPile " << (int)SIZE_ALPHA-1 << "\n";
    if (freq[SIZE_ALPHA-1] > 0) {
        sprintf (fnOutBWT,"%s%s%d",fileOutput.c_str(),".ebwt_",SIZE_ALPHA-1);
        OutFileBWT = fopen(fnOutBWT, "wb");
        if (OutFileBWT==NULL) {
            std::cerr << "Error opening: " << fnOutBWT << std::endl;
            exit (EXIT_FAILURE);
        }
        #if BUILD_LCP == 1
            sprintf (fnOutLCP,"%s%s%d",fileOutput.c_str(),".lcp_",SIZE_ALPHA-1);
            OutFileLCP = fopen(fnOutLCP, "wb");
            if (OutFileLCP==NULL) {
                std::cerr << "Error opening: " << fnOutLCP << std::endl;
                exit (EXIT_FAILURE);
            }
        #endif
        #if BUILD_DA == 1
            sprintf (fnOutDA,"%s%s%d",fileOutput.c_str(),".da_",SIZE_ALPHA-1);
            OutFileDA = fopen(fnOutDA, "wb");
            if (OutFileDA==NULL) {
                std::cerr << "Error opening: " << fnOutDA << std::endl;
                exit (EXIT_FAILURE);
            }
        #endif
        /*
        #if BUILD_SA == 1
            sprintf (fnOutSA,"%s%s%d",fileOutput.c_str(),".sa_",SIZE_ALPHA-1);
            OutFileSA = fopen(fnOutSA, "wb");
            if (OutFileSA==NULL) {
                std::cerr << "Error opening: " << fnOutSA << std::endl;
                exit (EXIT_FAILURE);
            }
        #endif
         */
        //std::cerr << "bwt\tlcp\tpos\tnumSeq\tSA\tQS\n";
        
        for (dataTypeNChar j = 0 ; j < freq[SIZE_ALPHA-1]; j++) {
            
            numcharBWT = fread(&bwt,sizeof(dataTypedimAlpha),1, InBWT);
            if ((unsigned int)bwt == 0)
                c=TERMINATE_CHAR;
            else
                c=bwt;
            fwrite(&c, sizeof(uchar), 1, OutFileBWT);
            #if BUILD_LCP == 1
                numcharLCP = fread(&lcp,sizeof(dataTypelenSeq),1,InLCP);
                assert(numcharLCP==numcharBWT);
                fwrite (&lcp, sizeof(dataTypelenSeq), 1 , OutFileLCP);
            #endif
            #if BUILD_DA == 1
                numcharDA = fread(&idSeq,sizeof(dataTypeNSeq),1,InDA);
                assert(numcharDA==numcharBWT);
                fwrite (&idSeq, sizeof(dataTypeNSeq), 1 , OutFileDA);
            #endif
            /*
            #if BUILD_SA == 1
                fwrite (&suff, sizeof(dataTypelenSeq), 1 , OutFileSA);
            #endif
            */
            
            //std::cerr << "-->" << numEle << "\t" << GSA.text << "\t" << GSA.suff <<  "\t(" << GSA.lcp <<  "," << minLCP << ")\t" << c << "\n";
            
            numEle++;
        }  //end-for
    } //end-if
    
    
    std::cerr <<  "The total number of elements is " << numEle << "\n";
    
    fclose(InBWT);
    fclose(OutFileBWT);
    #if BUILD_LCP == 1
        fclose(InLCP);
        fclose(OutFileLCP);
    #endif
    #if BUILD_DA == 1
        fclose(OutFileDA);
        fclose(InDA);
    #endif
    /*
    #if BUILD_SA == 1
        fclose(OutFileSA);
        fclose(InSA);
    #endif
    */
    
    std::cout << "Store info files\n";
	char *fnAux = new char[fileOutput.size()+100];
	sprintf (fnAux,"%s%s",fileOutput.c_str(),".info\0");
	FILE* OutFile = fopen(fnAux, "wb");
	if (OutFile==NULL) {
		std::cerr << "(lengthBWT+NSequences+sizeAlpha) Error opening " << fnAux << std::endl;
		exit (EXIT_FAILURE);
	}
	
	fwrite(&lenghtBWT,sizeof(dataTypeNChar),1,OutFile);
	fwrite(&numTotSeq,sizeof(dataTypeNSeq),1,OutFile);
	fwrite(&sizeAlpha,sizeof(dataTypedimAlpha),1,OutFile);
	fclose(OutFile);
	delete [] fnAux;
	
	
	std::cout << "Store tableOcc files\n";
	char *fnAuxTable = new char[fileOutput.size()+100];
	sprintf (fnAuxTable,"%s%s",fileOutput.c_str(),".table\0");
	FILE* OutFileTable = fopen(fnAuxTable, "wb");
	if (OutFileTable==NULL) {
		std::cerr << "(tableOcc) Error opening " << fnAuxTable << std::endl;
		exit (EXIT_FAILURE);
	}
	std::cerr << "Writing tableOcc:\n";
	for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++) {
		std::cerr << (int)j << ": ";
		for (dataTypedimAlpha h = 0 ; h < sizeAlpha; h++) {
			dataTypeNChar num = tableOcc[j][h];
			fwrite(&num,sizeof(dataTypeNChar),1,OutFileTable);
			std::cerr <<   tableOcc[j][h] << " \t";
		}
		std::cerr << std::endl;
	}
	
	fclose(OutFileTable);
	
	
	for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++) {
		delete [] tableOcc[j];
		tableOcc[j] = NULL;
	}
	delete [] tableOcc;
	
	delete [] fnAuxTable;
    
    return 0;
}

bool findLengthNseq( const string& input, const string& fileOutput, dataTypeNSeq *nSeq, dataTypeNChar *lengthTexts)
{
	
	//dataTypeNSeq nSeq = 0;
    //dataTypeNChar lengthTexts = 0;       //Total length of all texts without $-symbols
	
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
    
    
    while (getline(infile, bufChar))  {
        
        if ( (bufChar.length() > 0) && ((bufChar[bufChar.length()-2] == '\r') || (bufChar[bufChar.length()-2] == '\n')) )
            bufChar[bufChar.length()-2] = '\0';
        else if ( (bufChar.length() > 0) && ((bufChar[bufChar.length()-1] == '\r') || (bufChar[bufChar.length()-1] == '\n')) )
            bufChar[bufChar.length()-1] = '\0';
        
        //std::cerr << "Nuova lettura: seq:" << bufChar << ". Length " << bufChar.length() << " the last char is "<< (char)bufChar[bufChar.length()-1] << std::endl;
        
        if( bufChar[0] == '>' ) {     //it is the title of fasta file
            (*nSeq)++;
            //std::cerr << "--Length Sequence N. " << (int)(nSeq) << " is " << (int)charsNumber << std::endl;
            if (charsNumber != 0) {
                //charsNumber--;
                if (lengthRead != charsNumber) {
                    if (lengthRead < charsNumber)
                        lengthRead = charsNumber;
                    if (*nSeq > 2)
                        lenSeq = true;
                }
                
                //dataTypeNSeq numchar = fwrite (&charsNumber, sizeof(dataTypelenSeq), 1 , OutFileLen);
                //assert( numchar == 1); // we should always read the same number of characters
                lengthSeqVector.push_back (charsNumber);
                *lengthTexts += charsNumber;
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
    *lengthTexts += charsNumber;
    
    if (lenSeq == true)
        cerr << "The  (new and-or old) reads have a different length." << endl;
    else
        cerr << "The (new and-or old) reads have the same length." << endl;
    
    //Store vector in file
    dataTypeNSeq numchar = fwrite (&lengthSeqVector[0], sizeof(dataTypelenSeq), lengthSeqVector.size(), OutFileLen);
    
    std::cerr << "lengthTexts= " << (*lengthTexts) << " lengthRead= " << (int)lengthRead << std::endl;
    std::cerr << "lengthSeqVector.size= " << (lengthSeqVector.size()) << " nSeq= " << *nSeq << std::endl;
    assert (lengthSeqVector.size() == numchar);
    assert (lengthSeqVector.size() == *nSeq);
    
    lengthSeqVector.clear();
    //lengthSeqVector.shrink_to_fit();
    
    delete [] fileLen;
    fclose(OutFileLen);
    infile.close();
    
    return true;
}
