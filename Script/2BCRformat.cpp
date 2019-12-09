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
    fnLCP = fileFasta+".lcp\0";
    fnDA =fileFasta+".da\0";
    fnBWT = fileFasta+ ".ebwt";
    
	#if (OUTPUT_FORMAT == 4)
		ElementTypeEBWTda	ele;						
	#endif
	#if (OUTPUT_FORMAT == 5)
		ElementTypeLCPda ele;
	#endif
	#if (OUTPUT_FORMAT == 6)
		ElementTypeGSA ele;
	#endif
	
	 uchar bwt;
    uchar c;
    dataTypelenSeq lcp;
    dataTypeNSeq idSeq;
	
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
        fseek(InDA,0, SEEK_SET);
    #endif
    
    dataTypeNChar numcharLCP, numcharDA, numcharBWT;
    dataTypeNChar numEle=0, numRead=0, numWrite=0;
    
    /*
    dataTypelenSeq *bufferLCP = new dataTypelenSeq[BUFFERLCPSIZE];
    dataTypeNSeq *bufferEle = new dataTypeNSeq[BUFFERLCPSIZE];
    dataTypedimAlpha *BWTbuffer= new dataTypedimAlpha[BUFFERLCPSIZE];
    */
    
    
    
    
    
    //BCR 
	
	//std::cerr << "file BWT : "  << fnBWT <<  "." << std::endl;
	//std::cerr << "file LCP: "  << fnLCP <<  "." << std::endl;
	//std::cerr << "file PairSA: "  << fnPairSA <<  "." << std::endl;
	//std::cerr << "file EGSA: "  << fnEGSA <<  "." << std::endl;	
	
	//Open BCR files
	#if ( (OUTPUT_FORMAT == 3) || ( (OUTPUT_FORMAT == 4) || (OUTPUT_FORMAT == 5) || (OUTPUT_FORMAT == 6) ) )
		std::cerr << "Read the entire four files for BWT/LCP/DA/SA.\n";	
		
		#if ( (OUTPUT_FORMAT == 3) || ( (OUTPUT_FORMAT == 5) || (OUTPUT_FORMAT == 6) ) )
			char *fnBWT = new char[fileOutput.size()+8];
			static FILE *OutFileBWT; 
			sprintf (fnBWT,"%s%s",fileOutput.c_str(),".ebwt");
			OutFileBWT = fopen(fnBWT, "wb");
			if (OutFileBWT==NULL) {
				std::cerr << "Error opening " << std::endl;
				exit (EXIT_FAILURE);
			}
		#endif 
		
		#if BUILD_LCP == 1
			#if ( (OUTPUT_FORMAT == 3) || (OUTPUT_FORMAT == 4) || (OUTPUT_FORMAT == 6) ) 
				char *fnLCP_out = new char[fileOutput.size()+8];
				static FILE *OutFileLCP; 
				sprintf (fnLCP_out,"%s%s",fileOutput.c_str(),".lcp");
				OutFileLCP = fopen(fnLCP_out, "wb");
				if (OutFileLCP==NULL) {
					std::cerr << "Error opening " << std::endl;
					exit (EXIT_FAILURE);
				}
			#endif
		#endif
			
		#if BUILD_SA == 1
			#if ( (OUTPUT_FORMAT == 3) || (OUTPUT_FORMAT == 4) || (OUTPUT_FORMAT == 5) )
				char *fnSA_out = new char[fileOutput.size()+8];
				static FILE *OutFileSA; 
				sprintf (fnSA_out,"%s%s",fileOutput.c_str(),".posSA");
				OutFileSA = fopen(fnSA_out, "wb");
				if (OutFileSA==NULL) {
					std::cerr << "ls Error opening " << std::endl;
					exit (EXIT_FAILURE);
				}
			#endif
		#endif
		
		#if BUILD_DA == 1
			#if (OUTPUT_FORMAT == 3) 
				char *fnDA_out = new char[fileOutput.size()+8];
				static FILE *OutFileDA; 
				sprintf (fnDA_out,"%s%s",fileOutput.c_str(),".da");
				OutFileDA = fopen(fnDA_out, "wb");
				if (OutFileDA==NULL) {
					std::cerr << "Error opening " << std::endl;
					exit (EXIT_FAILURE);
				}
			#endif
		#endif

		#if BUILD_DA == 1
			#if (OUTPUT_FORMAT == 4)
				char *fnebwtDa = new char[fileOutput.size()+8];
				static FILE *OutFileebwtDa; 
				sprintf (fnebwtDa,"%s%s",fileOutput.c_str(),".ebwtDa");
				OutFileebwtDa = fopen(fnebwtDa, "wb");
				if (OutFileebwtDa==NULL) {
					std::cerr << "Error opening " << std::endl;
					exit (EXIT_FAILURE);
				}
			#endif
		#endif
		
		#if (BUILD_DA == 1) &&  (BUILD_LCP == 1)
			#if (OUTPUT_FORMAT == 5)
				char *fnlcpDa = new char[fileOutput.size()+8];
				static FILE *OutFileLCPDa; 
				sprintf (fnlcpDa,"%s%s",fileOutput.c_str(),".lcpDa");
				OutFileLCPDa = fopen(fnlcpDa, "wb");
				if (OutFileLCPDa==NULL) {
					std::cerr << "Error opening " << std::endl;
					exit (EXIT_FAILURE);
				}
			#endif
		#endif
		
		#if (BUILD_DA == 1) && (BUILD_SA == 1)
			#if (OUTPUT_FORMAT == 6)
				char *fnGSA = new char[fileOutput.size()+8];
				static FILE *OutFileGSA; 
				sprintf (fnGSA,"%s%s",fileOutput.c_str(),".gsa");
				OutFileGSA = fopen(fnGSA, "wb");
				if (OutFileGSA==NULL) {
					std::cerr << "Error opening " << std::endl;
					exit (EXIT_FAILURE);
				}
			#endif		
		#endif
	#endif
	
	dataTypeNChar val;
	dataTypeNSeq numSeq=0;

	while ( fread(&bwt,sizeof(dataTypedimAlpha),1, InBWT) > 0 )   {
        //std::cerr << "bwt\tlcp\tpos\tnumSeq\tSA\tQS\n";
		#if BUILD_LCP == 1
			numcharLCP = fread(&lcp,sizeof(dataTypelenSeq),1,InLCP);
		#endif
		#if BUILD_DA == 1
			numcharDA = fread(&idSeq,sizeof(dataTypeNSeq),1,InDA);
		#endif
		#if BUILD_SA == 1
			numcharPairSA = fread(&val,sizeof(dataTypeNChar),1,InSA);
		#endif
        //std::cerr << "-->" << numEle << "\t" << bwt << "\t" << idSeq<<  "\t" << GSA.lcp <<  "," << minLCP << ")\t" << c << "\n";
        if ((unsigned int)bwt == 0) {
            c=TERMINATE_CHAR;
			numSeq++;
		}
        else
            c=bwt;
		
			#if ( (OUTPUT_FORMAT == 3) || (OUTPUT_FORMAT == 5) || (OUTPUT_FORMAT == 6) )
					fwrite(&c, sizeof(uchar), 1, OutFileBWT);
			#endif

			#if  (BUILD_LCP == 1)
				#if ( (OUTPUT_FORMAT == 3) || (OUTPUT_FORMAT == 4) || (OUTPUT_FORMAT == 6) )
					fwrite (&lcp, sizeof(dataTypelenSeq), 1 , OutFileLCP);
					//assert(numcharLCP == numcharWrite); 
				#endif
			#endif

			#if (BUILD_SA == 1)
				#if ( (OUTPUT_FORMAT == 3) || (OUTPUT_FORMAT == 4) || (OUTPUT_FORMAT == 5) ) 
					fwrite (&val, sizeof(dataTypelenSeq), 1 , OutFileSA);
					//assert(numcharLCP == numcharWrite); 
				#endif
			#endif
			
			#if (BUILD_DA == 1) 
				#if (OUTPUT_FORMAT == 3) 
					fwrite (&idSeq, sizeof(dataTypeNSeq), 1 , OutFileDA);
					//assert(numcharLCP == numcharWrite); 
				#endif
			#endif
			
			#if ( (OUTPUT_FORMAT == 4) || (OUTPUT_FORMAT == 5) || (OUTPUT_FORMAT == 6) )
				
					#if (OUTPUT_FORMAT == 4) && (BUILD_DA == 1) 
						ele.bwt=c;							
						ele.da = idSeq;
						fwrite(&ele, sizeof(ElementTypeEBWTda), 1, OutFileebwtDa);
					#endif
					#if (OUTPUT_FORMAT == 5) && (BUILD_LCP == 1)  && (BUILD_DA == 1) 
						ele.lcp = GSA.lcp;
						ele.da = idSeq;
						fwrite(&ele, sizeof(ElementTypeLCPda), 1, OutFileLCPDa);
					#endif
					#if (OUTPUT_FORMAT == 6) && (BUILD_SA == 1)  && (BUILD_DA == 1) 
						ele.sa = val;
						ele.da = idSeq;
						fwrite(&ele, sizeof(ElementTypeGSA), 1, OutFileGSA);
					#endif
				
			#endif
				
			
			//std::cerr << "-->" << numEle << "\t" << GSA.text << "\t" << GSA.suff <<  "\t(" << GSA.lcp <<  "," << minLCP << ")\t" << c << "\n";
				
			numEle++;
		
	}

    
    std::cerr <<  "The total number of elements is " << numEle << "\n";
    
    fclose(InBWT);
	#if ( (OUTPUT_FORMAT == 3) || (OUTPUT_FORMAT == 5) || (OUTPUT_FORMAT == 6) )
		fclose(OutFileBWT);
	#endif
    
    #if BUILD_LCP == 1
        fclose(InLCP);
		#if ( (OUTPUT_FORMAT == 3) || (OUTPUT_FORMAT == 4) || (OUTPUT_FORMAT == 6) )
			fclose(OutFileLCP);
			delete [] fnLCP_out;
		#endif        
    #endif
	
    #if BUILD_DA == 1
		#if (OUTPUT_FORMAT == 3)
			fclose(OutFileDA);
			delete [] fnDA_out;
		#endif
        fclose(InDA);
		#if (OUTPUT_FORMAT == 4)
			fclose(OutFileebwtDa);
			delete [] fnebwtDa;
		#endif 
    #endif
	
	#if (OUTPUT_FORMAT == 5) && (BUILD_LCP == 1)  && (BUILD_DA == 1) 
		fclose(OutFileLCPDa);
		delete []  fnlcpDa;
	#endif
	
	#if (OUTPUT_FORMAT == 6) && (BUILD_SA == 1)  && (BUILD_DA == 1) 
		fclose(OutFileGSA);
		delete []  fnGSA;
	#endif
    
    #if BUILD_SA == 1
        fclose(OutFileSA);
		#if ( (OUTPUT_FORMAT == 3) || (OUTPUT_FORMAT == 4) || (OUTPUT_FORMAT == 5) ) 
			fclose(OutFileSA);
			delete [] fnSA_out;
		#endif
        fclose(InSA);
    #endif
    
    
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
