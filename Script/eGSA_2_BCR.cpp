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


//See Felipe Louza tool: eGSA
	typedef unsigned int int_text;
	typedef unsigned int int_suff; 	//-2^31 to 2^31
	typedef unsigned int int_lcp;
	typedef unsigned char int8; //0 to 2^8

	#define PREFIX_SIZE 		10
	
	typedef struct{
		dataTypeNChar	text;
		dataTypelenSeq	suff;
		dataTypelenSeq 	lcp;
		
		uchar		bwt;	
	} t_GSA;
	

int main(int argc, char **argv) {
	 
	#if (OUTPUT_FORMAT == 4)
		ElementTypeEBWTda	ele;						
	#endif
	#if (OUTPUT_FORMAT == 5)
		ElementTypeLCPda ele;
	#endif
	#if (OUTPUT_FORMAT == 6)
		ElementTypeGSA ele;
	#endif
	 
	if( argc != 2 )
    {
      std::cerr << "usage: " << argv[0] << " input" << std::endl;
	  std::cerr << "where input is the filename without 0.gesa " << std::endl;
       exit(1);
    }
	
	string fileInput = argv[1];
	std::cerr << "ReadEGSA: File: " << fileInput << std::endl;
	
	dataTypeNChar nShortRead, nTarget;
	

	//Open EGSA file	
	FILE 		*f_ESA;			// pointer to the ESA input file
	string fnEGSA;
	fnEGSA = fileInput + "." + "0" + ".gesa\0";
	
	f_ESA = fopen(fnEGSA.c_str(), "rb");
	if (f_ESA == NULL) {
		std::cerr << "readEGSA: Error opening: " << fnEGSA << std::endl;
		exit (EXIT_FAILURE);
	}	
	fseek(f_ESA , 0, SEEK_SET);			
	t_GSA GSA;
	
	
	//BCR 
	
	//std::cerr << "file BWT : "  << fnBWT <<  "." << std::endl;
	//std::cerr << "file LCP: "  << fnLCP <<  "." << std::endl;
	//std::cerr << "file PairSA: "  << fnPairSA <<  "." << std::endl;
	//std::cerr << "file EGSA: "  << fnEGSA <<  "." << std::endl;	
	
	//Open BCR files
	#if ( (OUTPUT_FORMAT == 3) || ( (OUTPUT_FORMAT == 4) || (OUTPUT_FORMAT == 5) || (OUTPUT_FORMAT == 6) ) )
		std::cerr << "Build the entire four files for BWT/LCP/DA/SA.\n";	
		
		#if ( (OUTPUT_FORMAT == 3) || ( (OUTPUT_FORMAT == 5) || (OUTPUT_FORMAT == 6) ) )
			char *fnBWT = new char[fileInput.size()+8];
			static FILE *OutFileBWT; 
			sprintf (fnBWT,"%s%s",fileInput.c_str(),".ebwt");
			OutFileBWT = fopen(fnBWT, "wb");
			if (OutFileBWT==NULL) {
				std::cerr << "storeEGSAcomplete: Error opening " << std::endl;
				exit (EXIT_FAILURE);
			}
		#endif 
		
		#if ( (OUTPUT_FORMAT == 3) || (OUTPUT_FORMAT == 4) || (OUTPUT_FORMAT == 6) ) 
			char *fnLCP = new char[fileInput.size()+8];
			static FILE *OutFileLCP; 
			sprintf (fnLCP,"%s%s",fileInput.c_str(),".lcp");
			OutFileLCP = fopen(fnLCP, "wb");
			if (OutFileLCP==NULL) {
				std::cerr << "storeEGSAcomplete: Error opening " << std::endl;
				exit (EXIT_FAILURE);
			}
		#endif
		#if ( (OUTPUT_FORMAT == 3) || (OUTPUT_FORMAT == 4) || (OUTPUT_FORMAT == 5) )
			char *fnSA = new char[fileInput.size()+8];
			static FILE *OutFileSA; 
			sprintf (fnSA,"%s%s",fileInput.c_str(),".posSA");
			OutFileSA = fopen(fnSA, "wb");
			if (OutFileSA==NULL) {
				std::cerr << "storeEGSAcomplete: Error opening " << std::endl;
				exit (EXIT_FAILURE);
			}
		#endif
		#if (OUTPUT_FORMAT == 3) 
			char *fnDA = new char[fileInput.size()+8];
			static FILE *OutFileDA; 
			sprintf (fnDA,"%s%s",fileInput.c_str(),".da");
			OutFileDA = fopen(fnDA, "wb");
			if (OutFileDA==NULL) {
				std::cerr << "storeEGSAcomplete: Error opening " << std::endl;
				exit (EXIT_FAILURE);
			}
		#endif
		
		
		#if (OUTPUT_FORMAT == 4)
			char *fnebwtDa = new char[fileInput.size()+8];
			static FILE *OutFileebwtDa; 
			sprintf (fnebwtDa,"%s%s",fileInput.c_str(),".ebwtDa");
			OutFileebwtDa = fopen(fnebwtDa, "wb");
			if (OutFileebwtDa==NULL) {
				std::cerr << "storeEGSAcomplete: Error opening " << std::endl;
				exit (EXIT_FAILURE);
			}
		#endif
		#if (OUTPUT_FORMAT == 5)
			char *fnlcpDa = new char[fileInput.size()+8];
			static FILE *OutFileLCPDa; 
			sprintf (fnlcpDa,"%s%s",fileInput.c_str(),".lcpDa");
			OutFileLCPDa = fopen(fnlcpDa, "wb");
			if (OutFileLCPDa==NULL) {
				std::cerr << "storeEGSAcomplete: Error opening " << std::endl;
				exit (EXIT_FAILURE);
			}
		#endif
		#if (OUTPUT_FORMAT == 6)
			char *fnGSA = new char[fileInput.size()+8];
			static FILE *OutFileGSA; 
			sprintf (fnGSA,"%s%s",fileInput.c_str(),".gsa");
			OutFileGSA = fopen(fnGSA, "wb");
			if (OutFileGSA==NULL) {
				std::cerr << "storeEGSAcomplete: Error opening " << std::endl;
				exit (EXIT_FAILURE);
			}
		#endif		
	#endif
	
	
	dataTypeNChar numcharBWT, numcharPairSA, numcharLCP, numEle=0, numRead=0, numWrite=0;

	uchar c;
	
	dataTypelenSeq minLCP=0;
	//std::cout << "nSeq" << "\t" << "Pos" << "\t" <<  "LCP" << "\t" <<  "minLCP" << "\t" <<  "BWT" << "\t" << "\n";
	while ( !feof(f_ESA) )   {
			//std::cerr << "bwt\tlcp\tpos\tnumSeq\tSA\tQS\n";
							
			numRead=fread(&GSA.text, sizeof(int_text), 1, f_ESA);	
			if (numRead > 0) {
				fread(&GSA.suff, sizeof(int_suff), 1, f_ESA);	
				fread(&GSA.lcp, sizeof(int_lcp), 1, f_ESA);	
				fread(&GSA.bwt, sizeof(int8), 1, f_ESA);	
				
				if (GSA.bwt == '\0')
					c='$';
				else
					c=GSA.bwt;
				
				#if ( (OUTPUT_FORMAT == 3) || (OUTPUT_FORMAT == 5) || (OUTPUT_FORMAT == 6) )
						fwrite(&c, sizeof(uchar), 1, OutFileBWT);
				#endif

				
				#if ( (OUTPUT_FORMAT == 3) || (OUTPUT_FORMAT == 4) || (OUTPUT_FORMAT == 6) )
					fwrite (&GSA.lcp, sizeof(dataTypelenSeq), 1 , OutFileLCP);
					//assert(numcharLCP == numcharWrite); 
				#endif
				#if ( (OUTPUT_FORMAT == 3) || (OUTPUT_FORMAT == 4) || (OUTPUT_FORMAT == 5) ) 
					fwrite (&GSA.suff, sizeof(dataTypelenSeq), 1 , OutFileSA);
					//assert(numcharLCP == numcharWrite); 
				#endif
				#if (OUTPUT_FORMAT == 3) 
					fwrite (&GSA.text, sizeof(dataTypeNSeq), 1 , OutFileDA);
					//assert(numcharLCP == numcharWrite); 
				#endif
				
				#if ( (OUTPUT_FORMAT == 4) || (OUTPUT_FORMAT == 5) || (OUTPUT_FORMAT == 6) )
					
						#if (OUTPUT_FORMAT == 4)
							ele.bwt=c;							
							ele.da = GSA.text;
							fwrite(&ele, sizeof(ElementTypeEBWTda), 1, OutFileebwtDa);
						#endif
						#if (OUTPUT_FORMAT == 5)
							ele.lcp = GSA.lcp;
							ele.da = GSA.text;
							fwrite(&ele, sizeof(ElementTypeLCPda), 1, OutFileLCPDa);
						#endif
						#if (OUTPUT_FORMAT == 6)
							ele.sa = GSA.suff;
							ele.da = GSA.text;
							fwrite(&ele, sizeof(ElementTypeGSA), 1, OutFileGSA);
						#endif
					
				#endif
					
				
				//std::cerr << "-->" << numEle << "\t" << GSA.text << "\t" << GSA.suff <<  "\t(" << GSA.lcp <<  "," << minLCP << ")\t" << c << "\n";
					
				numEle++;
			}

	}
	
	

	std::cerr <<  " The total number of elements is " << numEle << "\n";
	
	fclose(f_ESA);
	
	return 0;
}

