/**
 ** BCR is part of:
 ** BEETL: Burrows-Wheeler Extended Tool Library
 ** Documentation in: doc/BEETL.md
 **
 ** Copyright (c) 2011-2014 Illumina, Inc. **
 ** BEETL software package is
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** Citations: 
 ** 
 ** Markus J. Bauer, Anthony J. Cox and Giovanna Rosone
 ** Lightweight BWT Construction for Very Large String Collections.
 ** Proceedings of CPM 2011, pp.219-231
 ** 
 ** Markus J. Bauer, Anthony J. Cox, Giovanna Rosone and Marinella Sciortino
 ** Lightweight LCP Construction for Next-Generation Sequencing Datasets. 
 ** Proceedings of WABI 2012, pp 326-337, 2012
 
 ** Markus J. Bauer, Anthony J. Cox, Giovanna Rosone 
 ** Lightweight algorithms for constructing and inverting the BWT of string collections. 
 ** Theoretical Computer Science 483: 134-148 (2013)
 **  
 ** Anthony J. Cox, Fabio Garofalo, Giovanna Rosone, Marinella Sciortino
 ** Lightweight LCP construction for very large collections of strings. 
 ** Journal of Discrete Algorithms (2016) 37: 17-33
 **
 ** By Giovanna Rosone
 **
 **/

// Test driver for bwt collection
#include <iostream>
#include <assert.h>
#include <string.h>     // std::string, std::to_string
#include <sstream>
#include "Timer.hh"
#include <stdio.h>
#include <math.h> 
#include <fstream>

using std::cout;
using std::endl;

#include "BWTCollection.h"
//#include "Timer.hh"

using SXSI::BWTCollection;

int main(int argc, char *argv[])
{

    string inputPrevBCR="";
    #if BUILD_BCR_FROM_BCRpartials==0
        if( argc != 3 )  {
            std::cerr << "usage: " << argv[0] << " input output" << std::endl;
            exit(1);
        }
    #elif BUILD_BCR_FROM_BCRpartials==1
        if( argc != 4 )  {
            std::cerr << "usage: " << argv[0] << " input output inputPrevBCR" << std::endl;
            std::cerr << "where: " << std::endl;
            std::cerr << "\tinputPrevBCR --> prefix filename BCR partial files" << std::endl;
            exit(1);
        }
        inputPrevBCR = argv[3];
        std::cerr << "WARNING!: BUILD_BCR_FROM_BCRpartials==1 - The types must be equal to previous BCR partial files (see BUILD_BCR_FROM_BCRpartials in Parameters.h)" << endl;
    #endif
    
	
	std::cout << "BWTCollection: The input is " << argv[1] << std::endl;
	ifstream f(argv[1]);
    	if (f.good() == 0) {
		std::cerr  << argv[1] << " does not exist." << std::endl;
        	exit(1);
	} 
	std::cout << "BWTCollection: The output is " << argv[2] << std::endl;
	#if BUILD_BCR_FROM_BCRpartials==1
		std::cout << "BWTCollection: inputPrevBCR is set to " << inputPrevBCR << std::endl;
		ifstream f2(argv[2]);
		if (f2.good() == 0) {
			std::cerr  << argv[2] << " does not exist." << std::endl;
			exit(1);
		} 
	#endif

	#if ( FASTQ==1) 
        if (USE_QS==0) {
            std::cerr << "Error: if FASTQ==1 then USE_QS==1 (see USE_QS in Tools.h)" << endl;
            exit (EXIT_FAILURE);
        }
    #endif

    #if ( (USE_QS==1) &&  (FASTQ==0) )
        std::cout << "The title of the quality score file must start with '{'. Please use the script.\n";
    #endif

    #if BCR_SET == 0        //Build BCR for 1 sequence
        std::cout << "BCR of 1 sequence" << endl;

        if ( (BUILD_DA_bit == 1) || (BUILD_DA == 1)  )  {
            std::cout << "Error! Since BCR_SET == 0, then BUILD_DA or BUILD_DA_bit must be set to 0 (see Parameters.h).\n";
            exit (EXIT_FAILURE);
        }
        
        #if BCR_FROMCYC==1
            std::cerr << "Error: BCR_FROMCYC==1 (cyc files in input) is not implemented! (see BCR_FROMCYC in Parameters.h)" << endl;
            exit (EXIT_FAILURE);
        #endif
        
        #if USE_QS==1
            std::cerr << "Error: USE_QS==1 (quality score is not implemented! (see USE_QS in Tools.h)" << endl;
            exit (EXIT_FAILURE);
        #endif
            

        #if BUILD_BCR_FROM_BCRpartials==1
            std::cerr << "Error: BUILD_BCR_FROM_BCRpartials==1 is not implemented! (see BUILD_BCR_FROM_BCRpartials in Parameters.h)" << endl;
            exit (EXIT_FAILURE);
        #endif

        #if BCR_FROMCYC==1
            std::cerr << "Error: BCR_FROMCYC==1 is not implemented! (see BCR_FROMCYC in Parameters.h)" << endl;
            exit (EXIT_FAILURE);
        #endif

        #if BCR_INPUT_IN_MEMORY==1      // BCR reads from string
            std::cout << "loads the input file in a string and compute the BWT of the string." << std::endl;
        #else                           // BCR reads from file
            std::cout << "reads from file and computes the BWT of the reverse string." << std::endl;
        #endif
    #else
        std::cout << "BCR of multi-sequences" << endl;
        #if BCR_INPUT_IN_MEMORY==1      // BCR reads from string
            std::cerr << "Error BCR_SET == 0, so that BCR_INPUT_IN_MEMORY must be 0." << std::endl;
            exit (EXIT_FAILURE);
        #endif
        #if BCR_FROMCYC==1 && KSEQ_PARSER==1
            std::cerr << "Error: BCR_FROMCYC==1 (cyc files in input) is not implemented! (Please, set KSEQ_PARSER == 0)" << endl;
            exit (EXIT_FAILURE);
        #endif 
    #endif
    
    
    #if OUTPUT_linear_SuffixArray ==  1
        std::cout << "dataTypeNChar: sizeof(type of #characters): " << sizeof(dataTypeNChar) << "\n";
    #endif


        if ( (BCR_SET!=1) && (BCR_SET_ALN_RH) && (BUILD_BCR_FROM_BCRpartials!=0) && (BUILD_BCR_ALTERNATE!=0) && (BUILD_LCP!=1) )
        {
            std::cerr << "Error! The input is a set. BCR_SET must be set to 1, BUILD_BCR_FROM_BCRpartials must be set to 0, BUILD_BCR_ALTERNATE must be set to 1  (see Parameters.h).\n";
            exit (EXIT_FAILURE);
        }
        
        if(BCR_SET_ALN_RH!=1 && (BUILD_SAP || BUILD_RED_SAP))
        {
            std::cerr << "Error! SAP array is computed only if BCR_SET_ALN_RH is set to 1 (see Parameters.h).\n";
            exit (EXIT_FAILURE);
        }
		
		if(BUILD_SAP && BUILD_RED_SAP)
        {
            std::cerr << "Warning! To compute the reduced SAP array BUILD_SAP must be set to 0 (see Parameters.h).\n";
        }
         
        if (OUTPUT_FORMAT==0) {
            std::cout << "The output format of BCR is at most 5 files (ebwt, lcp, da, posSA, ebwt.qs) - built one after the other.\n";
        }
        
        if (OUTPUT_FORMAT==1) {
            if ((BUILD_LCP==1) && (BUILD_DA==1) && (BUILD_SA==1) )
                std::cout << "The output format of BCR is as the output of eGSA.\n";
            else {
                std::cerr << "Error! The output format of BCR is as the output of EGSA. BUILD_LCP and BUILD_SA and BUILD_DA must be set to 1 (see Parameters.h).\n";
                exit (EXIT_FAILURE);
            }
        }

        if (OUTPUT_FORMAT==2) {
            if ((BUILD_LCP==1) )
                std::cout << "The output format of BCR is a unique file .egsa. Order: bwt, lcp, da, sa.\n";
            else {
                std::cerr << "Error! The output format of BCR is a unique file .egsa (we do not use a struct). BUILD_LCP must be set to 1  (see Parameters.h), BUILD_DA and BUILD_SA could be set to a either 0 or 1.  Order: bwt, lcp, da, sa.\n";
                exit (EXIT_FAILURE);
            }
        }

        if (OUTPUT_FORMAT==3) {
            if ( ((BUILD_LCP == 1) || (BUILD_DA==1) || (BUILD_SA==1) || (BUILD_SAP==1) || (BUILD_RED_SAP==1) || KEEP_eBWT_IN_EXT_MEMORY==1) )
                std::cout << "The output format of BCR is at most 5 files (ebwt, lcp, da, posSA, SAP-array) at the same time.\n";
            else {
                std::cerr << "Error! The output format of BCR is at most 4 files (ebwt, lcp, da, posSA) at the same time\n";
                std::cerr << " BUILD_LCP or BUILD_DA==1 or BUILD_SA==1 or KEEP_eBWT_IN_EXT_MEMORY must be 1 (see Parameters.h).\n";
                exit (EXIT_FAILURE);
            }
        }
        
        if (OUTPUT_FORMAT==4) {
            if (BUILD_DA==1) 
                std::cout << "The output format of BCR is at most 3 files (ebwt, da), lcp, sa.\n";
            else {
                std::cerr << "Error! The output format of BCR is at most 3 files (ebwt, da), lcp, sa. BUILD_DA must be set to 1 (see Parameters.h).\n";
                exit (EXIT_FAILURE);
            }
        }
        
        if (OUTPUT_FORMAT==5) {
            if ( (BUILD_DA==1) && (BUILD_LCP==1) )
                std::cout << "The output format of BCR is at most 3 files ebwt, (lcp, da), sa.\n";
            else {
                std::cerr << "Error! The output format of BCR is at most 3 files ebwt, (lcp, da), sa. BUILD_DA and BUILD_LCP must be set to 1 (see Parameters.h).\n";
                exit (EXIT_FAILURE);
            }
        }
        
        if (OUTPUT_FORMAT==6) {
            if ( (BUILD_DA==1) && (BUILD_SA==1) )
                std::cout << "The output format of BCR is at most 3 files ebwt), lcp, (sa, da).\n";
            else {
                std::cerr << "Error! The output format of BCR is at most 3 files ebwt, lcp, (sa, da). BUILD_DA and BUILD_SA must be set to 1 (see Parameters.h).\n";
                exit (EXIT_FAILURE);
            }
        }

    std::cout << "********************" << endl;
    std::cout << "Setting:\n";
    std::cout << "\tSpecial symbols: " << std::endl;
    std::cout << "\t\tTERMINATE_CHAR: " << TERMINATE_CHAR << " Ascii: " << (int)TERMINATE_CHAR << std::endl;
    std::cout << "\t\tTERMINATE_CHAR_LEN: " << TERMINATE_CHAR_LEN << " Ascii: " << (int)TERMINATE_CHAR_LEN <<  std::endl;

    std::cout << "\tTypes:\n";
    std::cout << "\t\tdataTypedimAlpha: sizeof(type size of alpha): " << sizeof(dataTypedimAlpha) << " bytes \n";
    std::cout << "\t\tdataTypelenSeq: sizeof(type of seq length): " << sizeof(dataTypelenSeq) << " bytes \n";
    std::cout << "\t\tdataTypeNSeq: sizeof(type of #sequences): " << sizeof(dataTypeNSeq) << " bytes\n";
    std::cout << "\t\tdataTypeNChar: sizeof(type of #sequences): " << sizeof(dataTypeNChar) << " bytes\n";
    std::cout << "\t\sortElement: sizeof(type of sortElement): " << sizeof(sortElement) << " bytes\n";
    std::cout << "\t\vectTriple: sizeof(type of vectTriple): " << sizeof(vectTriple) << "bytes\n";
    #if OUTPUT_FORMAT == 1
        std::cerr << "EGSA: sizeof(type of t_GSA): " << sizeof(t_GSA) << "\n";
        std::cerr << "file EGSA: " << lengthTot_plus_eof * sizeof(t_GSA) << "\n";
    #endif


    std::cout << "\tAlignment of the strings in cyc files:\n";
    #if BCR_SET_ALN_RH == 1
        std::cout << "\t\tright" << endl;
    #else
        std::cout << "\t\tleft" << endl;
    #endif
	
    std::cout << "\tAlphabet order:\n";
    #if BUILD_BCR_ALTERNATE == 0
        std::cout << "\t\tLexicographic order" << endl;
    #else
        std::cout << "\t\tAlternate lexicographic order" << endl;
    #endif

    #if BUILD_SAP == 1 || BUILD_RED_SAP == 1        
	std::cout << "\t\tBuild SAP-array and use the folloging ";
    #endif
    std::cout << "\tStrings order:\n";
    #if RLO==0 && SAP_PLUS==0 && SAP_INVERSE==0 && SAP_RANDOM==0
	std::cout << "\t\tInput order" << endl;
    #endif
    #if RLO == 1
	std::cout << "\t\tReverse lexicographic order" << endl;
    #endif
    #if SAP_PLUS == 1
	std::cout << "\t\tSAP-PLUS" << endl;
    #endif
    #if SAP_INVERSE == 1
	std::cout << "\t\tAlternating sap (SAP-ALT) order" << endl;
    #endif
    #if SAP_RANDOM == 1
	std::cout << "\t\tSAP_RANDOM order" << endl;
    #endif
	
    std::cout << "\tEBWT computation (.ebwt file)" << endl;

    #if (BUILD_LCP==1)
        std::cout << "\t\ttogether LCP array (.lcp file)" << endl;
    #endif

    #if (BUILD_SA==1)
        std::cout << "\t\ttogether SA array, i.e. positions in the sequence, so they go from 0 to sequence length (.posSA file)" << endl;
    #endif

    #if (BUILD_DA==1)
        std::cout << "\t\ttogether DA array, i.e. ID/Colour of the sequence (.DA file)" << endl;
    #endif

    #if (BUILD_DA_bit==1) 
        std::cout << "\t\ttogether DA (bit vector) of the sets, i.e. ID/Colour of two sets (.bitDA file)" << endl;
    #endif

    #if USE_QS == 1
        std::cout << "\t\ttogether QS permutation" << endl;
    #endif

    #if BUILD_SAP == 1
        std::cout << "\t\ttogether SAP array" << endl;	
    #elif BUILD_RED_SAP == 1
		std::cout << "\t\ttogether reduced SAP array" << endl;
	#endif
	
    std::cout << "\tBWT partial:\n";
    #if  KEEP_eBWT_IN_EXT_MEMORY == 1    //BCR uses the internal memory for the BWT partial
        std::cout << "\t\tin external memory" << endl;
    #else
        std::cout << "\t\tin internal memory" << endl;
    #endif

    std::cout << "\tInput file:\n";
    #if  BCR_INPUT_IN_MEMORY == 0    //BCR reads cyc files
        std::cout << "\t\tin external memory, uses .cyc files by using ";
        if(BCR_SET_ALN_RH==0 )
            std::cout << "LEFT alignment" ;
        else
            std::cout << "RIGHT alignment" ;
        if(BCR_FROMCYC==0 )
            std::cout << " and builds the cyc files" << endl;
        else
            std::cout << " and using pre-computed cyc files" << endl;
    #else
        std::cout << "\t\tin internal memory (no cyc files)" << endl;
    #endif

    std::cout << "End Setting" << endl;
    std::cout << "********************" << endl;

	BWTCollection *BCRexternalBWT = BWTCollection::InitBWTCollection(argv[1], argv[2], inputPrevBCR );
	
	//cout << "finished iteration, usage: " << timer << endl;
	delete BCRexternalBWT;

	std::cerr << "\nThe BWT et al. is ready! " << std::endl;


	std::cerr << "The End!\n";
}
