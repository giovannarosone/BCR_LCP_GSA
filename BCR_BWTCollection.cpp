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
    #endif
    
	
	std::cout << "BWTCollection: The input is " << argv[1] << std::endl;
	std::cout << "BWTCollection: The output is " << argv[2] << std::endl;
	#if BUILD_BCR_FROM_BCRpartials==1
		std::cout << "BWTCollection: inputPrevBCR is set to " << inputPrevBCR << std::endl;
	#endif


	BWTCollection *BCRexternalBWT = BWTCollection::InitBWTCollection(argv[1], argv[2], inputPrevBCR );
	
	//cout << "finished iteration, usage: " << timer << endl;
	delete BCRexternalBWT;

	std::cerr << "\nThe BWT et al. is ready! " << std::endl;


	std::cerr << "The End!\n";
}
