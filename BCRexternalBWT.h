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

#ifndef _BCRexternalBWT_H_
#define _BCRexternalBWT_H_

#include <map>
#include "BWTCollection.h"
#include "Parameters.h" // Defines ulong and uchar.
#include <fstream>
#include <iostream>

//2022-10-28
#include <algorithm>
#include "Sorting.h"

class BCRexternalBWT : public SXSI::BWTCollection {
public:
    /**
     * Constructor
     */
        explicit BCRexternalBWT(char*, char*, string);

	~BCRexternalBWT();

	int buildBCR(char const *, char const *, char const *, string);

  void storeBWTFilePartial(uchar const *, dataTypelenSeq);
	#if BUILD_LCP == 1
	    void storeBWTandLCP(uchar const *, dataTypelenSeq);
		void storeEntireLCP(const char*);
    #endif
	void storeEntireBWTFilePartial(const char*);
	
	#if ((BUILD_LCP == 1) || (BUILD_DA==1) || (BUILD_SA==1) || KEEP_eBWT_IN_EXT_MEMORY==1)
		virtual int storeEGSAcomplete( const char* );
	#endif

	#if OUTPUT_FORMAT == 1
		virtual int storeEGSAoutputFromEntireFiles (string );
	#endif
	#if ( (BUILD_DA==1) || (BUILD_SA==1) )
		void storeEntirePairSA(const char*);
		//void storeEntireSAfromPairSA(const char*);
	#endif
	
	#if KEEP_eBWT_IN_EXT_MEMORY==0
		void storeEntireBWTIntMem(const char*);
		void  storeBWTIntMem(uchar const *, dataTypelenSeq) ;
		dataTypeNChar rankManySymbolsIntMem(dataTypedimAlpha , dataTypeNChar *,  dataTypeNChar, dataTypeNChar , uchar *);
	#endif
	
	dataTypeNChar rankManySymbolsFilePartial(FILE &, dataTypeNChar *, dataTypeNChar, uchar *);
			
	#if RLO==1
		void sapSort(std::vector<sortElement> &v, dataTypeNSeq start, dataTypeNSeq end);
		static bool cmpSapSort (sortElement a,sortElement b);
	#endif
	
private:
	#if BCR_SET_ALN_RH ==1
		dataTypeNSeq numToRemove=0;
		dataTypeNSeq contToRemove=0;
	#endif
	
	#if BUILD_BCR_FROM_BCRpartials == 1
		dataTypeNChar readPreviousBCR(string);
	#endif
	void printSegments();
	void printOutput(char *);
	void InsertNsymbols(uchar const *, dataTypelenSeq);
	void InsertFirstsymbols(uchar *); //Added/Modified/Removed 2016-02-23


	int createFilePartialBWT();
	FILE * openWriteFilePartialBWT_0();
	
	dataTypeNChar writeFilePartial(uchar * , FILE * ) ;
	FILE * openFilePartialIn( dataTypedimAlpha );
	FILE * openFilePartialOut(dataTypedimAlpha );
	int closeFilePartial(FILE * InFile);
	int renameFilePartial(dataTypedimAlpha );
	dataTypeNChar readOnFilePartial(uchar *, dataTypeNChar , FILE * );
	dataTypeNChar writeOnFilePartial(uchar *, dataTypeNChar , FILE * );
	dataTypeNChar writeSymbolOnFilePartial(uchar , dataTypeNChar , FILE * );
};

#endif
