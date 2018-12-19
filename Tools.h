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
 
 /*
 * Collection of basic tools and defines
 */

#ifndef _TOOLS_H_
#define _TOOLS_H_

#include <iostream>
#include <fstream>
#include "Parameters.h"

//#include <ctime>

#define BackByVector 1


#if ((BUILD_BCR_FROM_EGSA == 1) || (OUTPUT_FORMAT_EGSA == 1))

	typedef struct{
		dataTypelenSeq	suff;
		dataTypelenSeq 	lcp;
		uchar		bwt;
		dataTypeNSeq	text;
	} t_GSA;
	
	
/* 	struct __attribute__((__packed__)) t_GSA {
		dataTypeNSeq	text;
		dataTypelenSeq	suff;
		dataTypelenSeq 	lcp;

		uchar		bwt;
	}; */
	

#endif


//struct __attribute__((__packed__)) ElementTypeEBWTda {
struct ElementTypeEBWTda {
	uchar bwt;          //It is the position in the sequence, so it goes from 0 a length read
	dataTypeNSeq da;		//It is the number of the sequence.
};

//struct __attribute__((__packed__)) ElementTypeLCPda {
struct ElementTypeLCPda {
	dataTypelenSeq lcp;          //It is the position in the sequence, so it goes from 0 a length read
	dataTypeNSeq da;		//It is the number of the sequence.
};

//struct __attribute__((__packed__)) ElementTypeGSA {
struct ElementTypeGSA {
	dataTypelenSeq sa;          //It is the position in the sequence, so it goes from 0 a length read
	dataTypeNSeq da;		//It is the number of the sequence.
};


class Tools
{
private:
    static double startTime;
public:
    static void StartTimer();
    static double GetTime();
    static uchar * GetRandomString(unsigned, unsigned, unsigned &);
 	static uchar * GetFileContents(char *, ulong =0);
    static unsigned FloorLog2(ulong);
    static unsigned CeilLog2(ulong);
    static unsigned* MakeTable();
    static unsigned FastFloorLog2(unsigned);
};

#endif
