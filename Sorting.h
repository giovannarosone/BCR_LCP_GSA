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
 
#ifndef SORTED_INCLUDED
#define SORTED_INCLUDED
		
#include <stdlib.h>
#include <vector>
#include "Parameters.h" // Defines ulong and uchar.
#include "Tools.h"

#if BUILD_SAP==1 || BUILD_RED_SAP==1 || RLO==1 || SAP_INVERSE || SAP_PLUS || SAP_RANDOM
		//struct __attribute__((__packed__)) sortElement {
		struct sortElement {	
		sortElement() {};
		#if BUILD_LCP == 0
			sortElement( dataTypedimAlpha z, dataTypedimAlpha sap_value, dataTypeNChar x,dataTypeNSeq y) { pileN = z; sap = sap_value; posN = x; seqN = y;};
		#else
			sortElement( dataTypedimAlpha z, dataTypedimAlpha sap_value, dataTypeNChar x,dataTypeNSeq y, dataTypelenSeq l1, dataTypelenSeq l2) { pileN = z; sap = sap_value; posN = x; seqN = y; lcpCurN = l1; lcpSucN = l2; };
		#endif	
		~sortElement() {};
		dataTypedimAlpha pileN;
		bool sap;
		#if BUILD_LCP == 1
			dataTypelenSeq lcpCurN;
			dataTypelenSeq lcpSucN;
		#endif
		#if PI_PERM == 1
			#if PI_POS == 1
				dataTypeNSeq piPos; //extension of the triple for permutation array
			#endif
		#endif
		dataTypeNSeq seqN;
		dataTypeNChar posN;

		};
#else
	#if BUILD_LCP == 1
		//struct __attribute__((__packed__)) sortElement {
		struct sortElement {	
		sortElement() {};
		sortElement( dataTypedimAlpha z, dataTypeNChar x,dataTypeNSeq y, dataTypelenSeq l1, dataTypelenSeq l2) { pileN = z; posN = x; seqN = y; lcpCurN = l1; lcpSucN = l2; };
		~sortElement() {};
		dataTypedimAlpha pileN;
		dataTypelenSeq lcpCurN;
		dataTypelenSeq lcpSucN;
		dataTypeNSeq seqN;	
		dataTypeNChar posN;
		};
	#else
		//struct __attribute__((__packed__)) sortElement {
		struct sortElement {	
		sortElement() {};
		sortElement( dataTypedimAlpha z, dataTypeNChar x,dataTypeNSeq y) { pileN = z; posN = x; seqN = y;};
		~sortElement() {};
		dataTypedimAlpha pileN;
		dataTypeNSeq seqN;
		dataTypeNChar posN;
		};
	#endif
#endif

void quickSort(std::vector< sortElement > &v);

#endif
