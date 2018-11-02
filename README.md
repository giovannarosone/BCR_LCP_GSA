# BCR_LCP_GSA

Multi-string eBWT/LCP/GSA(DA/SA) computation

    BCR is part of:
    BEETL: Burrows-Wheeler Extended Tool Library
    [https://github.com/BEETL]
    Documentation in: doc/BEETL.md
    Copyright (c) 2011-2014 Illumina, Inc.
    BEETL software package is
    covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
       
    By Giovanna Rosone
   
   
Note:

BCR_LCP_GSA can compute at the same time:

- the (extended) Burrows–Wheeler transform (multi-string BWT)
- the longest common prefix array (optional, set BUILD_LCP to 1)
- the generalize suffix array (optional):
    - document array (DA[i] corresponds to the ID of the sequence of the symbol ebwt[i]), set BUILD_DA to 1
    - suffix array (SA[i] corresponds to the position of the sequence of the suffixes associated to the symbol ebwt[i]), set BUILD_SA to 1
    
of a very large collection of strings having different or same length. 

The output format can be:
- EGSA (OUTPUT_FORMAT_EGSA must be 1). The end-marker in .bwt file is the symbol '\0' [https://github.com/felipelouza/egsa]
- four files: .bwt, .lcp, .da, posSA (OUTPUT_FORMAT_EGSA must be 0). The end-marker in .bwt file is the symbol '#'

One should set the parameters in Tools.h

Default:
- Build eBWT/LCP/DA.
- Types: uchar for eBWT, uchar (max length of sequences 256) for LCP and uint for DA (BUILD_LCP==1 and BUILD_DA==1)  
- Store distinct files: OUTPUT_FORMAT_EGSA == 0

For instance, the data structures LCP and SA depends on the setting of dataTypeLengthSequences/dataTypelenSeq.
If your dataset contains sequences having a length greater than 256, you should set dataTypeLengthSequences to 1, so that dataTypelenSeq is set to uint. 

COMPILE


First: set the parameters (data structured that must be computed, types, and so on...) in parameters.h 


Second: make


RUN

./BCR_LCP_GSA inputFile outputFile 0


Example:

./BCR_LCP_GSA test/7seqsVar.fa test/7seqsVar.fa.out 0


References:

    Markus J. Bauer, Anthony J. Cox and Giovanna Rosone
    Lightweight BWT Construction for Very Large String Collections.
    Proceedings of CPM 2011, pp.219-231
    
    Markus J. Bauer, Anthony J. Cox, Giovanna Rosone and Marinella Sciortino
    Lightweight LCP Construction for Next-Generation Sequencing Datasets. 
    Proceedings of WABI 2012, pp 326-337
 
    Markus J. Bauer, Anthony J. Cox, Giovanna Rosone 
    Lightweight algorithms for constructing and inverting the BWT of string collections. 
    Theoretical Computer Science 483: 134-148 (2013)
     
    Anthony J. Cox, Fabio Garofalo, Giovanna Rosone, Marinella Sciortino
    Lightweight LCP construction for very large collections of strings. 
    Journal of Discrete Algorithms (2016) 37: 17-33



