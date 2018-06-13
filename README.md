# BCR_LCP_GSA

Multi-string eBWT/LCP/GSA computation

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
- the generalize suffix array  (optional)
- the longest common prefix array (optional)
of a very large collection of strings having different or same length. 


One can set the parameters in Tools.h

The output format can be:
- EGSA (OUTPUT_FORMAT_EGSA must be 1). The end-marker in .bwt file is the symbol '\0' [https://github.com/felipelouza/egsa]
- three files: .bwt, .lcp, .pairSA, (OUTPUT_FORMAT_EGSA must be 0). The end-marker in .bwt file is the symbol '#'


COMPILE

make

RUN

./BCR_LCP_GSA inputFile outputFile 0


References:

    Markus J. Bauer, Anthony J. Cox and Giovanna Rosone
    Lightweight BWT Construction for Very Large String Collections.
    Proceedings of CPM 2011, pp.219-231
    
    Markus J. Bauer, Anthony J. Cox, Giovanna Rosone and Marinella Sciortino
    Lightweight LCP Construction for Next-Generation Sequencing Datasets. 
    Proceedings of WABI 2012, pp 326-337, 2012
 
    Markus J. Bauer, Anthony J. Cox, Giovanna Rosone 
    Lightweight algorithms for constructing and inverting the BWT of string collections. 
    Theoretical Computer Science 483: 134-148 (2013)
     
    Anthony J. Cox, Fabio Garofalo, Giovanna Rosone, Marinella Sciortino
    Lightweight LCP construction for very large collections of strings. 
    Journal of Discrete Algorithms (2016) 37: 17-33



