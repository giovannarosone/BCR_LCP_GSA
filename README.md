Main author: Giovanna Rosone (giovanna.rosone (at) unipi.it)

# BCR_LCP_GSA

Multi-string BWT (and related data structures) computation
(Last Updated: November 15th 2022)

    BCR is part of:
    BEETL: Burrows-Wheeler Extended Tool Library
    [https://github.com/BEETL]
    Documentation in: doc/BEETL.md
    Copyright (c) 2011-2014 Illumina, Inc.
    BEETL software package is
    covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
       
    By Giovanna Rosone
   
   

BCR_LCP_GSA can compute at the same time:

- the (extended) Burrows-Wheeler transform (multi-string BWT)
- the longest common prefix array (optional, set BUILD_LCP to 1)
- the generalized suffix array (optional):
    - document array (DA[i] corresponds to the ID of the sequence of the symbol ebwt[i]), set BUILD_DA to 1
    - suffix array (SA[i] corresponds to the position of the suffixes of the sequence with id=DA[i] associated to the symbol ebwt[i]), set BUILD_SA to 1. You could not compute the DA array.
- the quality score permutation (see Install)
- the SAP-array and BWT-RLO (see Install)
    
of a very large collection of strings **having different or same length** and **any alphabet**. 

The output format for eBWT/LCP/GSA(DA/SA) can be:
- EGSA (OUTPUT_FORMAT must be 1). The end-marker in .bwt file is the symbol '\0' [https://github.com/felipelouza/egsa]
- at most four files: .ebwt, .lcp, .da, posSA (OUTPUT_FORMAT != 1). The end-marker in .bwt file is the symbol '#'
    - if OUTPUT_FORMAT == 0, the output format of BCR is at most 4 files - built one after the other
    - if OUTPUT_FORMAT == 1, the output format of BCR is as the output of EGSA (.gesa file). BUILD_LCP, BUILD_DA and BUILD_SA must be set to 1. Please, set the types as in eGSA
    - if OUTPUT_FORMAT == 2, the output format of BCR is a unique file .egsa. BUILD_LCP must be set to 1 (we do not use a struct), BUILD_DA and BUILD_SA could be set to either 0 or 1.  Order: ebwt, lcp, da, sa
    - if OUTPUT_FORMAT == 3, the output format of BCR is at most 4 files at the same time
    - if OUTPUT_FORMAT == 4, the output format of BCR is at most 3 files (ebwt, da), lcp, sa
    - if OUTPUT_FORMAT == 5, the output format of BCR is at most 3 files ebwt, (lcp, da), sa
    - if OUTPUT_FORMAT == 6, the output format of BCR is at most 3 files ebwt, lcp, (sa, da)


One should set the following size of the buffers. 
The two steps are consecutive, so you must check that their maximum does not exceed the total available memory.

Please check BUFFERSIZE in TransposeFasta.h. It should be set based on the total available memory.

It means that BCR needs (max length of the reads) * BUFFERSIZE bytes in order to build the cyc files.

Please check SIZEBUFFER in parameters.h. It should be set based on the total available memory:

It means that BCR could need to (sizeof(bwt[i]) + sizeof(lcp[i]) + sizeof(da[i]) + sizeof(sa[i])) * SIZEBUFFER bytes in order to keep the I/O buffers.

You also keep in mind that it could use (sizeof(bwt[i]) + 2*sizeof(lcp[i]) + sizeof(da[i]) + sizeof(x)) bytes per text, where x is the type necessary to index the characters of the (complete) ebwt.


Default:
- Build eBWT permutation
- Types: uchar (you have to reserve two characters) for eBWT, uchar (max length of sequences 256-1) for LCP (BUILD_LCP==1)  
- Store distinct partial files at the same time: OUTPUT_FORMAT == 3

For instance, the data structures LCP and SA depend on the setting of dataTypeLengthSequences/dataTypelenSeq.
If your dataset contains sequences having a length greater than 256-1, you should set dataTypeLengthSequences to 2, so that dataTypelenSeq is set to uint. 

BCR supports the inserting and deleting of new elements belonging to a sequence in the computed data structures. 
Now, BCR can add the elements of a new string collection in input by starting with the previously computed data structures. 
In this case, BCR takes in input the collection (.fasta) that to be inserted and the BCR partial files.
In order to do this, BUILD_BCR_FROM_BCRpartials must to be set to 1.
For instance:
```sh
./BCR_LCP_GSA test/2seqsVar.fa test/2seqsVar.fa.out test/part_7seqsVar
```
where part_7seqsVar is the prefix of the filenames of BCR partial files.

In order to get BCR partial files from eGSA [https://github.com/felipelouza/egsa], one can use the script eGSA_2_BCR_partials.
```sh
./eGSA_2_BCR_partials 7seqsVar.gesa part_7seqsVar
```


If BCR suddenly stops working, you could use the script to remove the cyc files.
```sh
./removeCycFile path #CycFiles
```

### Install

```sh
git clone https://github.com/giovannarosone/BCR\_LCP\_GSA
cd BCR_LCP_GSA
```
Open parameters.h file and, please, set the parameters (data structured that must be computed, types, and so on...).

```sh
make
```

If the input is a multi-line fastQ file (or a gz file containing a fastQ file) and you also want to build the permutation of the QS values then you should compile like this
```sh
make FASTQ=1
```
If you have a single-line fastQ file, you could use:
```sh
Fold long FASTA/Q lines and remove FASTA/Q comments:
./seqtk seq -l100 in.fa > out.fa
(https://github.com/lh3/seqtk.git)
```

If you want to compute the SAP-array along with the multi-string BWT of a string collection implicitly sorted in reverse lexicographical order (BWT-RLO), please compile using
```sh
make SAP=1
```

### Run
```sh
./BCR_LCP_GSA inputFile outputFile
```

### Examples
```sh
./BCR_LCP_GSA test/7seqsVar.fa test/7seqsVar.fa.out
```

```sh
./BCR_LCP_GSA test/test.fq test/test.fq.out
```

```sh
./BCR_LCP_GSA test/test.fq.gz test/test.fq.out
```


#### References:

    *** BWT and LCP
    
    Markus J. Bauer, Anthony J. Cox, Giovanna Rosone 
    Lightweight algorithms for constructing and inverting the BWT of string collections. 
    Theoretical Computer Science 483: 134-148 (2013)
     
    Anthony J. Cox, Fabio Garofalo, Giovanna Rosone, Marinella Sciortino
    Lightweight LCP construction for very large collections of strings. 
    Journal of Discrete Algorithms (2016) 37: 17-33
    
    *** Alternating BWT
    
    Raffaele Giancarlo, Giovanni Manzini, Giovanna Rosone, Marinella Sciortino: 
    A new class of searchable and provably highly compressible string transformations. 
    CPM 2019. Leibniz International Proceedings in Informatics (LIPIcs), 128, art. no. 12. 
    Schloss Dagstuhl- Leibniz-Zentrum fur Informatik GmbH, Dagstuhl Publishing.
    
    Raffaele Giancarlo, Giovanni Manzini, Antonio Restivo, Giovanna Rosone, Marinella Sciortino: 
    The Alternating BWT: An algorithmic perspective. 
    Theoretical Computer Science (2020), Volume 812, Pages 230-243, Elsevier B.V. 
    doi: 10.1016/j.tcs.2019.11.002.
    
    Ira M. Gessel, Antonio Restivo, Christophe Reutenauer,
    A bijection between words and multisets of necklaces, 
    European Journal of Combinatorics, Volume 33, Issue 7, 2012, Pages 1537-1546,
    https://doi.org/10.1016/j.ejc.2012.03.016.
    
    *** QS permutation
    
    Lilian Janin, Giovanna Rosone, and Anthony J. Cox: 
    Adaptive reference-free compression of sequence quality scores. 
    Bioinformatics (2014) 30 (1): 24-30, 
    doi: 10.1093/bioinformatics/btt257.
    
    *** SAP array
    
    Anthony J. Cox, Markus J. Bauer, Tobias Jakobi, Giovanna Rosone, 
    Large-scale compression of genomic sequence databases with the Burrows–Wheeler transform, 
    Bioinformatics (2012) 28 (11): 1415–1419, 
    doi: 10.1093/bioinformatics/bts173.

---
<small> Supported by the project Italian MIUR-SIR [CMACBioSeq][240fb5f5] ("_Combinatorial methods for analysis and compression of biological sequences_") grant n.~RBSI146R5L. P.I. Giovanna Rosone</small>

[240fb5f5]: http://pages.di.unipi.it/rosone/CMACBioSeq.html
