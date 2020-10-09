# BCR_LCP_GSA

Multi-string eBWT/LCP/GSA(DA/SA) computation
(Last Updated: December 9th 2019)

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
- the generalize suffix array (optional):
    - document array (DA[i] corresponds to the ID of the sequence of the symbol ebwt[i]), set BUILD_DA to 1
    - suffix array (SA[i] corresponds to the position of the suffixes of the sequence with id=DA[i] associated to the symbol ebwt[i]), set BUILD_SA to 1. You could not compute the DA array.
- the quality score permutation
    
of a very large collection of strings having different or same length. 

The output format can be:
- EGSA (OUTPUT_FORMAT must be 1). The end-marker in .bwt file is the symbol '\0' [https://github.com/felipelouza/egsa]
- at most four files: .ebwt, .lcp, .da, posSA (OUTPUT_FORMAT != 1). The end-marker in .bwt file is the symbol '#'
    - if OUTPUT_FORMAT == 0, the output format of BCR is at most 4 files - built one after the other
    - if OUTPUT_FORMAT == 1, the output format of BCR is as the output of EGSA (.gesa file). BUILD_LCP, BUILD_DA and BUILD_SA must be set to 1. Please, set the types as in eGSA
    - if OUTPUT_FORMAT == 2, the output format of BCR is a unique file .egsa. BUILD_LCP must be set to 1 (we do not use a struct), BUILD_DA and BUILD_SA could be set to a either 0 or 1.  Order: ebwt, lcp, da, sa
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
- Build eBWT/QS permutation.
- Types: uchar for eBWT, uchar (max length of sequences 256) for LCP (BUILD_LCP==1)  
- Store distinct partial files at the same time: OUTPUT_FORMAT == 3

For instance, the data structures LCP and SA depend on the setting of dataTypeLengthSequences/dataTypelenSeq.
If your dataset contains sequences having a length greater than 256, you should set dataTypeLengthSequences to 2, so that dataTypelenSeq is set to uint. 

In [Bauer el al, 2013] and [Cox et al, 2016], we explain the BCR supports the inserting and deleting of new elements belonging to a sequence in the computed data structures. 
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


---
<small> Supported by the project Italian MIUR-SIR [CMACBioSeq][240fb5f5] ("_Combinatorial methods for analysis and compression of biological sequences_") grant n.~RBSI146R5L. P.I. Giovanna Rosone</small>

[240fb5f5]: http://pages.di.unipi.it/rosone/CMACBioSeq.html
