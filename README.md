Main author: Giovanna Rosone (giovanna.rosone (at) unipi.it)

# BCR_LCP_GSA

Multi-string BWT (and related data structures) computation
(Last Updated: January 24th 2024)

    BCR is part of:
    BEETL: Burrows-Wheeler Extended Tool Library
    [https://github.com/BEETL]
    Documentation in: doc/BEETL.md
    Copyright (c) 2011-2014 Illumina, Inc.
    BEETL software package is
    covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
       
    By Giovanna Rosone

###  Idea

BCR (both the strategy described in the paper and the implementation) uses each string in circular way. 

A distinct end marker symbol is added (implicitly) to each string, and this allows the eBWT to be built by (implicitly) sorting the string suffixes, rather than cyclic rotations.

Note that BCR does not require linearizing the strings in the collection S by concatenating the strings and does not require explicitly computing the GSA.

BCR builds the BWT incrementally in m steps, where m is the length of the longest string in S (including the end-markers) by exploiting the basic properties of BWT and LF mapping.
The strings in S can be align either left or right and BCR scans all the strings in S from right to left at the same time, to incrementally construct the partial BWTs of  suffixes of S with length at most h \leq m (for more details see the paper).

Note that the BCR implementation prints in the output file (eBWT string) the same end-marker symbol, $\$$, for all strings (\ie, $\$_i=\$$ for all $i=0,\ldots,m-1$) so as not to increase the size of the alphabet, but one can decide whether to output also the list of end-marker indices in order of appearance in the BWT (please compile by setting STORE_INDICES_DOLLARS=1) or not (please compile by setting STORE_INDICES_DOLLARS=0).


###  The main features:

Given a very large collection of strings **having different or same length** and **any alphabet**, BCR_LCP_GSA can compute at the same time:
- the (extended) Burrows-Wheeler transform (multi-string BWT)
- the longest common prefix array (optional)
- the generalized suffix array (optional):
    - document array (DA[i] corresponds to the ID of the sequence of the symbol ebwt[i])
    - suffix array (SA[i] corresponds to the position of the suffixes of the sequence with id=DA[i] associated with the symbol ebwt[i]).
- the quality score permutation (see Install)
- the SAP-array (see Install)
- the BWT-RLO (see Install)
 
Other features not described here can be set directly in the Parameters.h file.

Default:
- Build only eBWT string
- Types: **uchar** for symbol in eBWT, **uint** for the length of the strings, *ulong* for the number of symbols in the input strings (length of the BWT) (see [Wiki](https://github.com/giovannarosone/BCR_LCP_GSA/wiki/Type-setting) for a description of how to change data types for your collection)
- Store distinct partial files at the same time: OUTPUT_FORMAT == 3 (see [Wiki](https://github.com/giovannarosone/BCR_LCP_GSA/wiki/Output-format))
  
For details on the output format, please, see [Wiki](https://github.com/giovannarosone/BCR_LCP_GSA/wiki/Output-format). 

One should set the following size of the buffers. This is important when reading the input file, please see Wiki](https://github.com/giovannarosone/BCR_LCP_GSA/wiki/Buffer-size), especially if you get errors.

BCR supports the inserting and deleting of new elements belonging to a sequence in the computed data structures. 
For details on the dynamic eBWT, please, see [Wiki](https://github.com/giovannarosone/BCR_LCP_GSA/wiki/Dynamic-EBWT). 

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

Alternatively, if you want to build the LCP and/or the DA data structures, you could compile using
```sh
make LCP=1 DA=1
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

To compute the reduced SAP-array, please compile using
```sh
make SAP=1
```

To compute the multi-string BWT of the string collection implicitly sorted in *reverse lexicographical order* (BWT-RLO), please compile using
```sh
make RLO=1
```

To output the indexes of the distinct end-marker symbols, , please compile using
```sh
make STORE_INDICES_DOLLARS=1
```

To print the output of BCR in the terminal, please compile using
```sh
make PRINT=1
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

### Some computed data structures
Let S={CG, CGAT, GGGAT, CGCT, AGCT, AGAT, AGGAT, GGCT}.
Note that in #seq appears the index (in [0,7]) of each end-markers in position order in the eBWT string. 

```sh
i	#seq	eBWT	LCP	DA	SAP-array	Sorted suffixes (only for clarity)
1	 	G	0	0	0		$_0
2	 	T	0	1	1		$_1
3	 	T	0	2	1		$_2
4	 	T	0	3	1		$_3
5	 	T	0	4	1		$_4
6	 	T	0	5	1		$_5
7	 	T	0	6	1		$_6
8	 	T	0	7	1		$_7
9	5	$	0	5	0		AGAT$_5
10	4	$	2	4	0		AGCT$_4
11	6	$	2	6	0		AGGAT$_6
12	 	G	1	1	0		AT$_1
13	 	G	2	2	1		AT$_2
14	 	G	2	5	1		AT$_5
15	 	G	2	6	1		AT$_6
16	0	$	0	0	0		CG$_0
17	1	$	2	1	0		CGAT$_1
18	3	$	2	3	0		CGCT$_3
19	 	G	1	3	0		CT$_3
20	 	G	2	4	1		CT$_4
21	 	G	2	7	1		CT$_7
22	 	C	0	0	0		G$_0
23	 	C	1	1	0		GAT$_1
24	 	G	3	2	1		GAT$_2
25	 	A	3	5	1		GAT$_5
26	 	G	3	6	1		GAT$_6
27	 	C	1	3	0		GCT$_3
28	 	A	3	4	1		GCT$_4
29	 	G	3	7	1		GCT$_7
30	 	G	1	2	0		GGAT$_2
31	 	A	4	6	1		GGAT$_6
32	7	$	2	7	0		GGCT$_7
33	2	$	2	2	0		GGGAT$_2
34	 	A	0	1	0		T$_1
35	 	A	1	2	1		T$_2
36	 	C	1	3	1		T$_3
37	 	C	1	4	1		T$_4
38	 	A	1	5	1		T$_5
39	 	A	1	6	1		T$_6
40	 	C	1	7	1		T$_7
```


### Wiki

See more details and additional features in [Wiki](https://github.com/giovannarosone/BCR_LCP_GSA/wiki). 


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
    doi:10.1016/j.ejc.2012.03.016.
    
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

    Davide Cenzato, Veronica Guerrini, Zsuzsanna Lipták, Giovanna Rosone,
    Computing the optimal BWT of very large string collections, 
    DCC 2023, Volume 2023-March, Pages 71-80,
    doi: 10.1109/DCC55655.2023.00015

---
<small> Supported by the project Italian MIUR-SIR [CMACBioSeq][240fb5f5] ("_Combinatorial methods for analysis and compression of biological sequences_") grant n.~RBSI146R5L. P.I. Giovanna Rosone</small>

[240fb5f5]: http://pages.di.unipi.it/rosone/CMACBioSeq.html
