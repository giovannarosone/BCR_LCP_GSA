CC = g++


FASTQ = 0  #Use fastQ file and handle quality score sequences
SAP = 0    #compute the reduced SAP array
RLO = 0    #compute the RLO-BWT

LCP = 0    #compute the LCP array
DA = 0     #compute the DA array

STORE_INDICES_DOLLARS = 1 #Store the indexes of each distint end-marker in eBWT string 
#For each end-marker, the binary file contains
#- indexes of the sequences
#- position in the eBWT string
#- symbol with which the associated suffix begins, i.e. the first symbol in the string.

PRINT = 0 #print in the terminal the BWT along with other data structures (if computed), such as LCP, DA, GSA, SAP array, and so on
#if STORE_INDICES_DOLLARS = 1 and PRINT = 1
# print the index of each end-markers symbol

DEFINES = -DFASTQ=$(FASTQ) -DSAP=$(SAP) -DRLO=$(RLO) -DLCP=$(LCP) -DDA=$(DA) -DSTORE_ENDMARKER_POS=$(STORE_INDICES_DOLLARS) -DprintFinalOutput=$(PRINT)

CPPFLAGS = -Wall -ansi -pedantic -g -O3 -std=c++11 $(DEFINES)


BCR_BWTCollection_obs = BCR_BWTCollection.o BWTCollection.o BCRexternalBWT.o Tools.o Sorting.o TransposeFasta.o Timer.o -lz
BCR_BWTCollection: $(BCR_BWTCollection_obs)
	$(CC) -o BCR_LCP_GSA $(BCR_BWTCollection_obs)

clean:
	rm -f core *.o *~ BCR_LCP_GSA

depend:
	$(CC) -MM *.cpp *.c > dependencies.mk

include dependencies.mk
