CC = g++


FASTQ = 0
SAP = 0
RLO = 0
SAP_INVERSE = 0
SAP_PLUS = 0
SAP_RANDOM = 0

LCP = 0
DA = 0

DEFINES = -DFASTQ=$(FASTQ) -DSAP=$(SAP) -DRLO=$(RLO)  -DSAP_INVERSE=$(SAP_INVERSE) -DSAP_PLUS=$(SAP_PLUS) -DSAPRANDOM=$(SAP_RANDOM) -DLCP=$(LCP) -DDA=$(DA) 

CPPFLAGS = -Wall -ansi -pedantic -g -O3 -std=c++11 $(DEFINES)


BCR_BWTCollection_obs = BCR_BWTCollection.o BWTCollection.o BCRexternalBWT.o Tools.o Sorting.o TransposeFasta.o Timer.o -lz
BCR_BWTCollection: $(BCR_BWTCollection_obs)
	$(CC) -o BCR_LCP_GSA $(BCR_BWTCollection_obs)

clean:
	rm -f core *.o *~ BCR_LCP_GSA

depend:
	$(CC) -MM *.cpp *.c > dependencies.mk

include dependencies.mk
