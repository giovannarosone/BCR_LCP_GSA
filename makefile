CC = g++


FASTQ = 0
SAP = 0
RLO = 0

LCP = 0
DA = 0

STORE_ENDMARKER_POS = 0

DEFINES = -DFASTQ=$(FASTQ) -DSAP=$(SAP) -DRLO=$(RLO) -DLCP=$(LCP) -DDA=$(DA) -DSTORE_ENDMARKER_POS=$(STORE_ENDMARKER_POS)

CPPFLAGS = -Wall -ansi -pedantic -g -O3 -std=c++11 $(DEFINES)


BCR_BWTCollection_obs = BCR_BWTCollection.o BWTCollection.o BCRexternalBWT.o Tools.o Sorting.o TransposeFasta.o Timer.o -lz
BCR_BWTCollection: $(BCR_BWTCollection_obs)
	$(CC) -o BCR_LCP_GSA $(BCR_BWTCollection_obs)

clean:
	rm -f core *.o *~ BCR_LCP_GSA

depend:
	$(CC) -MM *.cpp *.c > dependencies.mk

include dependencies.mk
