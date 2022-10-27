CC = g++


FASTQ = 0


DEFINES = -DFASTQ=$(FASTQ) 

CPPFLAGS = -Wall -ansi -pedantic -g -O3 -std=c++11 $(DEFINES)


BCR_BWTCollection_obs = BCR_BWTCollection.o BWTCollection.o BCRexternalBWT.o Tools.o Sorting.o TransposeFasta.o Timer.o -lz
BCR_BWTCollection: $(BCR_BWTCollection_obs)
	$(CC) -o BCR_LCP_GSA $(BCR_BWTCollection_obs)

clean:
	rm -f core *.o *~ BCR_LCP_GSA

depend:
	$(CC) -MM *.cpp *.c > dependencies.mk

include dependencies.mk
