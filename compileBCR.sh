
#!/bin/bash

mkdir out_BCRdata
mkdir out_table

make 
mv BCR_LCP_GSA BCR_inputBWT
make clean

make RLO=1
mv BCR_LCP_GSA BCR_rloBWT
make clean

make SAP_INVERSE=1
mv BCR_LCP_GSA BCR_altBWT
make clean

make SAP_PLUS=1
mv BCR_LCP_GSA BCR_plusBWT
make clean

make SAP_RANDOM=1
mv BCR_LCP_GSA BCR_randBWT
make clean

make SAP=1
mv BCR_LCP_GSA BCR_sapArray
make clean

make LCP=1
mv BCR_LCP_GSA BCR_LCP
make clean

make DA=1
mv BCR_LCP_GSA BCR_DA
make clean

make LCP=1 DA=1
mv BCR_LCP_GSA BCR_LCP_DA
make clean

g++ stastisticsSAPintervals.cpp -o stastisticsSAPintervals
