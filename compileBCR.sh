
#!/bin/bash

make 
mv BCR_LCP_GSA BCR_LCP_GSA_original
make clean

make RLO=1
mv BCR_LCP_GSA BCR_LCP_GSA_RLO
make clean

make SAP_INVERSE=1
mv BCR_LCP_GSA BCR_LCP_GSA_inverse
make clean

make BCR_LCP_GSA_plus=1
mv BCR_LCP_GSA BCR_LCP_GSA_plus
make clean
