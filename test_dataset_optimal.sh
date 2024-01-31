#!/bin/bash

dataset_directory=$1;
tools_directory=$2;
output_directory=$3;


dataset_array=("SRR12038588.fasta" "ERR7320.fasta" "ERR022075_1.fasta" "Saccharomyces_SRR327342.fasta" "SRR065390_1.fasta" "SRR065390_2.fasta")

for dataset in ${dataset_array[@]}; do

  (/usr/bin/time -v ./BCR_LCP_GSA_RLO ${dataset_directory}/${dataset} ${output_directory}/${dataset}_output/${dataset}.out) &> ${output_directory}/${dataset}_output/RLO.txt
  ${tools_directory}/number-runs ${output_directory}/${dataset}_output/${dataset}.out.ebwt >> ${output_directory}/${dataset}_output/RLO.txt

  (/usr/bin/time -v ./optimalBWT/permute $output_directory/${dataset}_output/${dataset}.out.ebwt $output_directory/${dataset}_output/${dataset}.out.bwt.sap 10) &> ${output_directory}/${dataset}_output/optimal.txt
  ${tools_directory}/number-runs ${output_directory}/${dataset}_output/${dataset}.out.ebwt.optbwt >> ${output_directory}/${dataset}_output/optimal.txt

  rm ${output_directory}/${dataset}_output/${dataset}.*

done
