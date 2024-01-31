
#!/bin/bash

dataset_directory=$1;
tools_directory=$2;
output_directory=$3;


dataset_array=("ERR7320.fasta")

for dataset in ${dataset_array[@]}; do
  (/usr/bin/time -v ./BCR_LCP_GSA_original ${dataset_directory}/${dataset} ${output_directory}/${dataset}_output/original/${dataset}.out) &> ${output_directory}/${dataset}_output/original/data.txt
  ${tools_directory}/number-runs ${output_directory}/${dataset}_output/original/${dataset}.out.ebwt >> ${output_directory}/${dataset}_output/original/data.txt
  rm ${output_directory}/${dataset}_output/original/${dataset}.*

  (/usr/bin/time -v ./BCR_LCP_GSA_RLO ${dataset_directory}/${dataset} ${output_directory}/${dataset}_output/RLO/${dataset}.out) &> ${output_directory}/${dataset}_output/RLO/data.txt
  ${tools_directory}/number-runs ${output_directory}/${dataset}_output/RLO/${dataset}.out.ebwt >> ${output_directory}/${dataset}_output/RLO/data.txt
  rm ${output_directory}/${dataset}_output/RLO/${dataset}.*

  (/usr/bin/time -v ${tools_directory}/BEETL/install/bin/beetl bwt --algorithm=ext --sap-ordering --output-format ASCII --intermediate-format ASCII --input ${dataset_directory}/${dataset} --output ${output_directory}/${dataset}_output/beetle/${dataset}.out) &> ${output_directory}/${dataset}_output/beetle/data.txt 
  cat ${output_directory}/${dataset}_output/beetle/${dataset}.out-B00 ${output_directory}/${dataset}_output/beetle/${dataset}.out-B01 ${output_directory}/${dataset}_output/beetle/${dataset}.out-B02 ${output_directory}/${dataset}_output/beetle/${dataset}.out-B03 ${output_directory}/${dataset}_output/beetle/${dataset}.out-B04 ${output_directory}/${dataset}_output/beetle/${dataset}.out-B05 >  ${output_directory}/${dataset}_output/beetle/${dataset}.out.ebwt 
  ${tools_directory}/number-runs ${output_directory}/${dataset}_output/beetle/${dataset}.out.ebwt >> ${output_directory}/${dataset}_output/beetle/data.txt
  rm ${output_directory}/${dataset}_output/beetle/${dataset}.*

  (/usr/bin/time -v ./BCR_LCP_GSA_inverse ${dataset_directory}/${dataset} ${output_directory}/${dataset}_output/inverse/${dataset}.out) &> ${output_directory}/${dataset}_output/inverse/data.txt
  ${tools_directory}/number-runs ${output_directory}/${dataset}_output/inverse/${dataset}.out.ebwt >> ${output_directory}/${dataset}_output/inverse/data.txt
  rm ${output_directory}/${dataset}_output/inverse/${dataset}.*

  (/usr/bin/time -v ./BCR_LCP_GSA_plus ${dataset_directory}/${dataset} $output_directory/${dataset}_output/plus/${dataset}.out) &>  $output_directory/${dataset}_output/plus/data.txt
  ${tools_directory}/number-runs ${output_directory}/${dataset}_output/plus/${dataset}.out.ebwt >> ${output_directory}/${dataset}_output/plus/data.txt

  ./optimalBWT/permute $output_directory/${dataset}_output/plus/${dataset}.out.ebwt $output_directory/${dataset}_output/plus/${dataset}.out.bwt.sap 10
  ${tools_directory}/number-runs ${output_directory}/${dataset}_output/plus/${dataset}.out.ebwt.optbwt > ${output_directory}/${dataset}_output/optimal/data.txt

  rm ${output_directory}/${dataset}_output/plus/${dataset}.*
done

