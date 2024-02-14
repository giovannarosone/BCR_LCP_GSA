
#!/bin/bash

dataset_directory=$1;
tools_directory=$2;
output_directory=$3;
optimal_directory=$4;


dataset_array=("10seqsVar.fa" "7seqsVar.fa")

echo -e "ORIGINAL\ndataset,%CPU,WALL_CLOCK,RAM,runs" &> ${output_directory}/table_original.txt
echo -e "RLO\ndataset,%CPU,WALL_CLOCK,RAM,runs" &> ${output_directory}/table_RLO.txt
#echo -e "beetle\ndataset,%CPU,WALL_CLOCK,RAM,runs" &> ${output_directory}/table_beetle.txt
echo -e "INVERSE\ndataset,%CPU,WALL_CLOCK,RAM,runs" &> ${output_directory}/table_inverse.txt
echo -e "PLUS\ndataset,%CPU,WALL_CLOCK,RAM,runs" &> ${output_directory}/table_plus.txt
echo -e "OPTIMAL\ndataset,%CPU,WALL_CLOCK,RAM,runs" &> ${output_directory}/table_optimal.txt
echo -e "SAP_INTERVAL\ndataset,sap_interval,sap_interval_type_two" &> ${output_directory}/table_temp.txt

for dataset in ${dataset_array[@]}; do

  (/usr/bin/time -v ./BCR_LCP_GSA_original ${dataset_directory}/${dataset} ${output_directory}/${dataset}.out) &> ${output_directory}/data_${dataset}.txt
  ${tools_directory}/number-runs ${output_directory}/${dataset}.out.ebwt >> ${output_directory}/data_${dataset}.txt
  rm ${output_directory}/${dataset}.*

  CPU=$(grep "CPU" ${output_directory}/data_${dataset}.txt | cut -f 7 -d " ")
  CLOCK=$(grep "Elapsed" ${output_directory}/data_${dataset}.txt | cut -f 8 -d " ")
  RAM=$(grep "Maximum" ${output_directory}/data_${dataset}.txt | cut -f 6 -d " ") 
  RUNS=$(grep "runs" ${output_directory}/data_${dataset}.txt | cut -f 4 -d " ") 
  echo "${dataset},${CPU},${CLOCK},${RAM},${RUNS}" >> ${output_directory}/table_original.txt
  rm ${output_directory}/data_*

  (/usr/bin/time -v ./BCR_LCP_GSA_RLO ${dataset_directory}/${dataset} ${output_directory}/${dataset}.out) &> ${output_directory}/data_${dataset}.txt
  ${tools_directory}/number-runs ${output_directory}/${dataset}.out.ebwt >> ${output_directory}/data_${dataset}.txt
  rm ${output_directory}/${dataset}.*

  CPU=$(grep "CPU" ${output_directory}/data_${dataset}.txt | cut -f 7 -d " ")
  CLOCK=$(grep "Elapsed" ${output_directory}/data_${dataset}.txt | cut -f 8 -d " ")
  RAM=$(grep "Maximum" ${output_directory}/data_${dataset}.txt | cut -f 6 -d " ") 
  RUNS=$(grep "runs" ${output_directory}/data_${dataset}.txt | cut -f 4 -d " ") 
  echo "${dataset},${CPU},${CLOCK},${RAM},${RUNS}" >> ${output_directory}/table_RLO.txt

  sap_counter=$(grep "total sap interval:" ${output_directory}/data_${dataset}.txt | cut -f 4 -d " ")
  rm ${output_directory}/data_*
  

  #(/usr/bin/time -v ${tools_directory}/BEETL/install/bin/beetl bwt --algorithm=ext --sap-ordering --output-format ASCII --intermediate-format ASCII --input ${dataset_directory}/${dataset} --output ${output_directory}/${dataset}.out) &> ${output_directory}/data_${dataset}.txt 
  #cat ${output_directory}/${dataset}.out-B00 ${output_directory}/${dataset}.out-B01 ${output_directory}/${dataset}.out-B02 ${output_directory}/${dataset}.out-B03 ${output_directory}/${dataset}.out-B04 ${output_directory}/${dataset}.out-B05 >  ${output_directory}/${dataset}.out.ebwt 
  #${tools_directory}/number-runs ${output_directory}/${dataset}.out.ebwt >> ${output_directory}/data_${dataset}.txt
  #rm ${output_directory}/${dataset}.*

  #CPU=$(grep "CPU" ${output_directory}/data_${dataset}.txt | cut -f 7 -d " ")
  #CLOCK=$(grep "Elapsed" ${output_directory}/data_${dataset}.txt | cut -f 8 -d " ")
  #RAM=$(grep "Maximum" ${output_directory}/data_${dataset}.txt | cut -f 6 -d " ") 
  #RUNS=$(grep "runs" ${output_directory}/data_${dataset}.txt | cut -f 4 -d " ") 
  #echo "${dataset},${CPU},${CLOCK},${RAM},${RUNS}" >> ${output_directory}/table_beetle.txt
  #rm ${output_directory}/data_*

  (/usr/bin/time -v ./BCR_LCP_GSA_inverse ${dataset_directory}/${dataset} ${output_directory}/${dataset}.out) &> ${output_directory}/data_${dataset}.txt
  ${tools_directory}/number-runs ${output_directory}/${dataset}.out.ebwt >> ${output_directory}/data_${dataset}.txt
  rm ${output_directory}/${dataset}.*

  CPU=$(grep "CPU" ${output_directory}/data_${dataset}.txt | cut -f 7 -d " ")
  CLOCK=$(grep "Elapsed" ${output_directory}/data_${dataset}.txt | cut -f 8 -d " ")
  RAM=$(grep "Maximum" ${output_directory}/data_${dataset}.txt | cut -f 6 -d " ") 
  RUNS=$(grep "runs" ${output_directory}/data_${dataset}.txt | cut -f 4 -d " ") 
  echo "${dataset},${CPU},${CLOCK},${RAM},${RUNS}" >> ${output_directory}/table_inverse.txt
  rm ${output_directory}/data_*

  (/usr/bin/time -v ./BCR_LCP_GSA_plus ${dataset_directory}/${dataset} $output_directory/${dataset}.out) &>  $output_directory/data_${dataset}.txt
  ${tools_directory}/number-runs ${output_directory}/${dataset}.out.ebwt >> ${output_directory}/data_${dataset}.txt

  CPU=$(grep "CPU" ${output_directory}/data_${dataset}.txt | cut -f 7 -d " ")
  CLOCK=$(grep "Elapsed" ${output_directory}/data_${dataset}.txt | cut -f 8 -d " ")
  RAM=$(grep "Maximum" ${output_directory}/data_${dataset}.txt | cut -f 6 -d " ") 
  RUNS=$(grep "runs" ${output_directory}/data_${dataset}.txt | cut -f 4 -d " ") 
  echo "${dataset},${CPU},${CLOCK},${RAM},${RUNS}" >> ${output_directory}/table_plus.txt
  rm ${output_directory}/data_*

  (/usr/bin/time -v ./BCR_LCP_GSA_sap ${dataset_directory}/${dataset} ${output_directory}/${dataset}.out) &> ${output_directory}/data_${dataset}.txt

  sap_counter_two=$(grep "total sap interval of type two:" ${output_directory}/data_${dataset}.txt | cut -f 7 -d " ")
  echo "${dataset},${sap_counter},${sap_counter_two}" >> ${output_directory}/table_temp.txt

  /usr/bin/time -v ${optimal_directory}/optimalBWT/permute ${output_directory}/${dataset}.out.ebwt ${output_directory}/${dataset}.out.bwt.red_sap 10 &> ${output_directory}/data_optimal_${dataset}.txt
  ${tools_directory}/number-runs ${output_directory}/${dataset}.out.ebwt.optbwt >> ${output_directory}/data_optimal_${dataset}.txt

  CPU=$(grep "CPU" ${output_directory}/data_optimal_${dataset}.txt | cut -f 7 -d " ")
  CLOCK=$(grep "Elapsed" ${output_directory}/data_optimal_${dataset}.txt | cut -f 8 -d " ")
  RAM=$(grep "Maximum" ${output_directory}/data_optimal_${dataset}.txt | cut -f 6 -d " ") 
  RUNS=$(grep "runs" ${output_directory}/data_optimal_${dataset}.txt | cut -f 4 -d " ") 
  echo "${dataset},${CPU},${CLOCK},${RAM},${RUNS}" >> ${output_directory}/table_optimal.txt
  
  rm ${output_directory}/data_*
  rm ${output_directory}/${dataset}.*

done

echo -e "\n" >> ${output_directory}/table_original.txt
echo -e "\n" >> ${output_directory}/table_RLO.txt
echo -e "\n" >> ${output_directory}/table_inverse.txt
echo -e "\n" >> ${output_directory}/table_plus.txt
echo -e "\n" >> ${output_directory}/table_optimal.txt
#echo -e "\n" >> ${output_directory}/table_beetle.txt


cat ${output_directory}/table_original.txt ${output_directory}/table_RLO.txt ${output_directory}/table_inverse.txt ${output_directory}/table_plus.txt ${output_directory}/table_optimal.txt ${output_directory}/table_temp.txt > ${output_directory}/table.txt
rm ${output_directory}/table_*
