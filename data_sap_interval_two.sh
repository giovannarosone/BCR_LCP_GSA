#!/bin/bash

dataset_directory=$1;
output_directory=$2;

dataset_array=("SRR065390_1.fasta")
#dataset_array=("10seqsVar.fa" "7seqsVar.fa" "8seq.fasta")

echo -e "Data sap interval di tipo 2\n" &> ${output_directory}/data_sap_interval_type_two

for dataset in ${dataset_array[@]}; do

    echo -e "\n" >> ${dataset_directory}/${dataset}.out.bwt.red_sap
    echo -e "\n" >> ${dataset_directory}/${dataset}.out.ebwt
    echo -e "$dataset" >> ${output_directory}/data_sap_interval_type_two

    sap_number=0
    found_sap=0

    sap_length=0
    total_sap_length=0

    diff_char_sap=0
    total_diff_char=0

    string=''

    while read -r -n10000000 sap_line && read -r -n10000000 bwt_line <&3; do

        #echo "$sap_line $bwt_line"
        for (( i=0; i<${#sap_line}; i++ )); do

            if [ ${sap_line:$i:1} -eq 1 ]; then

                if [ $found_sap -eq 0 ]; then
                    ((sap_number++))
                    echo "trovato sap $sap_number"
                    string="$character"
                    #echo "metto $character in string: $string"
                    found_sap=1
                fi

                ((sap_length++))

                c=${bwt_line:$i:1}
                if [[ $string != *$c* ]]; then
                    string="$string$c"
                    #echo "metto $c in string: $string"
                fi
            else 

                if [ $found_sap -eq 1 ]; then
                    total_sap_length=$((total_sap_length+sap_length+1))
                    #echo "lunghezza sap: $sap_length. totale: $total_sap_length"
                    sap_length=0
                    diff_char_sap=${#string}
                    total_diff_char=$((total_diff_char+diff_char_sap))
                    #echo "caratteri distinti sap $diff_char_sap. totale: $total_diff_char"
                    string=''
                    found_sap=0
                fi

                character=${bwt_line:$i:1}
            fi

        done
    done < ${dataset_directory}/${dataset}.out.bwt.red_sap 3< ${dataset_directory}/${dataset}.out.ebwt

    #echo "$sap_line $bwt_line"

    if [ $found_sap -eq 1 ]; then
        total_sap_length=$((total_sap_length+sap_length+1))
        #echo "lunghezza sap: $sap_length. totale: $total_sap_length"
        diff_char_sap=${#string}
        total_diff_char=$((total_diff_char+diff_char_sap))
        #echo "caratteri distinti sap $diff_char_sap. totale: $total_diff_char"
    fi

    echo -e "number of sap intervals: $sap_number\ntotal size of sap intervals: $total_sap_length\ntotal number of different caharacters in sap intervals: $total_diff_char\n\n" >> ${output_directory}/data_sap_interval_type_two
done
