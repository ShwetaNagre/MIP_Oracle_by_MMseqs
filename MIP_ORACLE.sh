#!/bin/bash

helpFunction()
{
   echo "NOTE: The FASTA definition line should match the CARD database format. Ex.>gb|EU555534|+|0-882|ARO:3002316|KPC-6 [Klebsiella pneumoniae]"
   echo "Usage: $0 -i KPC -o KPC_final_results -l mip_oracle -j /DATA/databases/blast/nt"
   echo -e "\t-i Name of the input FASTA file(There's no need to add the file extension)"
   echo -e "\t-o Name of the ouptut file(There's no need to add the file extension)"
   echo -e "\t-l The name of conda environment containing all the packages"
   echo -e "\t-j The location of the nt BLAST database"
   exit 1 # Exit script after printing help
}

while getopts "i:o:l:j:" opt
do
   case "$opt" in
      i ) parameterI="$OPTARG" ;;
      o ) parameterO="$OPTARG" ;;
      l ) parameterL="$OPTARG" ;;
      j ) parameterJ="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$parameterI" ] || [ -z "$parameterO" ] || [ -z "$parameterL" ] || [ -z "$parameterJ" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi
eval "$(conda shell.bash hook)"

echo "input filename=$parameterI.fasta"
echo "output filename=$parameterO.xlsx"
echo "Working conda environment: $parameterL"
echo "BLAST nr database path: $parameterJ"
now=$(date +"%T")
echo "Current time : $now"
echo "-------------------------------------"
echo "MIP ORACLE"
echo "-------------------------------------"
# Begin script in case all parameters are correct

conda activate $parameterL

python MIP_1.py $parameterI
echo
echo
now=$(date +"%T")
echo "Current time : $now"
echo "__________________________________________"
echo "All possible MIPs have been generated according to the target region length and arm lengths provided in the config file."
python MIP_2.py $parameterI
echo
now=$(date +"%T")
echo "Current time : $now"
echo "__________________________________________"
echo "The MIPs have been filtered for GC content, Temperature, and PolyNs as per the config file."
python MIP_3.py
echo
now=$(date +"%T")
echo "Current time : $now"
echo "__________________________________________"
echo "MMSEQS2 input files have been generated. Running MMSEQS2 now."
mmseqs easy-search whole_region.txt $parameterJ/mmseqs mmseqs2_result tmp --search-type 2
echo
now=$(date +"%T")
echo "Current time : $now"
echo "__________________________________________"
echo "MMSEQS2 results have been generated."
python MIP_4.py
echo
now=$(date +"%T")
echo "Current time : $now"
echo "__________________________________________"
echo "Retrieved the ids for sequences with no hits in the database"
python MIP_5.py
echo
now=$(date +"%T")
echo "Current time : $now"
echo "__________________________________________"
echo "The MIPs sequences have been obtained on the basis of no hits ids"
echo
echo"Current time : $now"
echo "-------------------------------------"
echo "MIP ORACLE"
echo "-------------------------------------"
