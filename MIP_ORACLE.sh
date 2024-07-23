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
echo
### Download the NCBI Taxonomy Database
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
mkdir taxonomy && tar -xxvf taxdump.tar.gz -C taxonomy
echo
now=$(date +"%T")
echo "Current time : $now"
echo "__________________________________________"
echo "Downloaded the NCBI Taxonomy Database."
echo
### Extract the FASTA and Taxonomy Mapping Files
blastdbcmd -db $parameterJ/nt -entry all > nt.fna
blastdbcmd -db $parameterJ/nt -entry all -outfmt "%a %T" > nt.fna.taxidmapping
echo
now=$(date +"%T")
echo "Current time : $now"
echo "__________________________________________"
echo "Extracted the FASTA and Taxonomy Mapping Files."
echo
### Create MMSeqs2 Database and taxonomic database
mmseqs createdb nt.fna nt.fnaDB
mmseqs createtaxdb nt.fnaDB tmp --ncbi-tax-dump taxonomy/ --tax-mapping-file nt.fna.taxidmapping
echo
now=$(date +"%T")
echo "Current time : $now"
echo "__________________________________________"
echo "Created MMSeqs2 Database and taxonomic database."
echo
### filter database to only include human sequences
mmseqs filtertaxseqdb nt.fnaDB seqTaxHumanDB --taxon-list 9606
echo
now=$(date +"%T")
echo "Current time : $now"
echo "__________________________________________"
echo "Filtered database to only include human sequences."
echo
### filter database to exclude human sequences
mmseqs filtertaxseqdb nt.fnaDB seqTaxNonHumanDB --taxon-list '!9606'
echo
now=$(date +"%T")
echo "Current time : $now"
echo "__________________________________________"
echo "Filtered database to exclude human sequences."
echo
now=$(date +"%T")
echo "Current time : $now"
echo "__________________________________________"
### Trying taxonomy search with whole_region.txt file
echo
# Creating a database for whole_region.txt
mmseqs createdb whole_region.txt whole_region_DB

mmseqs createtaxdb whole_region_DB tmp --ncbi-tax-dump taxonomy/ --tax-mapping-file /nfs_master/shweta/New_MIPS/shweta_mips/tax_mmseq/nt.fna.taxidmapping
echo
now=$(date +"%T")
echo "Current time : $now"
echo "__________________________________________"
echo "Created a database for whole_region.txt"
echo "__________________________________________"
# Search against the database of human sequences
mmseqs search whole_region_DB seqTaxHumanDB Human_Hits_for_QUERY tmp --search-type 2
echo
now=$(date +"%T")
echo "Current time : $now"
echo "__________________________________________"
echo "Searched against the database of human sequences."
echo "__________________________________________"
### Creating final tabular output
mmseqs convertalis whole_region_DB seqTaxHumanDB Human_Hits_for_QUERY Human_Hits_for_QUERY.tab 
echo
now=$(date +"%T")
echo "Current time : $now"
echo "__________________________________________"
echo "Created final tabular output for human hits."
echo "__________________________________________"
# Search against the database of Non-human sequences
mmseqs search whole_region_DB seqTaxNonHumanDB NonHuman_Hits_for_QUERY tmp --search-type 2 --max-seqs 10 --db-load-mode 2
echo
now=$(date +"%T")
echo "Current time : $now"
echo "__________________________________________"
echo "Searched against the database of Non-human sequences."
echo "__________________________________________"
### Creating final tabular output
mmseqs convertalis whole_region_DB seqTaxNonHumanDB NonHuman_Hits_for_QUERY NonHuman_Hits_for_QUERY.tab 
echo
now=$(date +"%T")
echo "Current time : $now"
echo "__________________________________________"
echo "Created final tabular output for Non-human hits."
echo "__________________________________________"
echo "MMSEQS2 results have been generated."
echo "__________________________________________"
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
