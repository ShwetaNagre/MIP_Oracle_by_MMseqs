# MIP_ORACLE
A software to filter and identify unique target regions with diagnostic significance in antimicrobial resistance genes, and various other pathogen genomes.

# Table of contents
1. Demo-Preview
2. Requirements
3. Design


# Molecular Inversion Probe(MIP)

Molecular Inversion Probes(MIPs) are single-stranded DNA molecules containing two complementary regions that flank the target DNA. These molecules often have a Fluorophore, DNA barcode, or Molecular tag for unique identification.

![MIP_example](https://github.com/SakshiPandey97/MIP_ORACLE/assets/59496870/9d92d545-ffe3-42c6-9125-0c3271ccd35f)


# Rough Design Outline -

1) Start with all possible MIPs by moving along the strand one base pair at a time.
2) Design MIPs for both the forward and reverse strands so that we have the highest probability of binding and then proceed to filter them according to three user-specified criteria: Temperature GC Content Nucleotide Repeats
3) Create a database in MMseqs2 format to only include the human (host) genome. Then filter the MIPs by searching them against the host genome(human).
4) To increase the probability of the MIP binding to the correct target region, search them against the non-redundant nucleotides database.
5) Filter out any MIPs which do not match any other organisms.

# Workflow

# Demo -

1) Obtain sequences of interest in a FASTA format, and make sure the organism name is present in the definition line of each sequence.
2) Following this download all the program files and store them in the same directory as the FASTA file.
3) Fill out the requirements to filter MIPs in the config file provided. The MIPs within the ranges given will be accepted. ex. all MIPs with 45<temp<70 will be taken.


# Run the shell script provided as so:

```bash
MIP_ORACLE.sh -i Trial_File -o trial_final_results -l mip_oracle -j /DATA/databases/blast/nt
```

nohup can also be used:

```bash
nohup MIP_ORACLE.sh -i Trial_File -o trial_final_results -l mip_oracle -j /DATA/databases/blast/nt > trial_log.out &
```

The following files will be generated:
1) The first file will contain all possible MIPs for the sequences provided.
2) The second and third files will contain Passable MIPs(The MIPs that met user requirements as per the config file), and Eliminated MIPs(MIPs that were filtered out).
3) The fourth file is the input file for MMseqs2 search containing arm1+target+arm2 sequences.
4) The fifth and sixth files are the MMseqs2 databases for only human sequences and other non-redundant nucleotides sequence format created from nt BLAST DB input .
5) The seventh file will contain the MIPs with no hits in the MMseqs2 DB, and the eighth file will have the filtered results.
Lastly, the final result file will be generated with filtered MIPs.

Add files image


# Requirements -

Nucleotide BLAST 2.12.0 + with the nt database.

Python 3.6 and the following python packages:

1. pandas=1.1.5
2. biopython=1.70
3. configparser
4. regex
5. xlsxwriter
6. openpyxl

Users can install the required packages through conda using the following command

```bash
conda create -n mip_oracle --file mip_oracle_env.txt
```
