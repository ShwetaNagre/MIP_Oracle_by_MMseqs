# MIP_ORACLE
A software to filter and identify unique target regions with diagnostic significance in antimicrobial resistance genes, and various other pathogen genomes.

# Table of contents
Demo-Preview
Requirements
Design


Molecular Inversion Probes(MIPs), are single-stranded DNA molecules containing two complementary regions which flank the target DNA. These molecules often have a Fluorophore, DNA barcode, or Molecular tag for unique identification.

# Molecular Inversion Probe(MIP)

# Rough Design Outline -

1) Start with all possible MIPs by moving along the strand one base pair at a time.
2) Design MIPs for both the forward and reverse strands so that we have the highest probability of binding and then proceed to filter them according to three user-specified criteria: Temperature GC Content Nucleotide Repeats
3) Create a database in mmseqs format to only include the human (host) genome. Then filter the MIPs by searching them against the host genome(human).
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
nohup MIP_ORACLE.sh -i Lactobacillus_fermentum_16S -o lacto16S_final_results -l mip_oracle -j /DATA/databases/blast/nt > lacto16S_log.out &
```

The following files will be generated(The first eight files will be in a folder called LOG_FILES):
1) The first file will contain all possible MIPs for the sequences provided.
2) The second and third file will contain Passable MIPs(The MIPs which met user requirements as per the config file), and Eliminated MIPs(MIPs which were filtered out).
3) The fourth file is the BLAST input containing arm1+target+arm2 sequences.
4) The fifth and sixth files are the .xml result files from BLAST.
5) The seventh file will contain the parsed BLAST results about each MIP, and the eighth file will have the filtered results.
Lastly the final result file will be generated in an excel format. image


# Requirements -

Nucleotide BLAST 2.12.0 + with the nt database.

Python 3.6 and the following python packages:

pandas=1.1.5
biopython=1.70
configparser
regex
xlsxwriter
openpyxl

Users can install the required packages through conda using the following command

```bash
conda create -n mip_oracle --file mip_oracle_env.txt
```
