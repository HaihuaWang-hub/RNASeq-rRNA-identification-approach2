#!/bin/bash -l


CPU="8"

#This is the code for generation of mulit- short database D1D2 region


###############################################################
#Download the raw LSU database #


#The database was downloaded from SILVA collected un-filtered database 
#As well, the sequences could be collected from International Nucleotide Sequence Database Collaboration (INSDC) databases

workdir="******************"
mkdir $work_dir/raw_database
cd $work_dir/raw_database

wget -c https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_LSUParc_tax_silva.fasta.gz $work_dir/raw_database
gunzip $work_dir/raw_database/SILVA_138.1_LSUParc_tax_silva.fasta.gz

#RNA to DNA
seqkit seq --rna2dna SILVA_138.1_LSUParc_tax_silva.fasta > SILVA_138.1_LSUParc_tax_silva_DNA.fasta

