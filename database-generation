#!/bin/bash -l

module load metaxa2
module load ITSx
module load seqkit

CPU="8"

#This is the code for generation of mulit- short database D1D2 region




###############################################################
#Download the raw LSU database #


#The database was downloaded from SILVA collected un-filtered database 
#As well, the sequences could be collected from International Nucleotide Sequence Database Collaboration (INSDC) databases

workdir="******************"
mkdir $work_dir/raw_database

wget -c https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_LSUParc_tax_silva.fasta.gz $work_dir/raw_database
gunzip $work_dir/raw_database/SILVA_138.1_LSUParc_tax_silva.fasta.gz

#RNA to DNA
seqkit seq --rna2dna 1.fasta > 2.fasta




###############################################
#pre-filter the LSU region by using the metaxa2


#https://microbiology.se/publ/metaxa2_users_guide_2.2.pdf
metaxa2 -i suillu_rhizopogon_database.fasta 
        -o suillu_rhizopogon_database_LSU 
        -f fasta \
        -t e \
        -g lsu \
        --complement T \
        --truncate T \
        --selection_priority domains \
        --not_found T \
        --preserve T \
        --cpu $CPU --multi_thread T 
        
# metaxa2 -i suillu_rhizopogon_database.fasta \     #input file name
#        -o suillu_rhizopogon_database_LSU \        #output
#        -f fasta \                                 #input format {a, auto; f, fasta; q, fastq; p, paired-end; pa, paired-fasta}
#        -t e \                                     #target organisms  {b, bacteria, a, archaea, e, eukaryota,m, mitochondrial, c,chloroplast, A, all}
#        -g lsu \                                   #identify and extract SSU or LSU rRNA genes (lsu/ssu)
#        --mode a \                                 # Controls the Metaxa2 operating mode {m,metagenome, g,genome, a, auto}
#        --complement T \                           # Metaxa2 checks both DNA strands for matches to HMM-profiles.  (T/F)
#        --truncate T \                             # Removes ends of SSU/LSU sequences if they are outside of the SSU/LSU region
#        --selection_priority domains \             # Determines what will be of highest priority when assessing the origin of the sequence {score, sum, domains, eval}
#        --not_found T \                            # If on, Metaxa outputs a list of entries that do notseem to be SSU/LSU sequences
#        --preserve T \
#        --cpu $CPU --multi_thread T 




        
###############################################
#pre-filter the LSU region by using the ITSx

#https://github.com/ncbi/ITSx/blob/master/ITSx%20User's%20Guide.pdf
ITSx -i suillu_rhizopogon_database_LSU.eukaryota.fasta \
     -o suillu_rhizopogon_database_ITS \
     -t Fungi \
     --save_regions SSU,ITS1,5.8S,ITS2,LSU \
     --only_full F \
     --selection_priority domains \
     --complement T \
     --truncate T \
     --not_found T \
     --preserve T \
     --cpu $CPU --multi_thread T




# ITSx -i suillu_rhizopogon_database_LSU.eukaryota.fasta \
#     -o suillu_rhizopogon_database_ITS \
#     -t Fungi \    #Organiams
#     --save_regions SSU,ITS1,5.8S,ITS2,LSU \ 
#     --only_full F \  
#     --selection_priority domains \  # Determines what will be of highest priority when assessing the origin of the sequence {score, sum, domains, eval}
#     --complement T \  checks both DNA strands for matches to HMM-profiles.  (T/F)
#     --truncate T \
#     --not_found T \
#     --preserve T \
#     --cpu $CPU --multi_thread T




###############################################
#filter the LSU region by using one of the D1D2 primers


#generate the sigle-line sequence fasta file
seqkit seq test.fa -w 0


fasta_file="test.fasta"
result_file="result.fasta"
forward_LR0R="ACCCGCTGAACTTAAGC"
reverse_LR3="CCCGTCTTGAAACACGG"
middle_LR21="AAAGGGAAACGCTTGA"
middle_LF402F="CCGATAGCG"
while read line
do
   var1=$(echo $line | grep ">" )
   var2=$(echo $line | grep "${forward_LR0R}" )
   var3=$(echo $line | grep "${reverse_LR3}" )
if [[ "$var1" != "" ]];
then
   echo $line >> $result_file
else
   if [[ "$var2" != "" ]];
   then
       if [[ "$var3" != "" ]];
       then
          echo $line| grep -i -o -P '(?<=ACCCGCTGAACTTAAGC).*(?=CCCGTCTTGAAACACGG)' >> $result_file
       else
          echo ${line#*${forward_LR0R}} >> $result_file
       fi
   else
       if [[ "$var3" != "" ]];
       then
          echo ${line%${reverse_LR3}*} >> $result_file
       fi
   fi
fi
done < test.fasta


