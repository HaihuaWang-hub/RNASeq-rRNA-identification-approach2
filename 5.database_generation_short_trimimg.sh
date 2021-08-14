#!/bin/bash -l

module load faSomeRecords
module load clustalo
module load 
module load 

CPU="8"
# conda install -c bioconda entrez-direct
# pip install -U ncbitax2lin /conda install -c bioconda taxonkit 

#This is the code for generation of mulit- short database D1D2 region

mv SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filter_final_uniq_assigned.fasta sequence_sep
mv SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filter_final_uniq_assigned_nmdump_removed.fasta sequence_sep
mv accession2taxid_assigned sequence_sep
mv taxid2lineage_assigned sequence_sep


######################################################################
#trim the seqences in database into 3 short fragments



#1. extract taxonomy names in Order level
grep -i -o -P '(?<=;o__).*(?=;f__)' taxid2lineage_assigned | sort -u > order_level_ids

for Order in $(cat order_level_ids)
do
  touch seq_num_in_level.txt
  printf ${Order} >> seq_num_in_level.txt
  printf "," >> seq_num_in_level.txt
  ###setup the work file named with Order name
  mkdir ${Order}_Order
  ###extract the order level's taxid
  grep -r "${Order}" taxid2lineage_assigned |awk '{print $1}' > ${Order}_Order/${Order}_taxid
  
  ###extract the order level's accession number 
  for taxid in $(cat ${Order}_Order/${Order}_taxid)
  do
     grep -r "${taxid}" accession2taxid_assigned | awk '{print $1}' |sort -u >>  ${Order}_Order/${Order}_accession_num
  done
  
  ###extract the order level's sequences
  for accession in $(cat ${Order}_Order/${Order}_accession_num)
  do
     /home/wang/genome_tools/faSomeRecords SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filter_final_uniq_assigned.fasta \
                                          ${Order}_Order/${Order}_accession_num \
                                          ${Order}_Order/${Order}.fasta
                                          
  done
  
  
  ### count the number of sequences in each order
  grep -c ">" ${Order}_Order/${Order}.fasta >> seq_num_in_level.txt
  
  clustalo -i ${Order}_Order/${Order}.fasta -o ${Order}_Order/${Order}.aln
  
  cp ${Order}_Order/${Order}.aln ${Order}_Order/order_seq_aligned.aln
  
  cd ${Order}_Order
  Rscript 6.seqlogo_plot.R
  mv sequence_logo.pdf ${Order}_seqLogo.pdf
  cd ../
  
  
  
  
  
  
done
   










rm -rf *_Order
















