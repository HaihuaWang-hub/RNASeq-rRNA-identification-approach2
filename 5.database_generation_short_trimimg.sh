#!/bin/bash -l
#SBATCH --job-name=LSU_database_generation_preprocess_trimming      # Job name
#SBATCH --mail-type=all            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=wanghaihua@ufl.edu       # Where to send mail	
#SBATCH --ntasks=1                      # Number of MPI ranks
#SBATCH --cpus-per-task=1               # Number of cores per MPI rank 
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks-per-node=1             # How many tasks on each node
#SBATCH --ntasks-per-socket=1           # How many tasks on each CPU or socket
#SBATCH --mem-per-cpu=7G             # Memory per core
#SBATCH --time=30-00:00:00                 # Time limit hrs:min:sec
#SBATCH --output=LSU_database_generation_preprocess_trimming_%j.log     # Standard output and error log

pwd; hostname; date


CPU="18"


#module load faSomeRecords
module load ucsc/20210803
module load clustalo
module load rstudio/1.4
module load seqkit/0.10.2
#module load mega/7.0.26


#This is the code for generation of mulit- short database D1D2 region
work_dir="/home/wanghaihua/all_data/database_generation/fungal_lsu_database/result_dir"
cd $work_dir
# SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filter_final_uniq_assigned.fasta 
# SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filter_final_uniq_assigned_nmdump_removed.fasta 
# accession2taxid_fungi_assigned 
# taxid2lineage_fungi 


######################################################################
#trim the seqences in database into 3 short fragments

#1. extract taxonomy names in Order level
grep -r ";o__" taxid2lineage_fungi   > taxid2lineage_fungi_0

while read line
do
echo $line  |cut -d ";" -f 4 |cut -d "_" -f 3  >> order_level_ids_0
done < taxid2lineage_fungi_0

#count the non-order taxonomy number
cat order_level_ids_0 | grep -v ^$ | wc -l

cat order_level_ids_0 |sort -u > order_level_ids
rm order_level_ids_0
rm taxid2lineage_fungi_0

mkdir $work_dir/sep_process
mkdir $work_dir/short_seqs
touch seq_num_in_level.txt

for Order in $(cat order_level_ids)
do

  printf "\n"${Order} >> seq_num_in_level.txt
  printf "," >> seq_num_in_level.txt
  ###setup the work file named with Order name
  mkdir $work_dir/sep_process/${Order}_Order
  
  ###extract the order level's taxid
  grep -r "${Order}" taxid2lineage_fungi |awk '{print $1}' > $work_dir/sep_process/${Order}_Order/${Order}_taxid
  
  ###extract the order level's accession number 
  for taxid in $(cat $work_dir/sep_process/${Order}_Order/${Order}_taxid)
  do
     grep -w "${taxid}" accession2taxid_fungi_assigned | awk '{print $1}'|cut -d ":" -f 2 |sort -u >>  $work_dir/sep_process/${Order}_Order/${Order}_accession_num
  done
  
  ###extract the order level's sequences
     faSomeRecords SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filter_final_uniq_assigned.fasta \
                                          $work_dir/sep_process/${Order}_Order/${Order}_accession_num \
                                          $work_dir/sep_process/${Order}_Order/${Order}.fasta

  
  ### count the number of sequences in each order
  printf $(cat $work_dir/sep_process/${Order}_Order/${Order}.fasta |grep -c ">" ) >> seq_num_in_level.txt
  
  clustalo --force  -i $work_dir/sep_process/${Order}_Order/${Order}.fasta -o $work_dir/sep_process/${Order}_Order/${Order}.aln
  
  cp $work_dir/sep_process/${Order}_Order/${Order}.aln $work_dir/sep_process/${Order}_Order/order_seq_aligned.aln
  
  cd $work_dir/sep_process/${Order}_Order
  Rscript  /home/wanghaihua/scripts_file/LSU_database_generation/6.seqlogo_plot.R
  mv sequence_logo.pdf ${Order}_seqLogo.pdf
  
  Rscript /home/wanghaihua/scripts_file/LSU_database_generation/7_sliding_window_comparation.r
  
  cat cut_start_1.txt |awk 'NR==2{print $2}' > cut_start_1
  cat cut_end_1.txt |awk 'NR==2{print $2}' > cut_end_1
  cat cut_start_2.txt |awk 'NR==2{print $2}' > cut_start_2
  cat cut_end_2.txt |awk 'NR==2{print $2}' > cut_end_2
  cat cut_start_3.txt |awk 'NR==2{print $2}' > cut_start_3
  cat cut_end_3.txt |awk 'NR==2{print $2}' > cut_end_3
  rm cut_start_1.txt cut_end_1.txt  cut_start_2.txt cut_end_2.txt cut_start_3.txt cut_end_3.txt
  
  cat order_seq_aligned.aln | seqkit subseq -r $(cat cut_start_1):$(cat cut_end_1) > database_1.aln
  cat order_seq_aligned.aln | seqkit subseq -r $(cat cut_start_2):$(cat cut_end_2) > database_2.aln
   cat order_seq_aligned.aln | seqkit subseq -r $(cat cut_start_3):$(cat cut_end_3) > database_3.aln
  cat database_1.aln | seqkit seq -g -l > ${Order}_short_1.fasta
  cat database_2.aln | seqkit seq -g -l > ${Order}_short_2.fasta
  cat database_3.aln | seqkit seq -g -l > ${Order}_short_3.fasta
  
  #Rscript ../8.shortdatabase_joint.R
  cp ${Order}_short_1.fasta $work_dir/short_seqs
  cp ${Order}_short_2.fasta $work_dir/short_seqs
   cp ${Order}_short_3.fasta $work_dir/short_seqs
  cd $work_dir
  
done
   










#rm -rf *_Order
















