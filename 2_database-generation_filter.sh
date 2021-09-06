#!/bin/bash -l
#SBATCH --job-name=LSU_database_generation_preprocess      # Job name
#SBATCH --mail-type=all            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=wanghaihua@ufl.edu       # Where to send mail
#SBATCH --ntasks=100                      # Number of MPI ranks
#SBATCH --cpus-per-task=9               # Number of cores per MPI rank 
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks-per-node=1             # How many tasks on each node
#SBATCH --ntasks-per-socket=1           # How many tasks on each CPU or socket
#SBATCH --mem-per-cpu=7G             # Memory per core
#SBATCH --time=30-00:00:00                 # Time limit hrs:min:sec
#SBATCH --output=LSU_database_generation_preprocess_%j.log     # Standard output and error log

pwd; hostname; date


CPU="18"


module load metaxa2
module load itsx/1.1b
#module load seqkit
module load clustalo
module load conda
conda activate RNASeq

#This is the code for generation of mulit- short database D1D2 region


###############################################
#pre-filter the LSU region by using the metaxa2

#https://zhuanlan.zhihu.com/p/345533172
#https://microbiology.se/publ/metaxa2_users_guide_2.2.pdf
metaxa2 -i SILVA_138.1_LSUParc_tax_silva_DNA.fasta \
        -o SILVA_138.1_LSUParc_tax_silva_DNA_database_LSU \
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
#Remove the non-LSU region by using the ITSx

#https://github.com/ncbi/ITSx/blob/master/ITSx%20User's%20Guide.pdf
#ITSx -i SILVA_138.1_LSUParc_tax_silva_DNA_database_LSU.eukaryota.fasta \
#     -o SILVA_138.1_LSUParc_tax_silva_DNA_database_LSU_ITSx \
#     -t Fungi \
#     --save_regions SSU,ITS1,5.8S \
#     --only_full F \
#     --selection_priority domains \
#     --complement T \
#     --truncate T \
#     --not_found T \
 #    --preserve T \
 #    --cpu $CPU --multi_thread T

# mv SILVA_138.1_LSUParc_tax_silva_DNA_database_LSU_ITSx_no_detections.fasta SILVA_138.1_LSUParc_tax_silva_LSU_DNA_pre-filtered_database.fasta


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
#filter the LSU region by using one of the D1D2 primers, meanwhile, trim the sequence based on the primers


#generate the sigle-line sequence fasta file
seqkit seq SILVA_138.1_LSUParc_tax_silva_DNA_database_LSU.eukaryota.fasta -w 0 > SILVA_138.1_LSUParc_tax_silva_LSU_DNA_pre-filtered_database_single_line.fasta


fasta_file="SILVA_138.1_LSUParc_tax_silva_LSU_DNA_pre-filtered_database_single_line.fasta"
result_file="SILVA_138.1_LSUParc_tax_silva_LSU_DNA_pre-filtered_database_result.fasta"
forward_LR0R="ACCCGCTGAACTTAAGC"
reverse_LR3_primer="CCGTGTTTCAAGACGGG"
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
    seq_name="$line"
    #echo $line >> $result_file
else
   if [[ "$var2" != "" ]];
   then
       if [[ "$var3" != "" ]];
       then
          echo -ne "$seq_name\n$line\n" |seqkit amplicon -F $forward_LR0R -R $reverse_LR3_primer  >> $result_file
#          echo $line|grep -B 1 "ACCCGCTGAACTTAAGC" |seqkit amplicon -F $forward_LR0R -R $reverse_LR3_primer  >> $result_file
#          echo $line| grep -i -o -P '(?<=ACCCGCTGAACTTAAGC).*(?=CCCGTCTTGAAACACGG)' >> $result_file
       else
          echo -ne "$seq_name\n$line\n" |seqkit amplicon -F $forward_LR0R -r 1:780 >> $result_file
#          echo $line|grep -B 1 "ACCCGCTGAACTTAAGC" |seqkit amplicon -F $forward_LR0R -r 1:780 >> $result_file
#          echo ${line#*${forward_LR0R}}  >> $result_file
       fi
   else
       if [[ "$var3" != "" ]];
       then
          echo -ne "$seq_name\n$line\n" |seqkit amplicon -R $reverse_LR3_primer -r -780:-1 >> $result_file
#          echo $line|grep -B 1 "CCCGTCTTGAAACACGG" |seqkit amplicon -R $reverse_LR3_primer -r -780:-1 >> $result_file
#          echo ${line%${reverse_LR3}*} >> $result_file
       fi
   fi
fi
done < $fasta_file


#remove the accession NO. without sequence and filter the sequence shorter than 250bp
seqkit seq -g -w 0 -m 250 -M 800  $result_file > SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filter_single_line.fasta





awk  -F '.'  '{print $1}' SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filter_single_line.fasta > SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filter_final.fasta






 # clustalo -i SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filter_final.fasta -o SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filter_final.aln




