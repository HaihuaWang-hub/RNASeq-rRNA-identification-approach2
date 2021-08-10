
#######################################################
##fastqc to show the duplications
#######################################################
fastqc 
multiqc

#######################################################
##fastuniq remove duplicates 
#######################################################

#merge paired end fastq with seqtk
ls *.gz |cut -d"_" -f 1,2,3,4,5,6|sort -u |while read id; do
/home/wang/genome_tools/merge-paired-reads.sh ${id}_R1.fastq.gz ${id}_R2.fastq.gz  ${id}.fastq.gz
done

 ls *.gz |cut -d"_" -f 1,2,3,4,5,6 |while read id; do
 seqkit rmdup --by-seq  --threads 5 ${id}.fastq.gz > dup_removed_data/${id}_dup_removed_fastq
 done


 seqkit rmdup --by-seq --threads 5 AA376R_S97_L003_paired_fungi_rRNA.fastq.gz > AA376R_S97_L003_paired_fungi_rRNA_dup_removed.fastq.gz


ls *.gz |cut -d"_" -f 1,2,3,4,5,6 |while read id; do
mv ${id}_dup_removed_fastq.gz ${id}_dup_removed.fastq
done

ls *.fastq |cut -d"." -f 1 |while read id; do
/home/wang/genome_tools/unmerge-paired-reads.sh ${id}.fastq \
${id}_R1.fastq ${id}_R2.fastq
done


fastuniq -i file_list -t q -o fungi_rRNA_R1.fastq -p fungi_rRNA_R2.fastq -c 0





Trinity --CPU 6 \
    --seqType fq \
    --left  S16_11_clean_pinus_sunny_unaligned_R1.fastq.gz   \
    --right S16_11_clean_pinus_sunny_unaligned_R2.fastq.gz  \
    --max_memory 20G \
    --output Trinity_output &
    
    
    Trinity --CPU 6 \
    --seqType fq \
    --samples_file file_list  \
    --max_memory 20G \
    --output denovo_trinity_output
    
    
    
/home/wang/miniconda3/envs/genomic/opt/trinity-2.12.0/util/misc/get_longest_isoform_seq_per_trinity_gene.pl \
  Trinity.fasta > trinity_unigene.fasta
  
  




#estimate abundance with fungal_D1D2 files
######################################################
mkdir /home/wang/usb/fungal_rRNA_D1D2/alignment_dir/
align_and_estimate_abundance.pl \
--transcripts /home/wang/usb/fungal_rRNA_D1D2/trinity_output/trinity_unigene.fasta \
--seqType fq \
--samples_file file_list \
--est_method RSEM \
--aln_method bowtie2 \
--trinity_mode \
--prep_reference  \
--output_dir /home/wang/usb/fungal_rRNA_D1D2/alignment_dir/ \
--thread_count 6


find * -name '*.isoforms.results'> quant.file 
abundance_estimates_to_matrix.pl --est_method RSEM \
--gene_trans_map /home/wang/usb/fungal_rRNA_D1D2/trinity_output/Trinity.fasta.gene_trans_map \
--name_sample_by_basedir \
--quant_files quant.file








#estimate abundance with bacteria files
######################################################
ls *.gz |cut -d"_" -f 1,2,3,4,5,6,7|sort -u > 1
ls *_R1.fastq.gz > 2
ls *_R2.fastq.gz > 3
paste 1 1 2 3 > file_list

ls *.gz|cut -d"_" -f 1,2,3|sort -u|while read id ;do
if [ -f "/home/wang/usb/2_bacteria+Archae/alignment_dir/${id}/RSEM.isoforms.results.ok" ]; then
   echo "${id} has been analyzed"
   else
align_and_estimate_abundance.pl \
--transcripts /home/wang/usb/2_bacteria+Archae/trinity_genome_guided_output/Trinity-GG_unigene.fasta \
--seqType fq \
--left ${id}_paired_bacteria+Archae_rRNA_R1.fastq.gz \
--right ${id}_paired_bacteria+Archae_rRNA_R2.fastq.gz \
--est_method RSEM \
--aln_method bowtie2 \
--trinity_mode \
--prep_reference  \
--output_dir /home/wang/usb/2_bacteria+Archae/alignment_dir/${id} \
--thread_count 1
fi
done

find * -name '*.isoforms.results'> quant.file 
abundance_estimates_to_matrix.pl --est_method RSEM \
--gene_trans_map /home/wang/usb/bacteria_rRNA/trinity_output/Trinity.fasta.gene_trans_map \
--name_sample_by_basedir \
--quant_files quant.file























grep -r "Boletales" fixrank_fungal_LSU_D1D2.fasta_classified.txt > Boletales.txt
grep -r "Hysterangiales" fixrank_fungal_LSU_D1D2.fasta_classified.txt > Hysterangiales.txt
grep -r "Atractiellales" fixrank_fungal_LSU_D1D2.fasta_classified.txt > Atractiellales.txt
cat Boletales.txt |cut -d ";" -f 1 > Boletales_id.txt
cat Hysterangiales.txt |cut -d ";" -f 1 > Hysterangiales_id.txt
cat Atractiellales.txt |cut -d ";" -f 1 > Atractiellales_id.txt
/home/wang/genome_tools/faSomeRecords trinity_unigene.fasta Boletales_id.txt Boletales_seq.fasta
/home/wang/genome_tools/faSomeRecords trinity_unigene.fasta Hysterangiales_id.txt Hysterangiales_seq.fasta
/home/wang/genome_tools/faSomeRecords trinity_unigene.fasta Atractiellales_id.txt Atractiellales_seq.fasta




#re-assembly
###############################################################
dir="/home/wang/usb/fungal_rRNA_D1D2/first_3_order"
mkdir $dir/database
mkdir $dir/Boletales_assembly_files
mkdir $dir/Atractiellales_assembly_files
mkdir $dir/Hysterangiales_assembly_files
mkdir $dir/Boletales_assemblies
mkdir $dir/Atractiellales_assemblies
mkdir $dir/Hysterangiales_assemblies


bowtie2-build --threads 6  $dir/Hysterangiales_seq.fasta $dir/database/Hysterangiales_seq
bowtie2-build --threads 6  $dir/Boletales_seq.fasta $dir/database/Boletales_seq
bowtie2-build --threads 6  $dir/Atractiellales_seq.fasta $dir/database/Atractiellales_seq


cd /home/wang/usb/fungal_rRNA_D1D2/rawdata
#########Boletales
ls *.gz |cut -d"_" -f 1,2,3 |sort -u |while read id; do
mkdir $dir/Boletales_assembly_files/${id}
bowtie2 -p 5 -x $dir/database/Boletales_seq \
-1 ${id}_fungi_rRNA_R1.fastq.gz \
-2 ${id}_fungi_rRNA_R2.fastq.gz \
-S $dir/Boletales_assembly_files/${id}/${id}.sam \
--al-conc-gz $dir/Boletales_assembly_files/${id}/${id}.fastq.gz \
2>$dir/Boletales_assembly_files/${id}/${id}.bowtie2.log

mv $dir/Boletales_assembly_files/${id}/${id}.fastq.1.gz $dir/Boletales_assembly_files/${id}/${id}_1.fastq.gz
mv $dir/Boletales_assembly_files/${id}/${id}.fastq.2.gz $dir/Boletales_assembly_files/${id}/${id}_2.fastq.gz
rm $dir/Boletales_assembly_files/${id}/${id}.sam

Trinity --CPU 6 \
    --seqType fq \
    --left  $dir/Boletales_assembly_files/${id}/${id}_1.fastq.gz  \
    --right $dir/Boletales_assembly_files/${id}/${id}_2.fastq.gz  \
    --max_memory 15G \
    --output $dir/Boletales_assembly_files/${id}/Trinity_output

/home/wang/miniconda3/envs/genomic/opt/trinity-2.12.0/util/misc/get_longest_isoform_seq_per_trinity_gene.pl \
  $dir/Boletales_assembly_files/${id}/Trinity_output/Trinity.fasta > $dir/Boletales_assemblies/${id}_trinity_unigene.fasta
  
mv $dir/Boletales_assembly_files/${id}/Trinity_output/Trinity.fasta $dir/Boletales_assemblies/${id}_Trinity.fasta

done



######Atractiellales
dir="/home/wang/usb/fungal_rRNA_D1D2/first_3_order"
cd /home/wang/usb/fungal_rRNA_D1D2/rawdata
ls *.gz |cut -d"_" -f 1,2,3 |sort -u |while read id; do
mkdir $dir/Atractiellales_assembly_files/${id}

bowtie2 -p 5 -x $dir/database/Atractiellales_seq \
-1 ${id}_fungi_rRNA_R1.fastq.gz \
-2 ${id}_fungi_rRNA_R2.fastq.gz \
-S $dir/Atractiellales_assembly_files/${id}/${id}.sam \
--al-conc-gz $dir/Atractiellales_assembly_files/${id}/${id}.fastq.gz \
2>$dir/Atractiellales_assembly_files/${id}/${id}.bowtie2.log

mv $dir/Atractiellales_assembly_files/${id}/${id}.fastq.1.gz $dir/Atractiellales_assembly_files/${id}/${id}_1.fastq.gz
mv $dir/Atractiellales_assembly_files/${id}/${id}.fastq.2.gz $dir/Atractiellales_assembly_files/${id}/${id}_2.fastq.gz
rm $dir/Atractiellales_assembly_files/${id}/${id}.sam

Trinity --CPU 3 \
    --seqType fq \
    --left  $dir/Atractiellales_assembly_files/${id}/${id}_1.fastq.gz  \
    --right $dir/Atractiellales_assembly_files/${id}/${id}_2.fastq.gz  \
    --max_memory 15G \
    --output $dir/Atractiellales_assembly_files/${id}/Trinity_output

/home/wang/miniconda3/envs/genomic/opt/trinity-2.12.0/util/misc/get_longest_isoform_seq_per_trinity_gene.pl \
  $dir/Atractiellales_assembly_files/${id}/Trinity_output/Trinity.fasta > $dir/Atractiellales_assemblies/${id}_trinity_unigene.fasta
  
mv $dir/Atractiellales_assembly_files/${id}/Trinity_output/Trinity.fasta $dir/Atractiellales_assemblies/${id}_Trinity.fasta

done






######Hysterangiales
dir="/home/wang/usb/fungal_rRNA_D1D2/first_3_order"
cd /home/wang/usb/fungal_rRNA_D1D2/rawdata
ls *.gz |cut -d"_" -f 1,2,3 |sort -u |while read id; do
mkdir $dir/Hysterangiales_assembly_files/${id}

bowtie2 -p 5 -x $dir/database/Hysterangiales_seq \
-1 ${id}_fungi_rRNA_R1.fastq.gz \
-2 ${id}_fungi_rRNA_R2.fastq.gz \
-S $dir/Hysterangiales_assembly_files/${id}/${id}.sam \
--al-conc-gz $dir/Hysterangiales_assembly_files/${id}/${id}.fastq.gz \
2>$dir/Hysterangiales_assembly_files/${id}/${id}.bowtie2.log

mv $dir/Hysterangiales_assembly_files/${id}/${id}.fastq.1.gz $dir/Hysterangiales_assembly_files/${id}/${id}_1.fastq.gz
mv $dir/Hysterangiales_assembly_files/${id}/${id}.fastq.2.gz $dir/Hysterangiales_assembly_files/${id}/${id}_2.fastq.gz
rm $dir/Hysterangiales_assembly_files/${id}/${id}.sam

Trinity --CPU 6 \
    --seqType fq \
    --left  $dir/Hysterangiales_assembly_files/${id}/${id}_1.fastq.gz  \
    --right $dir/Hysterangiales_assembly_files/${id}/${id}_2.fastq.gz  \
    --max_memory 15G \
    --output $dir/Hysterangiales_assembly_files/${id}/Trinity_output

/home/wang/miniconda3/envs/genomic/opt/trinity-2.12.0/util/misc/get_longest_isoform_seq_per_trinity_gene.pl \
  $dir/Hysterangiales_assembly_files/${id}/Trinity_output/Trinity.fasta > $dir/Hysterangiales_assemblies/${id}_trinity_unigene.fasta
  
mv $dir/Hysterangiales_assembly_files/${id}/Trinity_output/Trinity.fasta $dir/Hysterangiales_assemblies/${id}_Trinity.fasta

done


#################
#rename seq name with file name
##################
ls *_unigene.fasta|while read id;do
awk '/^>/{print ">" substr(FILENAME,1,length(FILENAME)-6); next} 1' ${id} > ${id}.renamed.fasta
done 

cat *.renamed.fasta > @@@_assembly_seqs.fasta

ls *.renamed.fasta|while read id;do
seqkit rename --by-name ${id} > ${id}.num.fasta

grep -r "Suillus;100%" allrank_assembly_seqs.fasta_classified.txt |cut -d ";" -f 1

grep -r "Suillus;100%" allrank_assembly_seqs.fasta_classified.txt |cut -d ";" -f 1 > selected_suillus_gene.txt


/home/wang/genome_tools/faSomeRecords assembly_seqs.fasta selected_suillus_gene.txt selected_suillus_gene.fasta























/home/wang/genome_tools/merge-paired-reads.sh merged_R1.fastq merged_R2.fastq  merged_paried.fastq

 seqkit rmdup --by-seq  --threads 8 merged_paried.fastq > dup_removedmerged_paried.fastq

/home/wang/genome_tools/unmerge-paired-reads.sh dup_removedmerged_paried.fastq dup_removedmerged_paried_R1.fastq dup_removedmerged_paried_R2.fastq

bowtie2-build --threads 6 Step5_D1D2_29342_plus_PMI93_82_plus_Scoth_Tr210_database29539.fasta fungal_rRNA_D1D2

ls *.gz |cut -d"_" -f 1,2,3 |while read id; do
bowtie2 -p 10 -x /home/wang/usb/fungal_rRNA_D1D2/database/fungal_rRNA_D1D2 \
-1 dup_removedmerged_paried_R1.fastq.gz \
-2 dup_removedmerged_paried_R2.fastq.gz \
-S dup_removed_merged.sam 
done


ls *.sam|while read id ;do (samtools sort -o bam -@ 8 -o $(basename ${id} ".sam").bam ${id});done

Trinity --genome_guided_bam dup_removed_merged.bam \
        --max_memory 50G \
        --genome_guided_max_intron 10000 \
        --genome_guided_min_coverage  300 \
        --CPU 1 --output Trinity_output_merged

/home/wang/miniconda3/envs/genomic/opt/trinity-2.12.0/util/misc/get_longest_isoform_seq_per_trinity_gene.pl Trinity-GG.fasta > Trinity-GG_unigene.fasta




mkdir /home/wang/usb/fungal_rRNA_D1D2/alignment_dir_ref_guide_assembly/


align_and_estimate_abundance.pl \
--transcripts /home/wang/usb/fungal_rRNA_D1D2/ref_guide_assembly/Trinity_output_merged/Trinity-GG_unigene.fasta \
--seqType fq \
--samples_file file_list \
--est_method RSEM \
--aln_method bowtie2 \
--trinity_mode \
--prep_reference  \
--output_dir /home/wang/usb/fungal_rRNA_D1D2/alignment_dir_ref_guide_assembly/ \
--thread_count 3





find * -name '*.isoforms.results'> quant.file 

abundance_estimates_to_matrix.pl --est_method RSEM \
--gene_trans_map /home/wang/usb/fungal_rRNA_D1D2/ref_guide_assembly/Trinity_output_merged/Trinity-GG.fasta.gene_trans_map \
--name_sample_by_basedir \
--quant_files quant.file

    
    
    
    
    
##################################################################  
### suillus+rhizopogon separately de novo assembly
##################################################################
bowtie2-build --threads 6 Rhizopogon.fasta rhizopogon
bowtie2-build --threads 6 rhizopogon+suillus.fasta rhizopogon+suillus
bowtie2-build --threads 6 suillus.fasta suillus


rhizopogon="/home/wang/usb/fungal_rRNA_D1D2/rawdata_for_seprate_denovo/ref/rhizopogon"
rhizopogon_suillus="/home/wang/usb/fungal_rRNA_D1D2/rawdata_for_seprate_denovo/ref/rhizopogon+suillus"
suillus="/home/wang/usb/fungal_rRNA_D1D2/rawdata_for_seprate_denovo/ref/suillus"


ls *.gz |cut -d"_" -f 1,2,3,4,5 |sort -u|while read id;do
time bowtie2 -p 2 -x $suillus \
-1 ${id}_R1.fastq.gz \
-2 ${id}_R2.fastq.gz \
-S aligned_seq/${id}.sam \
--al-conc-gz aligned_seq/${id}_aligned.fastq.gz
mv aligned_seq/${id}_aligned.fastq.1.gz aligned_seq/${id}_aligned_R1.fastq.gz
mv aligned_seq/${id}_aligned.fastq.2.gz aligned_seq/${id}_aligned_R2.fastq.gz
rm aligned_seq/${id}.sam
done



ls *.gz |cut -d"_" -f 1,2,3,4,5 |sort -u|while read id;do
time bowtie2 -p 3 -x $rhizopogon_suillus \
-1 ${id}_R1.fastq.gz \
-2 ${id}_R2.fastq.gz \
-S aligned_seq/${id}.sam \
--al-conc-gz aligned_seq/${id}_aligned.fastq.gz
mv aligned_seq/${id}_aligned.fastq.1.gz aligned_seq/${id}_aligned_R1.fastq.gz
mv aligned_seq/${id}_aligned.fastq.2.gz aligned_seq/${id}_aligned_R2.fastq.gz
rm aligned_seq/${id}.sam
done


    
    
    
    
########################################
#D1D2 de novo assembly for each sample
###########################################

ls *.gz|cut -d"_" -f 1,2,3,4,5 |sort -u |while read id;do
trim_galore -q 30 --phred33 --stringency 3 --length 140 \
--paired ${id}_R1.fastq.gz    ${id}_R2.fastq.gz \
--gzip \
--cores 6 \
-o /home/wang/usb/fungal_rRNA_D1D2/rawdata
done



dir=/home/wang/usb/fungal_rRNA_D1D2/cleandata
mkdir $dir/trinity_de_novo_for_each_sample




ls $dir/fastq_data/*.gz |cut -d"/" -f 8 |cut -d"_" -f 1,2,3 |sort -u |while read id;do

if [ -f "$dir/trinity_de_novo_for_each_sample/${id}_Trinity_unigene.fasta" ]; then
   echo "${id} has been analyzed"
   else
time Trinity --CPU 1 --seqType fq \
--left  $dir/fastq_data/${id}_R1.fastq.gz \
--right $dir/fastq_data/${id}_R2.fastq.gz \
--max_memory 2G \
--output $dir/trinity_de_novo_for_each_sample/${id}_trinity_output

/home/wang/miniconda3/envs/genomic/opt/trinity-2.12.0/util/misc/get_longest_isoform_seq_per_trinity_gene.pl $dir/trinity_de_novo_for_each_sample/${id}_trinity_output/Trinity.fasta > $dir/trinity_de_novo_for_each_sample/${id}_Trinity_unigene.fasta
fi
done


##############################
#rename the sequence with sample_name
ls *_unigene.fasta|cut -d"_" -f 1,2,3|while read id;do
mv ${id}_Trinity_unigene.fasta ${id}.fasta
done

ls *.fasta|cut -d "." -f 1| while read id;do
awk '/^>/{print ">" substr(FILENAME,1,length(FILENAME)-6); next} 1' ${id}.fasta > ${id}.renamed.fasta
done 

ls *.renamed.fasta|cut -d"." -f 1|while read id;do
seqkit rename --by-name ${id}.renamed.fasta > ${id}.rename_1.fasta
done

ls *.rename_1.fasta |cut -d"." -f 1 |while read id ; do
awk '{print $1}' ${id}.rename_1.fasta > ${id}.rename_2.fasta
done

cat *.rename_2.fasta > merged.fasta

#### ITSx -t F --cpu 4 --save_regions LSU -i merged_clean_2.fasta -o merged_clean_ITSx  --save_raw merged_clean_ITSx.fasta

/home/wang/genome_tools/faSomeRecords merged.fasta merged.fasta_RDP_classified_Boletales_id boletales_gene.fasta






########################################
#ITS de novo assembly for each sample
###########################################
dir=/home/wang/usb/2_rRNA_fungal_ITS
mkdir $dir/de_novo_for_each_sample
ls $dir/others/*.gz |cut -d"/" -f 7 |cut -d"_" -f 1,2,3,4,5 |sort -u |while read id;do

if [ -f "$dir/de_novo_for_each_sample/${id}_Trinity_unigene.fasta" ]; then
   echo "${id} has been analyzed"
   else
time Trinity --CPU 1 --seqType fq \
--left $dir/others/${id}_R1.fastq.gz \
--right $dir/others/${id}_R2.fastq.gz \
--max_memory 2G \
--output $dir/de_novo_for_each_sample/${id}_trinity_output

/home/wang/miniconda3/envs/genomic/opt/trinity-2.12.0/util/misc/get_longest_isoform_seq_per_trinity_gene.pl $dir/de_novo_for_each_sample/${id}_trinity_output/Trinity.fasta > $dir/de_novo_for_each_sample/${id}_Trinity_unigene.fasta
fi
done






##################################################################################
#approach 2
####################################################################################


#database
#filter the LSU sequence from silva download database 

#https://microbiology.se/publ/metaxa2_users_guide_2.2.pdf
metaxa2 -i suillu_rhizopogon_database.fasta \
        -o suillu_rhizopogon_database_LSU \
        -f fasta \
        -t e \
        -g lsu \
        --complement T \
        --truncate T \
        --selection_priority domains \
        --not_found T \
        --align a \
        --preserve T \
        --cpu 6 --multi_thread T 
        

#https://github.com/ncbi/ITSx/blob/master/ITSx%20User's%20Guide.pdf
ITSx -i suillu_rhizopogon_database_LSU.eukaryota.fasta \
     -o suillu_rhizopogon_database_ITS \
     -t Fungi \
     --save_regions SSU,ITS1,5.8S,ITS2 \
     --only_full F \
     --selection_priority domains \
     --complement T \
     --truncate T \
     --not_found T \
     --preserve T \
     --cpu 8 --multi_thread T




#fasta序列多行变单行
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


seqkit seq -g -w 0 -m 50 $result_file > 123.fasta
clustalo -i 123.fasta -o 123.2.aln
mafft 123.fasta > 123.mafft.aln





















        
clustalo -i suillu_rhizopogon_database_LSU.eukaryota.fasta -o suillu_rhizopogon_database_LSU.eukaryota.aln









#fasta序列多行变单行
seqkit seq test.fa -w 0
#awk '/^>/&&NR>1{print "";}{printf "%s",/^>/?$0"\n":$0}'     test.fa >test2.fa
#awk '{if($0~/>/) name=$0 ;else seq[name]=seq[name]$0;}END{for(i in seq) {if(length(seq[i])>len) print i"\n"seq[i]}}' test.fa >test2.fa

#RNA to DNA
seqkit seq --rna2dna 1.fasta > 2.fasta


##transfer the database to accession number 
awk  -F '.'  '{print $1}' SILVA_138.1_LSUParc_tax_silva_trunc.fasta > SILVA_138.1_LSUParc_tax_silva_trunc_assession.fasta


grep -r ">" SILVA_138.1_LSUParc_tax_silva_trunc_assession.fasta|cut -d">" -f 2  > accession

wc -l accession
##link the accession number to taxid and taxname
conda install -c bioconda entrez-direct
pip install -U ncbitax2lin /conda install -c bioconda taxonkit 
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz


touch accession2taxid
for i in $(cat accession)
do
esummary -db nuccore -id $i | xtract -pattern DocumentSummary -element Caption,TaxId >> accession2taxid
done

cat accession | epost -db nuccore | esummary -db nuccore | xtract -pattern DocumentSummary -element Caption,TaxId > accession2taxid

awk '{print $2}' accession2taxid > taxids_0
grep -v -w "0" taxids_0 > taxids
wc -l taxids_0
wc -l taxids

taxonkit lineage --data-dir ../taxdump taxids > lineage.txt

#taxonkit reformat -P --data-dir ../taxdump lineage.txt \
#--delimiter ";" \
#--fill-miss-rank \
#--format "{k};{p};{c};{o};{f};{g};{s}" \
#--miss-rank-repl "0" \
#--miss-rank-repl-prefix "unclassified " \
#| cut -f 1,3 > taxid2lineage



taxonkit reformat -P --data-dir ../taxdump lineage.txt \
--delimiter ";" \
--format "{k};{p};{c};{o};{f};{g};{s}" \
| cut -f 1,3 > taxid2lineage


 sed "/o__;f__;/d" taxid2lineage > taxid2lineage_final

awk '{print $1}' taxid2lineage_final > taxid_assigned

touch accession2taxid_assigned
for i in $(cat taxid_assigned)
do
grep -w "$i" accession2taxid >> accession2taxid_assigned
done


awk '{print $1}' accession2taxid_assigned > accession_assigned_0
sort -n accession_assigned_0 | uniq > accession_assigned

wc -l accession
wc -l accession_assigned

wc -l taxid_assigned

taxonkit lineage --data-dir ../taxdump taxid_assigned > lineage_assigned.txt

taxonkit reformat -P --data-dir ../taxdump lineage_assigned.txt \
--delimiter ";" \
--fill-miss-rank \
--format "{k};{p};{c};{o};{f};{g};{s}" \
--miss-rank-repl "0" \
--miss-rank-repl-prefix "unclassified " \
| cut -f 1,3 > taxid2lineage_assigned



taxid2lineage_assigned
########################################################################
#extract fasta with new accession ids

/home/wang/genome_tools/faSomeRecords ../SILVA_138.1_LSUParc_tax_silva_trunc_assession.fasta accession_assigned SILVA_138.1_LSUParc_tax_silva_trunc_assigned.fasta

seqkit rmdup SILVA_138.1_LSUParc_tax_silva_trunc_assigned.fasta  -n -o SILVA_138.1_LSUParc_tax_silva_trunc_assigned_dump_removed.fasta

seqkit seq SILVA_138.1_LSUParc_tax_silva_trunc_assigned_dump_removed.fasta -w 0 > SILVA_138.1_LSUParc_tax_silva_trunc_assigned_dump_removed_singleline.fasta



#################################################################################
#deduplication to remove redundant sequences with vsearch (100% similarity)
SILVA_138.1_LSUParc_tax_silva_trunc_assigned_dump_removed_singleline.fasta














#extract genes at family level
grep -A 1 -E "Suillus|Rhizopogon" SILVA_138.1_LSUParc_tax_silva_trunc.fasta > suillu_rhizopogon_database.fasta

seqkit rmdup input.fasta  -s -o output_dump_removed.fasta #去除重复的序列

1. alignment
clustalo -i input.fasta > output.aln


2. trimming


3. short database processing
seqkit seq -g -u input.aln > output.fasaln #remove gap 






fastq_to_fasta -i input.fq -o out.fa -Q 33

ls *.fasta |cut -d"_" -f 1,2,3 |sort -u|while read id ; do
blastn -subject /media/sf_linux_share/approach_2_test/database/suillu_rhizopogon_database_rmdup_2.fasta  \
-query ${id}_R1.fastq.fasta \
-out ${id}_shortdatabase_2.result \
-evalue 1e-5  -num_alignments 3 -num_threads 1 -parse_deflines  -outfmt 6 \
-ungapped -strand both -perc_identity 99 -qcov_hsp_perc 30 
done












blastn -db database/suillus_rhizopogon \
-query AA376R_S97_L003_R1.fastq.fasta \
-out AA376R_S97_L003_R1.fastq.fasta.result \
-evalue 1e-5  -num_alignments 3 -num_threads 1 -parse_deflines \
-ungapped -strand both -perc_identity 99 -qcov_hsp_perc 70




\
-template_length 121



-outfmt 6 \
-max_target_seqs 3
-subject database_name -subject_loc start-stop 






























    
