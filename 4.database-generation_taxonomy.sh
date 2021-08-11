#!/bin/bash -l

module load entrez-direct
module load ncbitax2lin
module load taxonkit 
module load 

CPU="8"
# conda install -c bioconda entrez-direct
# pip install -U ncbitax2lin /conda install -c bioconda taxonkit 

#This is the code for generation of mulit- short database D1D2 region


###############################################
#prepare the sequence taxonomy
#link the accession number to taxid and taxname

#1. download database
mkdir taxdump
wget -c https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
mv taxdump.tar.gz taxdump
tar -zxvf taxdump/taxdump.tar.gz -C taxdump

#2. extract the accession number
grep -r ">" SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filter_final_uniq.fasta|cut -d">" -f 2  > accession.number

wc -l accession.number > number_of_seq


#3. assign the accession number to taxid
touch accession2taxid
for i in $(cat accession.number)
do
esummary -db nuccore -id $i | xtract -pattern DocumentSummary -element Caption,TaxId >> accession2taxid
done

# cat accession | epost -db nuccore | esummary -db nuccore | xtract -pattern DocumentSummary -element Caption,TaxId > accession2taxid


#4. extract the taxid
awk '{print $2}' accession2taxid |grep -v -w "0" |sort -u > taxids
wc -l taxids

#5. link taxid to lineage
taxonkit lineage --data-dir taxdump taxids > lineage.txt


#6. formate the lineage
#taxonkit reformat -P --data-dir ../taxdump lineage.txt \
#--delimiter ";" \
#--fill-miss-rank \
#--format "{k};{p};{c};{o};{f};{g};{s}" \
#--miss-rank-repl "0" \
#--miss-rank-repl-prefix "unclassified " \
#| cut -f 1,3 > taxid2lineage

taxonkit reformat -P  lineage.txt \
--data-dir taxdump \
--delimiter ";" \
--format "{k};{p};{c};{o};{f};{g};{s}" \
| cut -f 1,3  > taxid2lineage



#7. refine the lineage (remove the lineage without "Order" level )
 sed "/o__;f__;/d" taxid2lineage > taxid2lineage_assigned

#8. extract assigned taxid
awk '{print $1}' taxid2lineage_assigned > taxid_assigned

#9. extract assigned accession2taxid
touch accession2taxid_assigned

for i in $(cat taxid_assigned)
do
grep -r "$i" accession2taxid >> accession2taxid_assigned
done

awk '{print $1}' accession2taxid_assigned > accession_assigned


wc -l accession.number
wc -l accession_assigned
wc -l taxid_assigned


########################################################################
#extract fasta with new accession ids

/home/wang/genome_tools/faSomeRecords SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filter_final_uniq.fasta \
                                       accession_assigned \
                                       SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filter_final_uniq_assigned.fasta


seqkit seq SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filter_final_uniq_assigned.fasta -w 0 >  SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filtered_assigned.fasta 


mv SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filtered_assigned.fasta 
mv accession2taxid_assigned
mv taxid2lineage_assigned





