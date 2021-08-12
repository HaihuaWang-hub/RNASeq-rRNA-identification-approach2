#!/bin/bash -l

module load vsearch


vsearch --cluster_size SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filter_final.fasta \
        -output SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filter_final_uniq.fasta \
        --id 0.99 \
        --centroids \
        --xsize \
        --threads 20 \
        --minseqlength 1 \
        --fasta_width 0
       
        
       
        
