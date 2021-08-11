#!/bin/bash -l

module load vsearch


vsearch --derep_fulllength SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filter_final.fasta \
        -output SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filter_final_uniq.fasta
        #-sizeout
        
        
        
        
