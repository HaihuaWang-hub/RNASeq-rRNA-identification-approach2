#!/bin/bash -l
#SBATCH --job-name=LSU_database_generation_preprocess_vsearch      # Job name
#SBATCH --mail-type=all            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=wanghaihua@ufl.edu       # Where to send mail	
#SBATCH --ntasks=10                      # Number of MPI ranks
#SBATCH --cpus-per-task=9               # Number of cores per MPI rank 
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks-per-node=1             # How many tasks on each node
#SBATCH --ntasks-per-socket=1           # How many tasks on each CPU or socket
#SBATCH --mem-per-cpu=7G             # Memory per core
#SBATCH --time=30-00:00:00                 # Time limit hrs:min:sec
#SBATCH --output=LSU_database_generation_preprocess_vsearch_%j.log     # Standard output and error log

pwd; hostname; date


CPU="18"



module load vsearch


vsearch --derep_fulllength SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filter_final.fasta \
        -output SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filter_final_uniq.fasta \
        --id 0.99999 \
        --centroids \
        --xsize \
        --threads 18 \
        --minseqlength 1 \
        --fasta_width 0
       
        
       
#vsearch --cluster_fast SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filter_final.fasta --id 0.99999 --centroids  SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filter_final_uniq.fasta
