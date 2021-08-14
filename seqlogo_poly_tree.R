#################################################################################
#gene logo :https://cloud.tencent.com/developer/article/1511065
#################################################################################
#install from CRAN
#install.packages("ggseqlogo")
#install from github
#devtools::install.github("omarwagih/ggseqlogo")
#install.packages("spiralize")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require('Biostrings')) BiocManager::install("Biostrings")
if (!require('ggplot2')) install.packages('ggplot2')
if (!require('ggseqlogo')) install.packages('ggseqlogo')
if (!require('seqinr')) install.packages('seqinr')
if (!require('gridExtra')) install.packages('gridExtra')

  
  #read the sequence
  my_fasta <- read.fasta(file = "order_seq_aligned.aln", seqtype = "DNA", as.string = T)
  
  #
  my_fasta_string = vector(mode = 'character')
  for (i in 1:length(my_fasta)){
    my_fasta_string[i] = toupper(c2s(my_fasta[[i]]))
  }
  
  
  seq_len=nchar(my_fasta_string)[1.]
  
  
  
  
  #seqlogo
  #seqlogo <- ggseqlogo(my_fasta_string,seq_type = "dna", method="bits") + theme(axis.text.x = element_blank()) 
  #ggsave("genelogo.pdf",dpi = 300, 
  #       width=300, height=15, unit = "in",limitsize = FALSE)
  #  annotate('rect', xmin = 0.5, xmax = 95.5, ymin = -0.05, ymax = 2.2, alpha = .1, col='black', fill='red')
  
  
  # get sequence based on the start and end site
  my_fasta_string_1_100 <-    subseq(my_fasta_string, start = 1,   end = 100)
  my_fasta_string_101_200 <-  subseq(my_fasta_string, start = 101, end = 200)
  my_fasta_string_201_300 <-  subseq(my_fasta_string, start = 201, end = 300)
  my_fasta_string_301_400 <-  subseq(my_fasta_string, start = 301, end = 400)
  my_fasta_string_401_500 <-  subseq(my_fasta_string, start = 401, end = 500)
  my_fasta_string_501_600 <-  subseq(my_fasta_string, start = 501, end = 600)
  my_fasta_string_601_700 <-  subseq(my_fasta_string, start = 601, end = seq_len)
  

  
  p1 <- ggseqlogo(my_fasta_string_1_100,seq_type = "dna", method="bits") +
    theme(axis.text.x = element_blank()) +
    labs(title = "Base site 1 - 100") 
 
  p2 <- ggseqlogo(my_fasta_string_101_200,seq_type = "dna") + 
    theme(axis.text.x = element_blank()) +
    labs(title = "Base site 101 - 200")

  p3 <- ggseqlogo(my_fasta_string_201_300,seq_type = "dna") + 
    theme(axis.text.x = element_blank())+
    labs(title = "Base site 201 - 300")

  p4 <- ggseqlogo(my_fasta_string_301_400,seq_type = "dna") + 
    theme(axis.text.x = element_blank())+
    labs(title = "Base site 301 - 400")

  p5 <- ggseqlogo(my_fasta_string_401_500,seq_type = "dna") + 
    theme(axis.text.x = element_blank())+
    labs(title = "Base site 401 - 500")

  p6 <- ggseqlogo(my_fasta_string_501_600,seq_type = "dna") + 
    theme(axis.text.x = element_blank())+
    labs(title = "Base site 501 - 600")

  p7 <- ggseqlogo(my_fasta_string_601_700,seq_type = "dna") + 
    theme(axis.text.x = element_blank())+
    labs(title = paste("Base site 601 -",seq_len,sep=" ")) +
    annotate('rect', xmin = seq_len-599, xmax = 100, ymin = -0.05, ymax = 2.2, alpha = .1, col='white', fill='white')
  
  
#  my_fasta_string_701_800 <-  subseq(my_fasta_string, start = 701, end = seq_len)
#  p8 <- ggseqlogo(my_fasta_string_701_800,seq_type = "dna") + 
#    theme(axis.text.x = element_blank())+
#    labs(title = paste("Base site 701 -",seq_len,sep=" "))
  
  
  
  p <- gridExtra::grid.arrange(p1,p2,p3,p4,p5,p6,p7,nrow=7)
  
  ggsave("sequence_logo.pdf",plot = p,dpi = 600, 
       width=10, height=10, unit = "in")
  
  
  
  
  
  
