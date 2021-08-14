
#################################################################################
#gene logo :https://cloud.tencent.com/developer/article/1511065
#################################################################################
#install from CRAN
#install.packages("ggseqlogo")
#install from github
#devtools::install.github("omarwagih/ggseqlogo")
#install.packages("spiralize")


library(ggplot2)
library(ggseqlogo)
library(seqinr)
library(Biostrings)


#read the sequence
my_fasta <- read.fasta(file = "Boletales.aln", seqtype = "DNA", as.string = T)


#将序列转为字符串向量
my_fasta_string = vector(mode = 'character')
for (i in 1:length(my_fasta)){
  my_fasta_string[i] = toupper(c2s(my_fasta[[i]]))
}


seq_len=nchar(my_fasta_string)[1.]

#seqlogo
ggseqlogo(my_fasta_string_1_100,seq_type = "dna", method="bits") + theme(axis.text.x = element_blank()) +
          annotate('rect', xmin = 0.5, xmax = 95.5, ymin = -0.05, ymax = 2.2, alpha = .1, col='black', fill='red')

# get sequence based on the start and end site
my_fasta_string_1_100 <-    subseq(my_fasta_string, start = 1,   end = 100)
p1 <- ggseqlogo(my_fasta_string_1_100,seq_type = "dna", method="bits") +
                 theme(axis.text.x = element_blank()) +
                 labs(title = "Base site 1 - 100")
my_fasta_string_101_200 <-  subseq(my_fasta_string, start = 101, end = 200)
p2 <- ggseqlogo(my_fasta_string_101_200,seq_type = "dna") + 
                theme(axis.text.x = element_blank()) +
                labs(title = "Base site 101 - 200")
my_fasta_string_201_300 <-  subseq(my_fasta_string, start = 201, end = 300)
p3 <- ggseqlogo(my_fasta_string_201_300,seq_type = "dna") + 
                theme(axis.text.x = element_blank())+
                labs(title = "Base site 201 - 300")
my_fasta_string_301_400 <-  subseq(my_fasta_string, start = 301, end = 400)
p4 <- ggseqlogo(my_fasta_string_301_400,seq_type = "dna") + 
                theme(axis.text.x = element_blank())+
                labs(title = "Base site 301 - 400")
my_fasta_string_401_500 <-  subseq(my_fasta_string, start = 401, end = 500)
p5 <- ggseqlogo(my_fasta_string_401_500,seq_type = "dna") + 
                theme(axis.text.x = element_blank())+
                labs(title = "Base site 401 - 500")
my_fasta_string_501_600 <-  subseq(my_fasta_string, start = 501, end = 600)
p6 <- ggseqlogo(my_fasta_string_501_600,seq_type = "dna") + 
                theme(axis.text.x = element_blank())+
                labs(title = "Base site 501 - 600")
my_fasta_string_601_700 <-  subseq(my_fasta_string, start = 601, end = seq_len)
p7 <- ggseqlogo(my_fasta_string_601_700,seq_type = "dna") + 
                theme(axis.text.x = element_blank())+
                labs(title = paste("Base site 601 -",seq_len,sep=" "))
my_fasta_string_701_800 <-  subseq(my_fasta_string, start = 701, end = seq_len)
p8 <- ggseqlogo(my_fasta_string_701_800,seq_type = "dna") + 
                theme(axis.text.x = element_blank())+
                 labs(title = paste("Base site 701 -",seq_len,sep=" "))



gridExtra::grid.arrange(p1,p2,p3,p4,p5,p6,p7,nrow=7)



rm(list=ls())

##########################################################
#Multi-gene polygenetic tree
##########################################################
library(ape)
library(ggplot2)
library(ggseqlogo)
library(seqinr)
library(Biostrings)
library(tidyr)


data1 <- read.fasta(file = "suillus_database_seq_from_sunny_short_1.fasta", seqtype = "DNA", as.string = T)
data2 <- read.fasta(file = "suillus_database_seq_from_sunny_short_2.fasta", seqtype = "DNA", as.string = T)
data3 <- read.fasta(file = "suillus_database_seq_from_sunny_short_3.fasta", seqtype = "DNA", as.string = T)




output_0 <- cbind(data1, data2, data3, fill.with.gaps=TRUE)
output<-output_0[ , !colnames(output_0) %in% c("fill.with.gaps")]
write.dna(output, file="suillus_database_seq_from_sunny_short_joint.fasta", 
          format="fasta",
          append = F,
          nbcol = 10,
          colsep = "",
          colw = 1000,
)








if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")


if (!requireNamespace(c('ape','phangorn','seqinr'))) {
  install.packages(c('ape','phangorn','seqinr'))
}



library(ape)
library(phangorn)
library(seqinr)
library(ggtree)
library(ggplot2)
library(patchwork)


test_dna = read.dna('Boletales.aln', format = 'fasta')
test_phyDat = phyDat(test_dna, type = 'DNA', levels = NULL)

# model evaluation
mt = modelTest(test_phyDat)
print(mt)

# calculate the distance
dna_dist = dist.ml(test_phyDat,model = 'JC69')

# NJ tree
tree_nj = NJ(dna_dist)
write.tree(tree_nj)

p_nj=ggtree(tree_nj) +
  geom_tiplab(size = 1, color = 'red')+
  labs(title = 'Neighbor Joining Tree')
ggsave("genetree.pdf",dpi = 300, 
       width=8, height=7, unit = "in")
p_nj


# UPGMA tree
tree_upgma = upgma(dna_dist)
p_UPGMA = ggtree(tree_upgma) +
  geom_tiplab(size = 1, color = 'black')+
  labs(title = 'UPGMA Tree')
ggsave("genetree.pdf",dpi = 300, 
       width=8, height=7, unit = "in")
p_UPGMA
​
# Maximal parsimony tree
parsimony(tree_upgma, data = test_phyDat)
parsimony(tree_nj, data = test_phyDat)
test_optim = optim.parsimony(tree_nj, data = test_phyDat)
tree_peatchet = pratchet(test_phyDat)

p_Maximum_Parsimony = ggtree(tree_peatchet) +
  geom_tiplab(size = 0.6, color = 'blue')+
  labs(title = 'Maximum Parsimony Tree')
p_Maximum_Parsimony
ggsave("genetree.pdf",dpi = 300, 
       width=8, height=7, unit = "in")

# Maximum Likelihood Method
tree_fit = pml(tree_nj, data = test_phyDat)
print(tree_fit)
fitJC = optim.pml(tree_fit, model = 'JC', rearrangement = 'stochastic')
logLik(fitJC)
tree_bs = bootstrap.pml(fitJC, bs = 1000, optNni = T,
                        multicore = F, 
                        control = pml.control(trace = 0))
tree_ml_bootstrap = plotBS(midpoint(fitJC$tree), tree_bs, p = 50, type = 'p')

p_Maximum_Likelihood = ggtree(tree_ml_bootstrap) +
  geom_tiplab(size = 3, color = 'purple')+
  labs(title = 'Maximum Likelihood Tree')
p_Maximum_Likelihood


#joint all method
p_all = p_nj + p_UPGMA + p_Maximum_Parsimony + p_Maximum_Likelihood +
  plot_layout(ncol = 2)
ggsave(p_all, filename = 'figures/all.pdf', width = 8, height = 8)


























