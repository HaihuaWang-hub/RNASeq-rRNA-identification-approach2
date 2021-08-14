if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


if (!require('seqinr')) install.packages('seqinr')
if (!require('ape')) install.packages('ape')


setwd("F:/linux_share/approach_2_test/database")
data1 <- read.fasta(file = "short_database_1.fasta", seqtype = "DNA", as.string = T)
data2 <- read.fasta(file = "short_database_2.fasta", seqtype = "DNA", as.string = T)
data3 <- read.fasta(file = "short_database_3.fasta", seqtype = "DNA", as.string = T)




output_0 <- cbind(data1, data2, data3, fill.with.gaps=TRUE)
output<-output_0[ , !colnames(output_0) %in% c("fill.with.gaps")]
write.dna(output, file="short_database_joint.fasta", 
          format="fasta",
          append = F,
          nbcol = 10,
          colsep = "",
          colw = 1000,
)
