

#################################################################################
#gene logo :https://cloud.tencent.com/developer/article/1511065
#################################################################################
#直接从CRAN中安装
#install.packages("ggseqlogo")
#从GitHub中安装
#devtools::install.github("omarwagih/ggseqlogo")

#加载包
library(ggplot2)
library(ggseqlogo)
library(seqinr)
library(Biostrings)

setwd("F:/linux_share/approach_2_test/database")

#加载个人数据
my_fasta <- read.fasta(file = "suillus_rhizopogon_database.fas", seqtype = "DNA", as.string = T)

my_protein = read.fasta(file = "proteins.fasta", seqtype = "AA",as.string = TRUE)


#输出fasta文件:https://blog.csdn.net/weixin_42960896/article/details/100339991
x1 = names(my_fasta)
res <- paste0(">",x1,"\n",my_fasta_string_1_100)
res
writeLines(res)

write.table(res,
            file = "test.fasta",
            row.names = F,
            quote = F)


#加载数据
data(ggseqlogo_sample)
seqs_dna
head(seqs_dna)[1]

#将序列转为字符串向量
my_fasta_string = vector(mode = 'character')
for (i in 1:length(my_fasta)){
  my_fasta_string[i] = toupper(c2s(my_fasta[[i]]))
}

# 使用start和end 参数获取子序列
my_fasta_string_1_100 <-  subseq(my_fasta_string, start = 1, end = 100)

#可视化
ggseqlogo(my_fasta_string_1_100)
ggplot()+geom_logo(seqs_dna$MA0001.1)+theme_logo()
ggseqlogo(seqs_dna$MA0001.1)
ggseqlogo(pfms_dna$MA0018.2)

#方法
#ggseqlogo通过method选项支持两种序列标志生成方法：bits和probability。
p1 <- ggseqlogo(seqs_dna$MA0001.1, method="bits")
p2 <- ggseqlogo(seqs_dna$MA0001.1, method="prob")
gridExtra::grid.arrange(p1,p2)



my_fasta_string_1_100 <-  subseq(my_fasta_string, start = 1, end = 100)
p1 <- ggseqlogo(my_fasta_string_1_100,seq_type = "dna")
my_fasta_string_101_200 <-  subseq(my_fasta_string, start = 101, end = 200)
p2 <- ggseqlogo(my_fasta_string_101_200,seq_type = "dna")
my_fasta_string_201_300 <-  subseq(my_fasta_string, start = 201, end = 300)
p3 <- ggseqlogo(my_fasta_string_201_300,seq_type = "dna")
my_fasta_string_301_400 <-  subseq(my_fasta_string, start = 301, end = 400)
p4 <- ggseqlogo(my_fasta_string_301_400,seq_type = "dna")
my_fasta_string_401_500 <-  subseq(my_fasta_string, start = 401, end = 500)
p5 <- ggseqlogo(my_fasta_string_401_500,seq_type = "dna")
my_fasta_string_501_600 <-  subseq(my_fasta_string, start = 501, end = 600)
p6 <- ggseqlogo(my_fasta_string_501_600,seq_type = "dna")
my_fasta_string_601_700 <-  subseq(my_fasta_string, start = 601, end = 700)
p7 <- ggseqlogo(my_fasta_string_601_700,seq_type = "dna")


gridExtra::grid.arrange(p1,p2,p3,p4,p5,p6,p7,nrow=7)

#序列类型
#ggseqlogo支持氨基酸、DNA和RNA序列类型，默认情况下ggseqlogo会自动识别数据提供的序列类型，也可以通过seq_type选项直接指定序列类型。
ggseqlogo(seqs_aa$AKT1, seq_type="aa")


#用数字来代替碱基
seqs_numeric <- chartr("ATGC", "1234", seqs_dna$MA0001.1)
ggseqlogo(seqs_numeric, method="prob", namespace=1:4)


#配色
#ggseqlogo可以使用col_scheme参数来设置配色方案，具体可参考?list_col_schemes
ggseqlogo(seqs_dna$MA0001.1, col_scheme="base_pairing")
#自定义配色
csl <- make_col_scheme(chars = c("A","T", "C", "G"), groups = c("gr1","gr1", "gr2","gr2"), cols = c("purple","purple","blue","blue"))
ggseqlogo(seqs_dna$MA0001.1,col_scheme=csl)
#连续配色
cs2 <- make_col_scheme(chars = c("A", "T", "C", "G"), values = 1:4)
ggseqlogo(seqs_dna$MA0001.1, col_scheme=cs2)


#同时绘制多个序列标志
ggseqlogo(seqs_dna, ncol = 4)







##########################################################
#多基因联合建树---序列串联/连接
##########################################################
library(ape)
library(ggplot2)
library(ggseqlogo)
library(seqinr)
library(Biostrings)
library(tidyr)
rm(list=ls())
setwd("D:/Desktop/bioinfomatic/2_polygenic_tree")


data1 <- read.fasta(file = "suillus_database_seq_from_sunny_short_1.fasta", seqtype = "DNA", as.string = T)
data2 <- read.fasta(file = "suillus_database_seq_from_sunny_short_2.fasta", seqtype = "DNA", as.string = T)
data3 <- read.fasta(file = "suillus_database_seq_from_sunny_short_3.fasta", seqtype = "DNA", as.string = T)




output_0 <- cbind(data1, data2, data3, fill.with.gaps=TRUE)
output<-output_0[ , !colnames(output_0) %in% c("fill.with.gaps")]
write.dna(output, file="suillus_database_seq_from_sunny_short_combined.fasta", 
          format="fasta",
          append = F,
          nbcol = 10,
          colsep = "",
          colw = 1000,
          )


#删除b,d列的处理方式,-which 可以用！代替
dataly[ , -which(colnames(dataly) %in% c("b","d"))]
dataly[ , !colnames(dataly) %in% c("b","d")]
#删除n行的处理方式
dataly[!rownames(dataly) %in% c("n") , ]
















library(ape)
library(phangorn)
library(seqinr)
library(ggtree)
library(ggplot2)
library(patchwork)

rm(list = ls())

setwd('F:/linux_share')

test_dna = read.dna('DOE_EXP3_D1D2_suillus_polygenrtic_tree.fasta', format = 'fasta')
test_phyDat = phyDat(test_dna, type = 'DNA', levels = NULL)

# 模型评估
mt = modelTest(test_phyDat)
print(mt)

# 计算距离
dna_dist = dist.ml(test_phyDat,model = 'JC69')

# NJ树
tree_nj = NJ(dna_dist)
write.tree(tree_nj)

p_nj=ggtree(tree_nj) +
  geom_tiplab(size = 1, color = 'red')+
  labs(title = 'Neighbor Joining Tree')
ggsave("genetree.pdf",dpi = 300, 
       width=8, height=7, unit = "in")
p_nj


# UPGMA树
tree_upgma = upgma(dna_dist)
p_UPGMA = ggtree(tree_upgma) +
  geom_tiplab(size = 1, color = 'black')+
  labs(title = 'UPGMA Tree')
ggsave("genetree.pdf",dpi = 300, 
       width=8, height=7, unit = "in")
p_UPGMA
​
# 最大简约树
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
ggt
# 最大似然法
tree_fit = pml(tree_nj, data = test_phyDat)
print(tree_fit)
fitJC = optim.pml(tree_fit, model = 'JC', rearrangement = 'stochastic')
logLik(fitJC)
tree_bs = bootstrap.pml(fitJC, bs = 1000, optNni = T,
                        multicore = F, 
                        control = pml.control(trace = 0))
tree_ml_bootstrap = plotBS(midpoint(fitJC$tree), tree_bs, p = 50, type = 'p')
​
p_Maximum_Likelihood = ggtree(tree_ml_bootstrap) +
  geom_tiplab(size = 3, color = 'purple')+
  labs(title = 'Maximum Likelihood Tree')
p_Maximum_Likelihood
​
p_all = p_nj + p_UPGMA + p_Maximum_Parsimony + p_Maximum_Likelihood +
  plot_layout(ncol = 2)
ggsave(p_all, filename = 'figures/all.pdf', width = 8, height = 8)




