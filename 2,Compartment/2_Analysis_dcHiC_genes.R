
## ----setup, echo=FALSE, message=FALSE, warning=FALSE----------------------------------
# Set up the environment
library(knitr)
opts_chunk$set(cache.path='cache/', fig.path='img/', cache=F, tidy=T, fig.keep='high', echo=F, dpi=100, warnings=F, message=F, comment=NA, warning=F, results='as.is', fig.width = 10, fig.height = 6) #out.width=700, 
library(pander)
panderOptions('table.split.table', Inf)
set.seed(1)
library(dplyr)
options(stringsAsFactors = FALSE)

#' 
#' # Libraries
#' 
## ----libraries------------------------------------------------------------------------
library(readxl)
library(readr)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(VennDetail)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("VennDetail")

library(annotables)
library("ggsci")
library(scales)
# scales::show_col(pal_lancet("lanonc")(8))
mycols = pal_lancet("lanonc")(8)
library(patchwork)

#' 
#' # Settings
#' 
## -------------------------------------------------------------------------------------
# General settings
# Cutoff for significant AB compartment changes，这些阈值设置都不是很确定！！！！！！！！！！
padj_compartment_cutoff <- 0.3
# Cutoff for significant KEGG enrichment
padj_kegg_cutoff <- 1
# Cutoff for significant MSigDb enrichment，分子特征数据库(MSigDB),这里的阈值设置问题，还是不清楚！！！！！！！！
padj_msigdb_cutoff <- 0.05

## -------------------------------------------------------------------------------------
# dcHiC_2021-09-03 analysis settings
# Mikhail's path 
# dir_data <- "/Users/mdozmorov/Documents/Data/GoogleDrive/HiC_files/results/AB_compartments/dcHiC_2021-09-03"
# fileNameIn1 <- "differential.intra_sample_group.Filtered.bedGraph" # Filtered results
# fileNameIn2 <- "differential.intra_sample_group.bedGraph" # Full results
# dcHiC_2021-12-03 analysis settings
setwd("/mnt/disk4/haitao/bysj_seu/geo_data/hic/script6/dchic/dchic_gene/")
dir_data <- "/mnt/disk4/haitao/bysj_seu/geo_data/hic/script6/HMEC_vs_BT549_100kb_dchic/DifferentialResult"
fileNameIn2 <- "HMEC_BT549_100Kb/viz/files/intra_compartment.bedGraph" # Full results
# Resolution
res_number <- 100000
res_text <- "100kb"
# Results folder
# Mikhail's path 
dir_results <- file.path(dir_data, "results")
# Maggie's path 
dir.create(dir_results, recursive = TRUE) # Create if does not exist
fileNameIn1 <- file.path(dir_results, paste0("AB_gene_summary_", res_text,"_", padj_compartment_cutoff, ".xlsx"))
#这里就是脚本1中选用了各种其他注释数据库的gene list

# Figures,这里从svg修改为pdf文件
fileNameOut1 <- file.path(dir_results, "HMECBT549_AA_RNAseq.pdf")
fileNameOut2 <- file.path(dir_results, "HMECBT549_BB_RNAseq.pdf")
fileNameOut3 <- file.path(dir_results, "HMECBT549_AB_RNAseq.pdf")
fileNameOut4 <- file.path(dir_results, "HMECBT549_BA_RNAseq.pdf")
fileNameOut5 <- file.path(dir_results, "HMECBT549_all_RNAseq.pdf")
fileNameOut6 <- file.path(dir_results, "HMECBT549_ABC_RNAseq.pdf")

#' 
## ----settings-------------------------------------------------------------------------
# Differentially expressed genes
# Single-cell data
# fileNameIn3 <- "/Users/mdozmorov/Documents/Data/GoogleDrive/HiC_files/scRNA-seq/DEGs_0.1_three.xlsx"
# Microarray publication
# fileNameIn3 <- "/Users/mdozmorov/Documents/Data/GoogleDrive/HiC_files/results/RNA-seq/RNA-seq_2019/UCD52_Human_MGT_vs_LiverMet.xlsx"
# Bulk UCD52CR vs PR
fileNameIn3 <- "/mnt/disk4/haitao/bysj_seu/geo_data/hic/script8/RNA/ens_fa_gtf/DEG/DEGs_edgeR_HMEC_BT549_annotated.xlsx"
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!注意这里的DEG是RNA-seq分析中03脚本部分做出来的，而03部分脚本有点难说，之后可以再换换


# Gene annotations基因注释部分与前面的一致
gene_annotations <- grch37[ !(grepl("_", grch37$chr) | grepl("GL", grch37$chr)), c("ensgene", "symbol", "biotype", "description")]
gene_annotations <- gene_annotations[ !duplicated(gene_annotations) & !is.na(gene_annotations$symbol) & gene_annotations$description != "", ]

#' 
#' # Load data
#' 
## ----data-----------------------------------------------------------------------------
# Compartment data
mtx_full <- read_tsv(file.path(dir_data, fileNameIn2)) # Full，还是和1一样使用viz中的intra的bedGraph文件，读取的就是区室的信息部分
# Process compartment data
mtx_filtered <- mtx_full[mtx_full$padj <= padj_compartment_cutoff, ] #筛选对应显著的区室，注意这里是筛选小于这个阈值的
#所以前面有关于区室的分析，只要是筛选区室的都要参考以及检验对应的阈值的信息数据，为什么是0.3
#！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！


mtx_filtered$compartment <- ifelse(mtx_filtered$HMEC > 0 & mtx_filtered$BT549 > 0, "AA", 
                                   ifelse(mtx_filtered$HMEC <= 0 & mtx_filtered$BT549 <= 0, "BB",
                                          ifelse(mtx_filtered$HMEC > 0 & mtx_filtered$BT549 <= 0, "AB", "BA")))




# Genes data，注意下面的数据吗，就是来源于03的RNA-seq的脚本所获取的数据！！！！！！！！！！！！！！！！！注意
#！！！！！！！！！！！！！！！！！！！！！！！！！！！！！所以之后如果RNA-seq那一块的结果修改以及矫正之后，这一块的内容同理操作
genes_AA_full <- read_xlsx(fileNameIn1, sheet = 1)
genes_BB_full <- read_xlsx(fileNameIn1, sheet = 2)
genes_AB_full <- read_xlsx(fileNameIn1, sheet = 3)
genes_BA_full <- read_xlsx(fileNameIn1, sheet = 4)
genes_all_full <- read_xlsx(fileNameIn1, sheet = 5)
# Process gene data
genes_AA <- unique(genes_AA_full$genes)
genes_BB <- unique(genes_BB_full$genes)
genes_AB <- unique(genes_AB_full$genes)
genes_BA <- unique(genes_BA_full$genes)

# DEGs data
degs_full <- read_xlsx(fileNameIn3)
# Filter significant
degs <- degs_full[degs_full$FDR < padj_msigdb_cutoff, ]
# Process DEGs data
degs_UP <- degs$genes[degs$logFC > 0]
degs_DN <- degs$genes[degs$logFC <= 0]
#总之就是两个阈值看看怎么设置的问题！！！！！！！！！！！！！！！！！！！！！！！！！包括RNA-seq那边也是一样，阈值可能会导致最后筛选出来的gene的数目上有些变化

#' 
#' # Genes and AB compartment change stats
#' 
## ----fig.height=2---------------------------------------------------------------------
compartment_gene_summary <- as.data.frame(table(mtx_filtered$compartment))
colnames(compartment_gene_summary) <- c("Type", "Compartment")
compartment_gene_summary$Gene <- c(length(genes_AA), length(genes_AB), length(genes_BA), length(genes_BB))
pander(compartment_gene_summary)
#---------------------------
#Type   Compartment   Gene 
#------ ------------- ------
#  AA       8064       869  

#AB       1431       3984 

#BA       2517       1691 

#BB       7930       408  
#---------------------------

compartment_gene_summary_long <- melt(compartment_gene_summary, id = "Type")
#将宽格式的数据框 compartment_gene_summary 转换为长格式。在长格式中，原始数据框的列名将被整合成两列：一列用于标识不同类别（这里使用了原始数据框中的 "Type" 列），另一列则包含了相应类别的值。长格式通常更适合用于数据的处理和分析。
colnames(compartment_gene_summary_long) <- c("Switch", "Type", "Number")
#Switch        Type Number
#1     AA Compartment   8064
#2     AB Compartment   1431
#3     BA Compartment   2517
#4     BB Compartment   7930
#5     AA        Gene    869
#6     AB        Gene   3984
#7     BA        Gene   1691
#8     BB        Gene    408


# display.brewer.pal(7, "Spectral")
# brewer.pal(7, "Spectral")
# ggplot(compartment_gene_summary_long, aes(x = Type, y = Number, group = Switch)) +
#   geom_bar(stat = "identity", position = "dodge", aes(fill = Switch)) +
#   theme_bw() +
#   scale_fill_manual(values = mycols[1:4])  +
#   facet_wrap(Type ~ .)


p1 <- ggplot(compartment_gene_summary_long %>% filter(Type == "Compartment"), aes(x = Switch, y = Number, group = Switch)) +
  geom_bar(stat = "identity", position = "dodge", aes(fill = Switch)) +
  theme_bw() +
  scale_fill_manual(values = mycols[1:4])  +
  ggtitle("Compartment switch counts")
#筛选出 Type 列为 "Compartment" 的行。然后指定了x轴为 Switch，y轴为 Number，并根据 Switch 列进行分组
#实际上筛选的就是对应类型转换的区室bin（显著）

p2 <- ggplot(compartment_gene_summary_long %>% filter(Type == "Gene"), aes(x = Switch, y = Number, group = Switch)) +
  geom_bar(stat = "identity", position = "dodge", aes(fill = Switch)) +
  theme_bw() +
  scale_fill_manual(values = mycols[c(1, 3, 2, 4)])  +
  ggtitle("Overlapping gene counts")
#同上

p1 + p2
ggsave(paste0("HEMC_BT549_dcHiC_genes_", padj_compartment_cutoff, ".pdf"), width = 6, height = 2)

#' 
#' # Gene overlap
#' 
## -------------------------------------------------------------------------------------
ven <- venndetail(list(AA = genes_AA, BB = genes_BB, AB = genes_AB, BA = genes_BA))
plot(ven, type = "upset")
plot(ven) #这里就是upset集合条形图+venn圆图，这里要清楚一下plot存储，然后才能显示单一的绘图
plot(ven, type = "vennpie") #多层pie图

#UpSet 样式的图表可以显示各个分类之间的交集和差集关系，以及它们的数量
#总的来说就是upset前面4条形都是对应的区室转换的数据，后面4条就是交集
#下面的分析就完全看自己的结果了，根据上面的结果，有交集的gene是
#BB-AB,AA-BA,AA-AB,BB-BA，除此之外其他部分组合就没有了
#factor(ven$Subset)  Levels: AA AA_AB AA_BA AB BA BB BB_AB BB_BA
#因为数据不同，所以分析的时候结果部分和源代码不同，需要进行修改

#' 
#' ## Genes in BB and BA
#' 
## -------------------------------------------------------------------------------------
left_join(data.frame(symbol = getSet(ven, subset = "BB_BA")$Detail %>% sort),
          gene_annotations, by = c("symbol")) %>% dplyr::select(c("symbol", "biotype", "description")) %>% pander


result_BB_BA <- left_join(data.frame(symbol = getSet(ven, subset = "BB_BA")$Detail %>% sort),
          gene_annotations, by = c("symbol")) %>% dplyr::select(c("symbol", "biotype", "description")) 
write.csv(result_BB_BA, "result_BB_BA.csv", row.names = FALSE)
#执行了一个左连接（left join）操作，将两个数据框进行合并。左侧的数据框是由 getSet(ven, subset = "BB_BA")$Detail %>% sort 得到的数据框，其中包含了 Venn 图中 "BB_BA" 部分的基因。右侧的数据框是 gene_annotations，假设包含了基因的注释信息。by = c("symbol") 指定了合并的关键字是基因的符号（symbol）
#使用管道操作符 %>% 将左连接的结果传递给 dplyr::select() 函数，选择感兴趣的列，这里是基因的符号、生物类型（biotype）、描述（description）
#使用 pander() 函数时，它会在控制台中输出漂亮格式的表格
#修改一下以上代码，把这些交集gene保存到指定文件中


#' 
#' ## Genes in BB and AB
#' 
## -------------------------------------------------------------------------------------
result_BB_AB <- left_join(data.frame(symbol = getSet(ven, subset = "BB_AB")$Detail %>% sort),
          gene_annotations, by = c("symbol")) %>% dplyr::select(c("symbol", "biotype", "description")) 
write.csv(result_BB_AB, "result_BB_AB.csv", row.names = FALSE)
#' 
#' ## Genes in AA and AB
#' 
## -------------------------------------------------------------------------------------
result_AA_AB <- left_join(data.frame(symbol = getSet(ven, subset = "AA_AB")$Detail %>% sort),
          gene_annotations, by = c("symbol")) %>% dplyr::select(c("symbol", "biotype", "description")) 
write.csv(result_AA_AB, "result_AA_AB.csv", row.names = FALSE)
#' 
#' ## Genes in AA and BA
#' 
## -------------------------------------------------------------------------------------
result_AA_BA <- left_join(data.frame(symbol = getSet(ven, subset = "AA_BA")$Detail %>% sort),
          gene_annotations, by = c("symbol")) %>% dplyr::select(c("symbol", "biotype", "description")) 
write.csv(result_AA_BA, "result_AA_BA.csv", row.names = FALSE)


#' 
#' # Overlap with DEGs，这里就是联系差异表达gene了
#' 
#' ## Upregulated genes
#' 
## -------------------------------------------------------------------------------------
#还是那句话，阈值设置的不太行，所以对应的数目上太多了，上调以及下调的太多了，当然可以取个交集，3个取交集之后说不定会变小
#总之还是那个RNA-seq脚本中的问题
print(paste("Total upregulated:", length(degs_UP)))  #4102
# Intersect with 
ven <- venndetail(list(AA = genes_AA, BB = genes_BB, AB = genes_AB, BA = genes_BA, UP = degs_UP))
#unique(factor(ven$Subset))
#实际分析过程中可以观察到组合是，一共是17个
#BB_BA_UP AA_BA_UP BA_UP BB_AB_UP AA_AB_UP AB_UP BB_UP AA_UP UP BB_BA
#AA_BA BA BB_AB AA_AB AB BB AA
#这里可以自己抓取数据子集进行分析，比如说如果是观察单一切换与up的关系，就可以获取
#AB BB AA BA UP以及BA_UP  AB_UP  BB_UP  AA_UP
#下面的代码是行不通的，因为这个绘图使用的数据就是venndatail产生的，不可避免的是只要输入数据部分有联系就会检测到
#所有只要有交集，是不可能人工指定数据的
# 选择 Subset 列中包含所需字符串的行
#subset_names <- c("AB", "BB", "AA", "BA", "UP", "BA_UP", "AB_UP", "BB_UP", "AA_UP")
#subset_indices <- grepl(paste(subset_names, collapse = "|"), ven$Subset)
# 从 ven 数据框中提取所需列
#ven_up <- ven[subset_indices, ]
#plot(ven_up,type = "upset")
#plot(getSet(ven, subset =c("AB", "BB", "AA", "BA", "UP", "BA_UP", "AB_UP", "BB_UP", "AA_UP")))
#除了不能直接在绘图中体现以外，其实可以自己提取数据再绘图

plot(ven,type="vennpie") #这里绘制的就是上面的AB分组与up的不同组合，所以有x-up，也有x-y-up
#或者改为绘制"vennpie"
plot(ven,type="upset")


#比较不同组/条件之间基因集合的交集和重叠情况
print(paste("Upregulated, no overlap with AB compartment change:", length(getSet(ven, subset = "UP")$Detail)))
#3007，即上调gene 4102中一共有3007都是和区室转换状态无关的，
#这个结果不太符合我的预期，还是要修改一下数据，无论是dchic还是DEG，最好越往区室联系上靠越好
#学会查看upset图之后，其实就会发现ven中标签up的其实是up全集中没有交集的差集
# Table summary of overlaps
ven_summary <- table(result(ven)$Subset) %>% as.data.frame() 
#将gene数目subset转换为数据框辅助
colnames(ven_summary) <- c("Subset", "Number")
# Compartments overlapping with DEGs，匹配子集中含有 "_UP" 的行
ven_summary[grepl("_UP", ven_summary$Subset), ] %>% pander()

AB_onlyUP <- ven_summary[grepl("_UP", ven_summary$Subset), ]
write.csv(AB_onlyUP, "AB_onlyUP.csv", row.names = FALSE)
#上面这里自己加了2行，直接提取对应up的数据

left_join(getSet(ven, subset = ven_summary$Subset[grepl("_UP", ven_summary$Subset)]), 
          gene_annotations, by = c("Detail" = "symbol")) %>% dplyr::select(c("Subset", "Detail", "biotype", "description")) %>% pander
AB_onlyUP_all <- left_join(getSet(ven, subset = ven_summary$Subset[grepl("_UP", ven_summary$Subset)]), 
                           gene_annotations, by = c("Detail" = "symbol")) %>% dplyr::select(c("Subset", "Detail", "biotype", "description"))
write.csv(AB_onlyUP_all, "AB_onlyUP_all.csv", row.names = FALSE)

#这里保存的是有注释的gene数据，上面是只有gene数据


#' 
#' ## Downregulated genes
#' 
## -------------------------------------------------------------------------------------
print(paste("Total downregulated:", length(degs_DN)))  #4331
# Intersect with 
ven <- venndetail(list(AA = genes_AA, BB = genes_BB, AB = genes_AB, BA = genes_BA, DN = degs_DN))

plot(ven,type="vennpie") 
plot(ven,type="upset")

print(paste("Downregulated, no overlap with AB compartment change:", length(getSet(ven, subset = "DN")$Detail))) #2975
#和上面同理，4331中有2975是与区室无关的

# Table summary of overlaps
ven_summary <- table(result(ven)$Subset) %>% as.data.frame()
colnames(ven_summary) <- c("Subset", "Number")
# Compartments overlapping with DEGs
ven_summary[grepl("_DN", ven_summary$Subset), ] %>% pander()
AB_onlyDN  <- ven_summary[grepl("_DN", ven_summary$Subset), ]
write.csv(AB_onlyDN, "AB_onlyDN.csv", row.names = FALSE)


left_join(getSet(ven, subset = ven_summary$Subset[grepl("_DN", ven_summary$Subset)]), 
          gene_annotations, by = c("Detail" = "symbol")) %>% dplyr::select(c("Subset", "Detail", "biotype", "description")) %>% pander
AB_onlyDN_all <- left_join(getSet(ven, subset = ven_summary$Subset[grepl("_DN", ven_summary$Subset)]), 
                           gene_annotations, by = c("Detail" = "symbol")) %>% dplyr::select(c("Subset", "Detail", "biotype", "description")) 
write.csv(AB_onlyDN_all, "AB_onlyDN_all.csv", row.names = FALSE)



#' 
#' # Correlations
#' 
#' Exploratory analysis. Expectation is that changes in eigenvectors ("D.EV") would correlate with changes in gene expression ("log2FoldChange"). And, it will hold for any compartment change
#' 前面的分析其实就是为什么在AB中会有上调DEG等情况，总体上是无法避免的，所以要用分析手段检测AB区室意义
#' 那就是检测AB与上调，BA与下调的实际相关性
#' ## AA
#' 
## ----fig.height=4, fig.width=5--------------------------------------------------------
DEGs_compartment <- left_join(genes_AA_full[, c("genes", "D.EV")], degs_full[, c("genes", "logFC")], by = c("genes"))
#统一gene，只要信息基因表达的染色质三维结构变化信息（D.EV）和基因表达变化信息（logFC）
DEGs_compartment <- DEGs_compartment[complete.cases(DEGs_compartment), ]  #删除含有缺失值的行
ggplot(DEGs_compartment, aes(y = D.EV, x = logFC, color = logFC)) +
  geom_point(size=2) +
  scale_color_gradient(low = "green", high = "red") +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("dcHiC AA vs. expression changes","r= -0.099 p= 2.581e-02")
#lm 使用线性模型拟合
# Calculate correlation
res <- Hmisc::rcorr(DEGs_compartment$D.EV, DEGs_compartment$logFC)
print(paste("Pearson correlation:", round(res$r[1, 2], digits = 3), "P-value", formatC(res$P[1, 2], format = "e", digits = 3) ))
#"Pearson correlation: -0.099 P-value 2.581e-02"，上面这一句代码的效果要加上去
#可以修改如下，再添加到上面的副标题中，反正是相关系数，和显著p值
print(paste("r=", round(res$r[1, 2], digits = 3), "p=", formatC(res$P[1, 2], format = "e", digits = 3) ))

ggsave(fileNameOut1, width = 7, height = 6)
#效果不是很好，总体是下调的
#！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！修改上游dchic中的gene以及DEG中基因，看看有没有反过来


#' 
#' ## BB
#' 
## ----fig.height=4, fig.width=5--------------------------------------------------------
DEGs_compartment <- left_join(genes_BB_full[, c("genes", "D.EV")], degs_full[, c("genes", "logFC")], by = c("genes"))
DEGs_compartment <- DEGs_compartment[complete.cases(DEGs_compartment), ]
ggplot(DEGs_compartment, aes(y = D.EV, x = logFC, color = logFC)) +
  geom_point(size=2) +
  scale_color_gradient(low = "green", high = "red") +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("dcHiC BB vs. expression changes","r= 0.177 p= 8.693e-03")
# Calculate correlation
res <- Hmisc::rcorr(DEGs_compartment$D.EV, DEGs_compartment$logFC)
print(paste("Pearson correlation:", round(res$r[1, 2], digits = 3), "P-value", formatC(res$P[1, 2], format = "e", digits = 3) ))
print(paste("r=", round(res$r[1, 2], digits = 3), "p=", formatC(res$P[1, 2], format = "e", digits = 3) ))

ggsave(fileNameOut2, width = 7, height = 6)

#' 
#' ## AB
#' 
## ----fig.height=4, fig.width=5--------------------------------------------------------
DEGs_compartment <- left_join(genes_AB_full[, c("genes", "D.EV")], degs_full[, c("genes", "logFC")], by = c("genes"))
DEGs_compartment <- DEGs_compartment[complete.cases(DEGs_compartment), ]
ggplot(DEGs_compartment, aes(y = D.EV, x = logFC, color = logFC)) +
  geom_point(size=2) +
  scale_color_gradient(low = "green", high = "red") +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("dcHiC AB vs. expression changes","r= -0.02 p= 3.150e-01")
# Calculate correlation
res <- Hmisc::rcorr(DEGs_compartment$D.EV, DEGs_compartment$logFC)
print(paste("Pearson correlation:", round(res$r[1, 2], digits = 3), "P-value", formatC(res$P[1, 2], format = "e", digits = 3) ))
print(paste("r=", round(res$r[1, 2], digits = 3), "p=", formatC(res$P[1, 2], format = "e", digits = 3) ))

ggsave(fileNameOut3, width = 7, height = 6)

#' 
#' ## BA
#' 
## ----fig.height=4, fig.width=5--------------------------------------------------------
DEGs_compartment <- left_join(genes_BA_full[, c("genes", "D.EV")], degs_full[, c("genes", "logFC")], by = c("genes"))
DEGs_compartment <- DEGs_compartment[complete.cases(DEGs_compartment), ]
ggplot(DEGs_compartment, aes(y = D.EV, x = logFC, color = logFC)) +
  geom_point(size=2) +
  scale_color_gradient(low = "green", high = "red") +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("dcHiC BA vs. expression changes","r= -0.02 p= 3.150e-01")
# Calculate correlation
res <- Hmisc::rcorr(DEGs_compartment$D.EV, DEGs_compartment$logFC)
print(paste("Pearson correlation:", round(res$r[1, 2], digits = 3), "P-value", formatC(res$P[1, 2], format = "e", digits = 3) ))
print(paste("r=", round(res$r[1, 2], digits = 3), "p=", formatC(res$P[1, 2], format = "e", digits = 3) ))

ggsave(fileNameOut4, width = 7, height = 6)

#' 
#' ## All
#' 
## ----fig.height=4, fig.width=5--------------------------------------------------------
# Collapse all genes by max 
genes_all_full_selected <- genes_all_full[, c("genes", "D.EV")]
maxabs <- function (x) max(abs(x))
genes_all_full_selected <- aggregate(genes_all_full_selected$D.EV, by = list(genes_all_full_selected$genes), "max")
colnames(genes_all_full_selected) <- c("genes", "D.EV") 

DEGs_compartment <- left_join(genes_all_full_selected, degs_full[, c("genes", "logFC")], by = c("genes"))
#DEGs_compartment <- DEGs_compartment[abs(DEGs_compartment$logFC) > 2, ] 
#实际文献绘图过程中是使用了logFC与2边界的数据，所以这句代码理应加上去，但是下面又有操作
DEGs_compartment <- DEGs_compartment[complete.cases(DEGs_compartment), ]
ggplot(DEGs_compartment, aes(y = D.EV, x = logFC, color = logFC)) +
  geom_point(size=2) +
  scale_color_gradient(low = "green", high = "red") +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("dcHiC all vs. expression changes")
# Calculate correlation
res <- Hmisc::rcorr(DEGs_compartment$D.EV, DEGs_compartment$logFC)
print(paste("Pearson correlation:", round(res$r[1, 2], digits = 3), "P-value", formatC(res$P[1, 2], format = "e", digits = 3) ))
print(paste("r=", round(res$r[1, 2], digits = 3), "p=", formatC(res$P[1, 2], format = "e", digits = 3) ))


# Remove unchanged genes with logFC +/- 1SD
dcHiC_D.EV_threshold <- sd(DEGs_compartment$D.EV)
Differential_expression_log2FC_threshold <- sd(DEGs_compartment$logFC)
DEGs_compartment_filtered <- DEGs_compartment %>% dplyr::filter(abs(D.EV) > dcHiC_D.EV_threshold & abs(logFC) > Differential_expression_log2FC_threshold)

ggplot(DEGs_compartment_filtered, aes(y = D.EV, x = logFC, color = logFC)) +
  geom_point(size=2, size = 3) +
  scale_color_gradient(low = "green", high = "red") +
  geom_smooth(method = "lm", se = FALSE) + # , max.overlaps = 10
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("dcHiC all vs. expression changes","r= 0.267 p= 1.158e-11")
ggsave(fileNameOut5, width = 7, height = 6)
# Calculate correlation
res <- Hmisc::rcorr(DEGs_compartment_filtered$D.EV, DEGs_compartment_filtered$logFC)
print(paste("Pearson correlation:", round(res$r[1, 2], digits = 3), "P-value", formatC(res$P[1, 2], format = "e", digits = 3) ))
print(paste("r=", round(res$r[1, 2], digits = 3), "p=", formatC(res$P[1, 2], format = "e", digits = 3) ))

#' 
#' ## ABC transporters only
#' 这部分分析比较奇怪，因为原文献中是用了WGS分析delete的某些区域的gene，某种程度上相当于使用了先验功能的gene（所谓先验，就是gene list不是在hic分析中定位到的）
#' 所以这部分我感觉可以用于CRC分析，分析一下TNBC的CRC是不是不变的
#' 
## -------------------------------------------------------------------------------------
# Collapse all genes by max 
genes_all_full_selected <- genes_all_full[, c("genes", "D.EV")]
maxabs <- function (x) max(abs(x))
genes_all_full_selected <- aggregate(genes_all_full_selected$D.EV, by = list(genes_all_full_selected$genes), "max")
colnames(genes_all_full_selected) <- c("genes", "D.EV") 

DEGs_compartment <- left_join(genes_all_full_selected, degs_full[, c("genes", "logFC")], by = c("genes"))
# DEGs_compartment <- DEGs_compartment[abs(DEGs_compartment$logFC) > 2, ]
DEGs_compartment <- DEGs_compartment[complete.cases(DEGs_compartment), ]
ggplot(DEGs_compartment[grepl("^ABC", DEGs_compartment$genes), ], aes(y = D.EV, x = logFC, color = logFC, label = genes)) +
  geom_point(size=2) +
  geom_text_repel(colour = "black", size = 3) +
  scale_color_gradient(low = "green", high = "red") +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("dcHiC vs. ABC transporter expression")
ggsave(fileNameOut6, width = 3.5, height = 3.5)
# Calculate correlation
res <- Hmisc::rcorr(DEGs_compartment$D.EV[grepl("^ABC", DEGs_compartment$genes)], DEGs_compartment$logFC[grepl("^ABC", DEGs_compartment$genes)])
print(paste("Pearson correlation:", round(res$r[1, 2], digits = 3), "P-value", formatC(res$P[1, 2], format = "e", digits = 3) ))


#' 
#' 
#' # Pathview
#' 
## ----eval = FALSE---------------------------------------------------------------------
## library(pathview)
## OrgDb <- "org.Hs.eg.db"
## species <- "hsa"
## # Differentially expressed genes
## degs <- data.frame(genes = genes_all_full$genes, logFC = genes_all_full$D.EV)
## # Convert to EntrezID
## degs.eg <-clusterProfiler::bitr(degs$genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)
## degs <- left_join(degs, degs.eg, by = c("genes" = "SYMBOL"))
## degs <- degs[!is.na(degs$ENTREZID), ]
## degs <- aggregate(x = degs$logFC, by = list(degs$ENTREZID), FUN = max )
## colnames(degs) <- c("ENTREZID", "logFC")
## # Construct vector of FCs
## degs.genes <- degs$logFC
## names(degs.genes) <- degs$ENTREZID
## 
## # hsa00350	Tyrosine metabolism
## # kegg_ids <- data.frame(ID = c("00350", "05217", "00980", "00982"),
## #                        Description = c("Tyrosine metabolism", "Basal cell carcinoma", "Metabolism of xenobiotics by cytochrome P450", "Drug metabolism - cytochrome P450"))
## kegg_ids <- read_xlsx("/Users/mdozmorov/Documents/Data/GoogleDrive/HiC_files/results/AB_compartments/dcHiC_2021-12-03/results/AB_gene_enrichment_250kb_0.3.xlsx", sheet = "GSEA.KEGG.ALL")
## 
## j <- 0
## # Cycle through each KEGG ID
## for (i in 1:nrow(kegg_ids)) {
##   print(kegg_ids$Description[i])
##   # Get KEGG pathway and overlay DEGs
##   pathview(gene.data = degs.genes, pathway.id = as.character(kegg_ids$ID[i]), species = species, gene.idtype = "ENTREZ", gene.annotpkg = OrgDb, out.suffix = make.names(kegg_ids$Description[i]))
##   # Rename PNG file
##   fileNamePngIn  <- paste0(kegg_ids$ID[i], ".", make.names(kegg_ids$Description[i]), ".png")
##   fileNamePngOut <- paste0(formatC(j, format="g", digits=2, flag = "0"), ".", kegg_ids$ID[i], ".", make.names(kegg_ids$Description[i]), ".png")
##   system(paste0("mv ", fileNamePngIn, " ", fileNamePngOut))
##   j <- j + 1 # Increase counter
##   system(paste0("rm ", kegg_ids$ID[i], ".*")) # Clean up temporary files
## }
## # brew install imagemagick
## system(paste0("convert ", "*", species, "*.png ", "WGS/results/pathways_dcHiC.pdf")) # Combine PNGs into one PDF
## system(paste0("rm ", "*", species, "*.png"))

#' 
#' ## Selected
#' 
## ----eval=FALSE-----------------------------------------------------------------------
## library(pathview)
## library(readxl)
## library(clusterProfiler)
## library(dplyr)
## OrgDb <- "org.Hs.eg.db"
## species <- "hsa"
## 
## dir_data <- "/Users/mdozmorov/Documents/Data/GoogleDrive/HiC_files/results/AB_compartments/dcHiC_2021-12-03/"
## fileNameIn1 <- file.path(dir_data, "results/AB_gene_summary_250kb_1.xlsx")
## # All genes ranked by eigenvector differences
## res <- read_xlsx(fileNameIn1, sheet = "GenesAll")
## degs <- data.frame(genes = res$genes, logFC = res$D.EV)
## # Convert to EntrezID
## degs.eg <-clusterProfiler::bitr(degs$genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)
## degs <- left_join(degs, degs.eg, by = c("genes" = "SYMBOL"))
## degs <- degs[!is.na(degs$ENTREZID), ]
## degs <- aggregate(x = degs$logFC, by = list(degs$ENTREZID), FUN = max )
## colnames(degs) <- c("ENTREZID", "logFC")
## # Construct vector of FCs
## degs.genes <- degs$logFC
## names(degs.genes) <- degs$ENTREZID
## 
## kegg_names = c("ABC transporters",
## "Adherens junction",
## "Breast cancer", # hsa05224
## "Cell cycle",
## "Drug metabolism - cytochrome P450",
## "Estrogen signaling pathway",
## "Focal adhesion",
## "Hedgehog signaling pathway",
## "mTOR signaling pathway",
## "p53 signaling pathway",
## "Pathways in cancer",
## "PI3K-Akt signaling pathway",
## "Ras signaling pathway",
## "Ribosome biogenesis in eukaryotes",
## "Ribosome",
## "TGF-beta signaling pathway",
## "TNF signaling pathway",
## "Toll-like receptor signaling pathway",
## "Wnt signaling pathway",
## "Chemical carcinogenesis",
## "Metabolism of xenobiotics by cytochrome P450",
## "Oxidative phosphorylation")
## 
## xx <- paths.hsa  # Text to ID mapping dataset
## setdiff(kegg_names, xx) # Names that did not map
## # Google and Manually map unmapped IDs
## # xx[grep("Alzheimer", xx, ignore.case = TRUE)]
## kegg_unmapped <- c("05224") # Manually map unmapped ones
## # kegg_unmapped <- NULL # Use if all mapped
## kegg_ids <- c(names(xx)[xx %in% kegg_names]) # All KEGG IDs
## kegg_ids <- kegg_ids[match(kegg_names, xx[xx %in% kegg_names])] # Match name order
## kegg_ids <- sub("hsa", "", kegg_ids) # Strip off "hsa"
## kegg_ids <- kegg_ids[!is.na(kegg_ids)] # Remove NAs
## kegg_ids <- c(kegg_ids, kegg_unmapped) # Attach unmapped
## 
## j <- 0
## # Cycle through each KEGG ID
## for (i in 1:length(kegg_ids)) {
##   # Get KEGG pathway and overlay DEGs
##   pathview(gene.data = degs.genes, pathway.id = as.character(kegg_ids[i]), species = species, gene.idtype = "ENTREZ", gene.annotpkg = OrgDb, out.suffix = "selected")
##   # Rename PNG file
##   fileNamePngIn  <- paste0(species, kegg_ids[i], ".selected.png")
##   fileNamePngOut <- paste0(formatC(j, format="g", digits=2, flag = "0"), ".", species, kegg_ids[i], ".", make.names(kegg_ids[i]), ".png")
##   system(paste0("mv ", fileNamePngIn, " ", fileNamePngOut))
##   j <- j + 1 # Increase counter
##   system(paste0("rm ", species, kegg_ids[i], ".*")) # Clean up temporary files
## }
## # brew install imagemagick
## system(paste0("convert ", "*.", species, "*.png ", file.path(dir_data, "results/pathways_dcHiC_selected.pdf"))) # Combine PNGs into one PDF
## system(paste0("rm ", "*.", species, "*.png"))

#' 
