
## ----setup, echo=FALSE, message=FALSE, warning=FALSE----------------------------------
# Set up the environment
library(knitr)
opts_chunk$set(cache.path = "cache/", fig.path = "img/", cache = F, tidy = T, fig.keep = "high", echo = F, dpi = 100, warnings = F, message = F, comment = NA, warning = F, results = "as.is", fig.width = 10, fig.height = 6, cache.lazy = FALSE) # out.width=700,
library(pander)
panderOptions("table.split.table", Inf)
set.seed(1)
library(dplyr)
options(stringsAsFactors = FALSE)

#' 
#' Processing differential compartment results from `dcHi-C` analysis. Significant changes may be from A to A (AA), B to B (BB), A to B (AB), and B to A (BA).
#' 
#' - Proportion of significant AB compartment changes, genome- and Chromosome-specific. Sortable tables and plots.
#' - Karyoplots
#' - Overlapping gene enrichment analysis, KEGG, MSigDb.
#' 
#' - BED files for each type of AB compartment change. Named as `BED_<compartment>_<resolution>_<padj_compartment_cutoff>.bed`, e.g., `BED_AA_250kb_0.01.xlsx`
#' - `AB_gene_summary_250kb.xlsx` Excel file with genes overlapping AB compartment switches (individual worksheets), other summaries.

#' # Libraries

## ----libraries------------------------------------------------------------------------
library(annotables)
library(rCGH)
library(GenomicRanges)
library(clusterProfiler)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
#这里安装包的时候遇到了一些非0退出status的问题，可以通过手动安装一些R包来解决问题
#https://github.com/YuLab-SMU/clusterProfiler/issues/665 主要是手动安装前面两个包，安装了下面这两个包之后再去安装上面的clusterprofiler
#但是后面的BSgenome会需要用该包的更高级，所以还是在下面下载的时候全部更新
download.file("https://cran.r-project.org/src/contrib/Archive/igraph/igraph_1.6.0.tar.gz", "igraph_1.6.0.tar.gz")
install.packages("igraph_1.6.0.tar.gz", repos = NULL, type = "source")
download.file("https://cran.r-project.org/src/contrib/Archive/tidygraph/tidygraph_1.3.0.tar.gz", "tidygraph_1.3.0.tar.gz")
install.packages("tidygraph_1.3.0.tar.gz", repos = NULL, type = "source")


# library(enrichR)
library(tidyr)
library(stringr)
library(writexl)
library(readr)
library(ggplot2)
library(reshape2)
library(ggsci)
library(grid)
library(gridExtra)
library(BSgenome)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("BSgenome")


library(msigdbr)  #检索基因集及其成员基因的数据框架
library(karyoploteR)
BiocManager::install("karyoploteR")
#事实上只要安装了一次biocmanager，后续就不要再执行上面那条代码了

# options(scipen = 999999999999999) 后面会有一系列的NAs introduced by coercion to integer range，不知道是不是这个，暂时不用
library(data.table)
library(MDmisc) # BiocManager::install("mdozmorov/MDmisc", update = FALSE)
#install_github('mdozmorov/MDmisc')
#remotes::install_github("mdozmorov/MDmisc")
#BiocManager::install("mdozmorov/MDmisc")
devtools::install_github('mdozmorov/MDmisc')

#！！！！！！！！！！！！！目前MDmisc这个包暂时安装不上，看到时候在其他的项目中使用能否安装，现在能够装了，就是在前面的9999999那一段要禁止掉重新运行




#' # Settings
# General settings
# Cutoff for significant AB compartment changes，此处的显著值看看在哪里用得到
padj_compartment_cutoff <- 0.3
# Rerun setting, affects overwriting files.
rerun <- TRUE
# How to sort the barchart. If TRUE, only a sum of AB and BA changes is used. 
# If FALSE, the total significant changes (AB, BA, AA, BB) is used，效果目前未知
AB_BA <- TRUE


## ## Dozmorov's Path
## # dir_data <- "/Users/mdozmorov/Documents/Data/GoogleDrive/HiC_files/04-28-2021_dcHiC-results-with-problematic-chrs/dchic_results_2pc/DifferentialCompartment"
## ## Maggie's path
## # dir_data <- "~/Google Drive/My Drive/HiC_files/04-28-2021_dcHiC-results-with-problematic-chrs/dchic_results_2pc/DifferentialCompartment"
## # fileNameIn1 <- "UCD52PR_vs_UCD52CR_differential_compartments.bedGraph" # Filtered results
## # fileNameIn2 <- "UCD52PR_vs_UCD52CR_full_compartment_details.bedGraph" # Full results
#上面的输出结果无法和dchic的输出对应上！！！！！！！！！！！！！！！！！！！！！！！！！！！！！

## # # Resolution
## # res_number <- 100000
## # res_text <- "100kb"
## # Results folder
## # Mikhail's path
## # dir_results <- file.path("/Users/mdozmorov/Documents/Data/GoogleDrive/HiC_files/results/04-28-2021_dcHiC-results-with-problematic-chrs", res_text)
## # Maggie's path
## # dir_results <- file.path("~/Google Drive/My Drive/HiC_files/results/04-28-2021_dcHiC-results-with-problematic-chrs", res_text)
## # dir.create(dir_results, recursive = TRUE) # Create if does not exist
## # fileNameOut1 <- file.path(dir_results, paste0("AB_gene_summary_", res_text, ".xlsx"))

#' 
## -------------------------------------------------------------------------------------
# dcHiC_2021-09-03 analysis settings
# Mikhail's path 
# dir_data <- "/Users/mdozmorov/Documents/Data/GoogleDrive/HiC_files/results/AB_compartments/dcHiC_2021-09-03"
# fileNameIn1 <- "differential.intra_sample_group.Filtered.bedGraph" # Filtered results
# fileNameIn2 <- "differential.intra_sample_group.bedGraph" # Full results
#这两个文件又是对应的什么！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！


dir_data <- "/mnt/disk4/haitao/bysj_seu/geo_data/hic/script6/HMEC_vs_BT549_100kb_dchic/DifferentialResult" #修改成自己dchic分析出来的结果目录
#/Users/mdozmorov/Documents/Data/GoogleDrive/HiC_files/results/AB_compartments/dcHiC_2021-09-03类比
fileNameIn2 <- "HMEC_BT549_100Kb/viz/files/intra_compartment.bedGraph" # Full results
#CR_vs_PR_After_PC_Reselction_PR_chr3_PC1/viz/files/intra_compartment.bedGraph类比，可以看出一下类似关系：
#HMEC_BT549_100Kb文件夹目录=CR_vs_PR_After_PC_Reselction_PR_chr3_PC1/，因为下一级都有viz
#DifferentialResult文件夹目录=/Users/mdozmorov/Documents/Data/GoogleDrive/HiC_files/results/AB_compartments/dcHiC_2021-09-03/差不多，因为是file的上一级文件
#当然其实可以使用绝对路径

# Resolution
res_number <- 100000
res_text <- "100kb"
# Results folder 存储结果的文件夹
dir_results <- file.path(dir_data, "results")
dir.create(dir_results, recursive = TRUE) # Create if does not exist
# All results
fileNameOut1 <- file.path(dir_results, paste0("AB_gene_summary_", res_text,"_", padj_compartment_cutoff, ".xlsx"))
# BED file name
fileNameOut2 <- file.path(dir_results, paste0("AB_", res_text, "_", padj_compartment_cutoff, ".bed"))
#参照第19，20行，主要是这个阈值看看后期要不要改


#' ## Organism selection
## ----organism-------------------------------------------------------------------------
chr <- paste0("chr", c(1:22)) # Chromosomes，这里其实是可以改成知道后续的chrX染色体的，但是因为实际上自己在dchic中已经移除了chrY,所以不确定这里是否会影响
#c(paste("chr", 1:22, sep = ""), "chrX", "chrY")，先看看chr1-22效果是什么
#这里估计就是之后要研究的对象list

# Import centromeric regions from rCGH 着丝粒中心粒区域
hg19_centro <- hg19
# Change chromosome indicator to include "chr"
hg19_centro$chrom <- paste0("chr", hg19_centro$chrom)  #这里其实可以改一下，就是chr23以及24改成X和Y
# Make GRanges from centromeric locations
hg19_centro.gr <- makeGRangesFromDataFrame(hg19_centro, seqnames.field = "chrom", start.field = "centromerStart", end.field = "centromerEnd", keep.extra.columns = TRUE)
#新的 GRanges 对象，它是通过调用 makeGRangesFromDataFrame 函数创建的。这个对象包含了从 hg19_centro 数据框中提取的染色体、起始位置和终止位置信息


library(org.Hs.eg.db)
OrgDb <- "org.Hs.eg.db"
species <- "hsa"
KEGG <- "KEGG_2019_Human"
#指定要使用的数据库、物种和 KEGG 数据集

# Annotables
gene_annotations <- grch37[!(grepl("_", grch37$chr) | grepl("GL", grch37$chr)), c("symbol", "description", "biotype")]
#从 grch37 数据框中选择了三列数据，分别是 "symbol"、"description" 和 "biotype"。grch37 数据框中的 "chr" 列包含了染色体信息，这行代码使用 grepl 函数排除了包含下划线 "_" 或 "GL" 的染色体。c 函数用于选择列
#应该就是获取规则染色体的这三列注释，可以通过grch37查看
gene_annotations <- gene_annotations[!duplicated(gene_annotations$symbol) & !is.na(gene_annotations$symbol) & gene_annotations$description != "", ]
#排除了重复的基因符号、缺失的基因符号和空白的基因描述

# BSgenome settings
#library(BiocManager)
#install("BSgenome.Hsapiens.UCSC.hg19")
bsgenome <- "BSgenome.Hsapiens.UCSC.hg19" #这里去不去掉引号都是不影响的
chrom.sizes <- data.frame(chr = chr, size = seqlengths(getBSgenome(genome = bsgenome, masked = FALSE))[chr])
#这里用上了chr，对应的是研究对象chr1-22，后续可以改成chrX看一下效果，参考第139行
#设置 BSgenome 对象的参数，并创建了一个数据框，其中包含了染色体的信息，获得了指定chr+长度信息

# Get all human genes
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
genomewide.genes <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
#single.strand.genes.only=FALSE,上面这里会有一些正义反义链上的gene丢失掉
genomewide.genes <- keepSeqlevels(genomewide.genes, c(paste0("chr", 1:22), "chrX"), pruning.mode = "tidy")
#调用了 keepSeqlevels 函数，用于保留指定的染色体。在这个例子中，保留了 1 到 22 号染色体和 X 染色体
gene_symbol <- bitr(genomewide.genes$gene_id, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = OrgDb)
#调用了 bitr 函数，用于将基因 ID 转换为基因符号。fromType 参数指定了输入的 ID 类型，toType 参数指定了输出的 ID 类型，OrgDb 参数指定了要使用的基因注释数据库
gene_symbol <- left_join(gene_symbol, gene_annotations, by = c("SYMBOL" = "symbol"))
#将基因符号与基因注释数据进行关联。by 参数指定了关联的列
gene_symbol_entrez <- left_join(data.frame(gene_id = genomewide.genes$gene_id), gene_symbol, by = c("gene_id" = "ENTREZID"))
#将基因 ID 与基因符号进行关联。by 参数指定了关联的列
genomewide.genes$symbol <- gene_symbol_entrez$SYMBOL #将基因符号添加到 genomewide.genes 数据框中
genomewide.genes$description <- gene_symbol_entrez$description #将基因描述添加到 genomewide.genes 数据框中
genomewide.genes$biotype <- gene_symbol_entrez$biotype  #将基因生物类型添加到 genomewide.genes 数据框中

# MSigDb organism
msigdbr_org <- "Homo sapiens" # species


#' # Load A/B data and Replace Centromeric Regions with NaN
# mtx_filtered <- read_tsv(file.path(dir_data, fileNameIn1)) # Filtered
# mtx_filtered <- as.data.frame(mtx_filtered)
mtx_full <- read_tsv(file.path(dir_data, fileNameIn2)) # Full,这里确证了最开始的dir和file2确实是路径一体的
#另外注意一下其他的file文件资料


# Log10-untransform the p-value,注意一下这里的padj在对应的dchic里的解释是怎么样的
#！！！！！！！！！！！！！！！！！！！！暂时就这样指数化
#但是注意这里因为取消了最开始的99999的设置，所以输出是指数化的，但是不确定这样对于后续的结果会不会有影响
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!主要是在数值上，如果有影响的话，就之后再修改回来
mtx_full$padj <- 10^(-mtx_full$padj)
# Chromosome sizes
chrom_sizes <- data.frame(chr = chr, size = seqlengths(getBSgenome(genome = bsgenome, masked = FALSE))[chr])
#注意这里和前面对应的169行是一样的
# Make GRanges
chrom_sizes_gr <- GRanges(seqnames = paste0(chrom_sizes$chr), ranges = IRanges(start = 0, end = chrom_sizes$size))
#创建一个GRanges对象，表示每个染色体的范围

# import centromeric regions from rCGH
hg19_centro<-hg19
# change chromosome indicator to include "chr"
hg19_centro$chrom <- paste0("chr", hg19_centro$chrom)
#但是这里需要修改，将chr23和24修改为x和y
hg19_centro$chrom <- ifelse(hg19_centro$chrom == "chr23", "chrX",
                            ifelse(hg19_centro$chrom == "chr24", "chrY", hg19_centro$chrom))
# make gRanges from centromeric locations 创建一个GRanges对象，表示每个染色体的中心区域的范围
hg19_centro.gr<-makeGRangesFromDataFrame(hg19_centro, seqnames.field = "chrom", start.field = "centromerStart", end.field = "centromerEnd", keep.extra.columns = TRUE)

# make gRanges from the resolution specific coordinates创建一个GRanges对象，表示基因组中的一些特定区域，这里使用的就是dchic结果中区室的特殊区域了
AB.gr<-GRanges(seqnames = mtx_full$chr, IRanges(start = mtx_full$start, end = mtx_full$end))

mcols(AB.gr) <- mtx_full[, c("BT549", "HMEC", "sample_maha", "padj")] #将选定的列设置为GRanges对象的元数据,这些都是原来的mtx中能够查看的
# Remove regions outside of chromosome bounds
olap_chrom <- findOverlaps(AB.gr, chrom_sizes_gr)  #查找两个GRanges对象之间的重叠区域，这里应该找的是我们的dchic数据和前面数据库中指定的染色体范围的数据
AB.gr <- AB.gr[ queryHits(olap_chrom) ] #索引查询查询区域重叠的区域，西安市染色体范围内
# Remove overlaps between the compartments and the centromeric regions
olap_centromere<-findOverlaps(AB.gr, hg19_centro.gr)  #再是中心区域范围内排除，设置为NA
AB.gr$BT549[ queryHits(olap_centromere) ] <- NA
AB.gr$HMEC[ queryHits(olap_centromere) ] <- NA

mtx_full <- as.data.frame(AB.gr)
colnames(mtx_full)[1] <- "chr"
#下面注释的原封不动的!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!值得回过头来观察
# summary(mtx_filtered$padj) # max: 0.01
# summary(mtx_full$padj) # max: 1
# dim(mtx_full[mtx_full$padj <= max(mtx_filtered$padj), ]) # Filtered Full = Filtered
# summary(mtx_filtered$sample_maha)
# summary(mtx_full$sample_maha)
# plot(density(mtx_filtered$glosh))
# plot(density(mtx_full$glosh))
# cor(mtx_filtered$sample_maha, mtx_filtered$UCD52CR - mtx_filtered$UCD52PR)
# View(data.frame(mtx_filtered$sample_maha, mtx_filtered$UCD52CR - mtx_filtered$UCD52PR))


#' # Genome-wide proportion of A/B compartments
# Genome-wide proportion of AB compartments
conditions <- c("HMEC", "BT549") # Conditions

# Matrix to store results创建一个矩阵，用于存储两种条件下的A/B区域的比例,所以是2x2
proportioNumber_AB_genomewide <- matrix(data = 0, ncol = 2, nrow = 2)
colnames(proportioNumber_AB_genomewide) <- conditions
rownames(proportioNumber_AB_genomewide) <- c("Proportion A", "Proportion B")
# For a given condition, calculate the proportion of AB compartments
for (j in conditions) {
  proportioNumber_A <- sum(mtx_full[, j] >  0, na.rm = TRUE) / length(mtx_full[, j][!is.na(mtx_full[, j])] ) #计算了mtx_full数据框中的某一列（由j指定）中大于0的值的数量，并除以该列中非缺失值的总数，以得到A区域的比例
  proportioNumber_B <- sum(mtx_full[, j] <= 0, na.rm = TRUE) / length(mtx_full[, j][!is.na(mtx_full[, j])] ) #计算了mtx_full数据框中的某一列（由j指定）中小于等于0的值的数量，并除以该列中非缺失值的总数，以得到B区域的比例
  proportioNumber_AB_genomewide["Proportion A", j] <- proportioNumber_A
  proportioNumber_AB_genomewide["Proportion B", j] <- proportioNumber_B
}
kable(round(proportioNumber_AB_genomewide, digits = 3))


#' # Chromosome-specific proportion of A/B compartments 对每一条染色体而言进行分析
#' - "UCD52PRA", "UCD52PRB", "UCD52CRA", "UCD52CRB" - the proportion of A and B compartments in each condition.
#' 对应就是"HMEC-A", "HMEC-B", "BT549-A", "BT549-B" - the proportion of A and B compartments in each condition.
#' - "A_log2FC", "B_log2FC" - log2 fold change in the A/B compartment proportions between CR and PR conditions.
#' PR对应HEMC，CR对应BT549
#' - "AB_log2FC" - difference between the A/B ratios in each condition.
#' 
## -------------------------------------------------------------------------------------
# Matrix to store results，存储每个染色体在两个条件下的A/B区域的比例
proportioNumber_AB_chromosome <- matrix(data = 0, ncol = 4, nrow = length(chr))
colnames(proportioNumber_AB_chromosome) <- c("HMECA", "HMECB", "BT549A", "BT549B")
rownames(proportioNumber_AB_chromosome) <- chr
#注意这里的chr是指的前面的1-22，所以没有使用x和y
# For a given condition, calculate the proportion of AB compartments
for (i in chr) {
  for (j in conditions) {
    mtx_full_subset <- mtx_full[mtx_full$chr == i, j]
    proportioNumber_AB_chromosome[i, paste0(j, "A")] <- sum(mtx_full_subset >  0, na.rm = TRUE) / length(mtx_full_subset[!is.na(mtx_full_subset)])
    proportioNumber_AB_chromosome[i, paste0(j, "B")] <- sum(mtx_full_subset <= 0, na.rm = TRUE) / length(mtx_full_subset[!is.na(mtx_full_subset)])
  }
}
#和前面全染色体的没有太大区别
# Convert to data frame
proportioNumber_AB_chromosome <- data.frame(Chromosome = rownames(proportioNumber_AB_chromosome), proportioNumber_AB_chromosome) #加了个染色体行名
# Append ratio of A compartment change添加了一个新的列"A_log2FC"，表示A区域的比例变化,在两个情况中
proportioNumber_AB_chromosome <- proportioNumber_AB_chromosome %>% mutate(A_log2FC = log2(BT549A) - log2(HMECA))
# proportioNumber_AB_chromosome$A_log2FC <- log2(proportioNumber_AB_chromosome$UCD52CRA / proportioNumber_AB_chromosome$UCD52PRA)
# Append ratio of B compartment change
proportioNumber_AB_chromosome <- proportioNumber_AB_chromosome %>% mutate(B_log2FC = log2(BT549B / HMECB))
# Append ratio of A/B ratio compartment change A/B比率在不同情况下的差异
proportioNumber_AB_chromosome <- proportioNumber_AB_chromosome %>% mutate(AB_log2FC = log2(HMECA / HMECB) - log2(BT549A / BT549B) )
# Display interactively
DT::datatable(round_df(proportioNumber_AB_chromosome, 5), options = list(pageLength = 22))

 
#' # Data preparation, filtered
## -------------------------------------------------------------------------------------
# Create data that resembles objects used in 05_AB_eigenvector，目前没有在脚本俩看到这个文件，
#但是后续下面好几行用的都是这个数据进行的绘图，还是跑一遍
#> summary(mtx_full$padj)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0000  0.1591  1.0000  0.6912  1.0000  1.0000 
#这里就是问题了，不知道如何选择padj，要修改的话就直接这里从第79行修改 
#！！！！！！！！！！！！！！！！！！！！！！！！到时候再回过头来修改一下对应的padj数据即可

AB_HMEC <- data.frame(chr   = mtx_full[mtx_full$padj <= padj_compartment_cutoff, "chr"] %>% unlist,
                    start = mtx_full[mtx_full$padj <= padj_compartment_cutoff, "start"] %>% unlist,
                    end   = mtx_full[mtx_full$padj <= padj_compartment_cutoff, "end"] %>% unlist,
                    EV    = mtx_full[mtx_full$padj <= padj_compartment_cutoff, "HMEC"] %>% unlist)
#从 mtx_full 数据框中选取 padj 列小于或等于 padj_compartment_cutoff 的行，并且只包含了 HMEC 和 BT549 列的值。%>% unlist 是用于将结果转换为向量
#总之都是提取特定阈值限制下的数据框，然后分别只提取特定的一列，
#然后注意这里的mtx_full用的是最原始的dchic中获得区室矩阵文件
#所以这里的EV就作为黑匣子，反正EV对应原始dchic区室矩阵文件中HMEC&BT549这两列
#！！！！！！！！！！！！！！！！！！！！！！！！！！！！

AB_BT549 <- data.frame(chr   = mtx_full[mtx_full$padj <= padj_compartment_cutoff, "chr"] %>% unlist,
                    start = mtx_full[mtx_full$padj <= padj_compartment_cutoff, "start"] %>% unlist,
                    end   = mtx_full[mtx_full$padj <= padj_compartment_cutoff, "end"] %>% unlist,
                    EV    = mtx_full[mtx_full$padj <= padj_compartment_cutoff, "BT549"] %>% unlist)
AB_metadata <- data.frame(sample_maha = mtx_full[mtx_full$padj <= padj_compartment_cutoff, "sample_maha"] %>% unlist,
                          padj  = mtx_full[mtx_full$padj <= padj_compartment_cutoff, "padj"] %>% unlist)

# run assignment function on the data为两个数据框 AB_PR 和 AB_CR 添加一个名为 compartment 的新列，并根据 EV 列的值将这些行分配到 "A" 或 "B" 区域
AB_HMEC_assigned <- AB_HMEC; AB_HMEC_assigned$compartment <- ifelse(AB_HMEC_assigned$EV >= 0, "A", "B")
AB_BT549_assigned <- AB_BT549; AB_BT549_assigned$compartment <- ifelse(AB_BT549_assigned$EV >= 0, "A", "B")
#数据框添加一个名为 compartment 的新列。这个新列的值根据 EV 列的值进行分配，如果 EV 列的值大于等于0，则 compartment 列的值为 "A"，否则为 "B
#推测这里EV即原矩阵中对应的文件应该是区室score值之类的，然后这里是新加了一列，依据EV值将对应的区室的类型定下来

# Add dataset signifier to each column for distinguishing the two data sets in combined data frame
#列名添加前缀，以便区分两个数据集在合并数据框中的来源
colnames(AB_HMEC_assigned) <- paste0("HMEC", colnames(AB_HMEC_assigned))
colnames(AB_BT549_assigned) <- paste0("BT549", colnames(AB_BT549_assigned))

# combine the two datasets
AB_conditions_assigned <- cbind(AB_HMEC_assigned, AB_BT549_assigned)

# create a column for the difference in EV for each region between the two datasets
#数据框创建了一个名为 D.EV 的新列，用于存储两个数据集之间 EV 列的差异
AB_conditions_assigned$D.EV <- AB_conditions_assigned$BT549EV - AB_conditions_assigned$HMECEV  #EV是score相减，变化
AB_conditions_assigned$compartment <- paste0(AB_conditions_assigned$HMECcompartment, AB_conditions_assigned$BT549compartment)  #这里是AB.BA等变化
AB_conditions_assigned <- cbind(AB_conditions_assigned, AB_metadata)
#前面从原始mtx数据中提取出来两个子数据，一个就是上面的ABcondition，要进行操作，另外一个是metadata不进行操作



#' # Data preparation, unfiltered，前面所谓filtered，实际上下面这个操作和上面是一模一样的
##上面是依据了阈值，这里没有用到阈值padcutoff，但是实际执行的操作是一模一样的
# Create data that resembles objects used in 05_AB_eigenvector，还是同理，依据上面的数据但是没有变化
AB_HMEC_all <- data.frame(chr   = mtx_full[, "chr"  ] %>% unlist,
                        start = mtx_full[, "start"] %>% unlist,
                        end   = mtx_full[, "end"  ] %>% unlist,
                        EV    = mtx_full[, "HMEC"   ] %>% unlist)
AB_BT549_all <- data.frame(chr   = mtx_full[, "chr"  ] %>% unlist,
                        start = mtx_full[, "start"] %>% unlist,
                        end   = mtx_full[, "end"  ] %>% unlist,
                        EV    = mtx_full[, "BT549"   ] %>% unlist)
AB_metadata_all <- data.frame(sample_maha = mtx_full[, "sample_maha"] %>% unlist,
                          padj  = mtx_full[, "padj"] %>% unlist)

# run assignment function on the data
AB_HMEC_all_assigned <- AB_HMEC_all; AB_HMEC_all_assigned$compartment <- ifelse(AB_HMEC_all_assigned$EV >= 0, "A", "B")
AB_BT549_all_assigned <- AB_BT549_all; AB_BT549_all_assigned$compartment <- ifelse(AB_BT549_all_assigned$EV >= 0, "A", "B")

# Add dataset signifier to each column for distinguishing the two data sets in combined data frame
colnames(AB_HMEC_all_assigned) <- paste0("HMEC", colnames(AB_HMEC_all_assigned))
colnames(AB_BT549_all_assigned) <- paste0("BT549", colnames(AB_BT549_all_assigned))

# combine the two datasets
AB_conditions_assigned_all <- cbind(AB_HMEC_all_assigned, AB_BT549_all_assigned)

# create a column for the difference in EV for each region between the two datasets
AB_conditions_assigned_all$D.EV <- AB_conditions_assigned_all$BT549EV - AB_conditions_assigned_all$HMECEV
AB_conditions_assigned_all$compartment <- paste0(AB_conditions_assigned_all$HMECcompartment, AB_conditions_assigned_all$BT549compartment)
AB_conditions_assigned_all <- cbind(AB_conditions_assigned_all, AB_metadata_all)

#' 
#' 
#' ## Exploratory analysis of eigenvector differences and other measures
#' 
## -------------------------------------------------------------------------------------
DT::datatable(AB_conditions_assigned)
#DT::datatable(AB_conditions_assigned_all)
#就是表格，不知道暂时如何探索


#' # RLE, Run Length Encoding确定连续区域中具有相同标识的区域（例如，从 A 到 A，从 A 到 B 等）
#' Rle for determining consecutive regions of compartments with same identity (e.g. A to A, A to B, etc )
# empty dataframe for the RLE information to be stored
ABrun_all <- data.frame()
# for each chromosome
for (i in chr) {
  # first subset by chromosome
  AB_conditions_assigned_chr <- subset(AB_conditions_assigned, HMECchr == i)
  # run RLE and extract the results into a dataframe. Contains the length of each continuous segment and the identity of the run (AB, BA, AA, BB). Length is in number of segments of each segment type, so multiply by resolution.
  #将运行长度编码的结果转换为数据框，其中包含了连续相同值的长度和这些值的标识
  ABrun <- data.frame(lengths = unclass(rle(AB_conditions_assigned_chr$compartment)$lengths * res_number), compartment = unclass(rle(AB_conditions_assigned_chr$compartment)$values))
  # Fix for the end coordinate of the chromosome not spanning full resolution bin size
  #将 ABrun 数据框中最后一行的 lengths 列的值设置为 ABrun$lengths 列的最后一个值减去分辨率的大小再加上染色体的长度
  ABrun[nrow(ABrun), "lengths"] <- last(ABrun$lengths) - res_number + (max(AB_conditions_assigned_chr$HMECend) - max(AB_conditions_assigned_chr$HMECstart) + 1)
  # rearrange the results of rle to get the coordinates of the AB compartments to correspond to resolution size
  # make the start coordinate zero
  ABrun$start <- 0
  # use the cumulative sum of the lengths at each row to determine the end coordinates of the region range
  ABrun$end <- cumsum(ABrun$lengths)
  # use a lag of one to offset the end coordinates
  lag <- lag(ABrun, 1L)
  # use the offset coordinates in the lag object to set the start coordinates
  ABrun$start <- lag$end
  # fix start coordinates of each range so that they end in 1; otherwise each consecutive start coordinate is the same as the previous region's end coordinate
  ABrun$start <- ABrun$start + 1
  # fix first region's start coordinate, which is an NA, by changing it simply to 1
  ABrun[1, "start"] <- 1
  # add a column with the chr
  ABrun$chr <- i
  # combine the data for each chromosome into a whole genome object
  ABrun_all <- rbind(ABrun_all, ABrun)
}
#这段代码的目的是对 AB_conditions_assigned 数据框中的连续区域进行运行长度编码（RLE），以确定连续区域中具有相同标识的区域（例如，从 A 到 A，从 A 到 B 等）。然后，将结果存储到 ABrun_all 数据框中，并且从 ABrun_all 数据框中提取出 A 区域和 B 区域的信息
# subset the master dataframe of all run lengths to just those for A to A and B to B for downstream gene analysis.
ABrun_AB <- subset(ABrun_all, compartment == "AB")
ABrun_BA <- subset(ABrun_all, compartment == "BA")


#' # Summary of AB compartments
#' - Only significant changes are considered.
#' - Changes within the same compartment type are possible.
#' - "Percent_XX" - percent of AB compartment changes with respect to chromosome length.
#' - "Number_XX" - number of resolution bins of AB compartment changes.
#' - "AB/BA_min/max/med" - min mean and max run lengths for AB and BA compartment changes for each chromosome, calculations are divided by 1000 to convert bases to kb
# empty data frame to contain the summary data
AB_summary <- data.frame(chr = chr)

# modify how the chromosomes in the chom.sizes object are formatted. Adds "chr" to make it compatible with the RLE and combined EV datasets
# chrom.sizes$chr <- paste0("chr", chrom.sizes$chr)

for (i in chr) {
  #print(chr)
  # subset the combined data by chromosome
  AB_conditions_assigned_chr <- subset(AB_conditions_assigned, HMECchr == i)
  # subset the computed RLE data by chromosome
  ABrun <- subset(ABrun_all, chr == i)
  # subset the per chromosome data into AB and BA compartments
  AB_chr <- subset(ABrun, compartment == "AB")
  BA_chr <- subset(ABrun, compartment == "BA")
  # percent of PR data that is in an A compartment  ( # of regions assigned as A times resolution divided by total chromosome size )
  AB_summary$HMEC_A[which(AB_summary$chr == i)] <- round(length(AB_conditions_assigned_chr$HMECcompartment[which(AB_conditions_assigned_chr$HMECcompartment == "A")]) * res_number / chrom.sizes$size[which(chrom.sizes$chr == i)] * 100, 2)
  # percent of PR data that is in an B compartment  ( # of regions assigned as B times resolution divided by total chromosome size )
  AB_summary$HMEC_B[which(AB_summary$chr == i)] <- round(length(AB_conditions_assigned_chr$HMECcompartment[which(AB_conditions_assigned_chr$HMECcompartment == "B")]) * res_number / chrom.sizes$size[which(chrom.sizes$chr == i)] * 100, 2)
  # percent of CR data that is in an A compartment  ( # of regions assigned as A times resolution divided by total chromosome size )
  AB_summary$BT549_A[which(AB_summary$chr == i)] <- round(length(AB_conditions_assigned_chr$BT549compartment[which(AB_conditions_assigned_chr$BT549compartment == "A")]) * res_number / chrom.sizes$size[which(chrom.sizes$chr == i)] * 100, 2)
  # percent of CR data that is in an B compartment  ( # of regions assigned as B times resolution divided by total chromosome size )
  AB_summary$BT549_B[which(AB_summary$chr == i)] <- round(length(AB_conditions_assigned_chr$BT549compartment[which(AB_conditions_assigned_chr$BT549compartment == "B")]) * res_number / chrom.sizes$size[which(chrom.sizes$chr == i)] * 100, 2)
  # Percent of chromosome length as one of four switches: AA, BB, AB, BA.  These are summed together in a fifth category of switch type
  # Percent of chromosome that was A to A compartment compartment
  #将 AB_summary 数据框中 chr 列等于 i 的行的 Percent_AB 和 Percent_BA 列的值设置为染色体为 i 的 PR 数据中 A-B 和 B-A 区域的百分比
  AB_summary$Percent_AA[which(AB_summary$chr == i)] <- round(sum(AB_conditions_assigned_chr$HMECend[AB_conditions_assigned_chr$compartment == "AA"] - AB_conditions_assigned_chr$HMECstart[AB_conditions_assigned_chr$compartment == "AA"]) / chrom.sizes$size[which(chrom.sizes$chr == i)] * 100, 2)
  # Percent of chromosome that was B to B compartment compartment
  AB_summary$Percent_BB[which(AB_summary$chr == i)] <- round(sum(AB_conditions_assigned_chr$HMECend[AB_conditions_assigned_chr$compartment == "BB"] - AB_conditions_assigned_chr$HMECstart[AB_conditions_assigned_chr$compartment == "BB"]) / chrom.sizes$size[which(chrom.sizes$chr == i)] * 100, 2)
  # Percent of chromosome that was A to B compartment compartment
  AB_summary$Percent_AB[which(AB_summary$chr == i)] <- round(sum(AB_conditions_assigned_chr$HMECend[AB_conditions_assigned_chr$compartment == "AB"] - AB_conditions_assigned_chr$HMECstart[AB_conditions_assigned_chr$compartment == "AB"]) / chrom.sizes$size[which(chrom.sizes$chr == i)] * 100, 2)
  # Percent of chromosome that was B to A compartment compartment
  AB_summary$Percent_BA[which(AB_summary$chr == i)] <- round(sum(AB_conditions_assigned_chr$HMECend[AB_conditions_assigned_chr$compartment == "BA"] - AB_conditions_assigned_chr$HMECstart[AB_conditions_assigned_chr$compartment == "BA"]) / chrom.sizes$size[which(chrom.sizes$chr == i)] * 100, 2)
    # Total number of bins with AB compartment changes
  # Total number of bins staying AA 
  AB_summary$Number_AA[which(AB_summary$chr == i)] <- sum(AB_conditions_assigned_chr$compartment == "AA")
  # Total number of bins staying BB 
  AB_summary$Number_BB[which(AB_summary$chr == i)] <- sum(AB_conditions_assigned_chr$compartment == "BB")
  # Total number of bins switching AB
  AB_summary$Number_AB[which(AB_summary$chr == i)] <- sum(AB_conditions_assigned_chr$compartment == "AB")
  # Total number of bins switching BA
  AB_summary$Number_BA[which(AB_summary$chr == i)] <- sum(AB_conditions_assigned_chr$compartment == "BA")
  
  # add min mean and max run lengths for AB and BA compartment changes for each chromosome
  # AB min, mean, and max of run lengths. Calculations are divided by 1000 to convert bases to kb.
  if (nrow(AB_chr) > 0) {
    AB_summary$ABmin[which(AB_summary$chr == i)] <- round(min(AB_chr$lengths) / res_number, 0)
    AB_summary$ABmean[which(AB_summary$chr == i)] <- round(mean(AB_chr$lengths) / res_number, 2)
    AB_summary$ABmax[which(AB_summary$chr == i)] <- round(max(AB_chr$lengths) / res_number, 0)
  } else {
    AB_summary$ABmin[which(AB_summary$chr == i)] <- 0
    AB_summary$ABmean[which(AB_summary$chr == i)] <- 0
    AB_summary$ABmax[which(AB_summary$chr == i)] <- 0
  }
  # BA min, mean, and max of run lengths
  if (nrow(BA_chr) > 0) {
    AB_summary$BAmin[which(AB_summary$chr == i)] <- round(min(BA_chr$lengths) / res_number, 0)
    AB_summary$BAmean[which(AB_summary$chr == i)] <- round(mean(BA_chr$lengths) / res_number, 2)
    AB_summary$BAmax[which(AB_summary$chr == i)] <- round(max(BA_chr$lengths) / res_number, 0)
  } else {
    AB_summary$BAmin[which(AB_summary$chr == i)] <- 0
    AB_summary$BAmean[which(AB_summary$chr == i)] <- 0
    AB_summary$BAmax[which(AB_summary$chr == i)] <- 0
  }
}
# Subset and rename columns
#AB_summary <- AB_summary[, c("chr", "Percent_AA", "Percent_BB", "Percent_AB", "Percent_BA", "Number_AA", "Number_BB", "Number_AB", "Number_BA")]
# Visualize
DT::datatable(AB_summary, options = list(pageLength = 22))
#这里对应的是分析出区室转换的图谱



#' # KaryoploteR Plots and Per Chromosome AB Compartments and EVs
#简单来说就是将每一条chr上的AB区室以及score值映射到基因组chr坐标上的图
# diagnostic plots of the AB compartments and changes. Each data set EV is plotted and the change in EV and color coded to identify switches
for (i in chr) {
  # subset by chromosome
  AB_conditions_assigned_chr <- subset(AB_conditions_assigned, HMECchr == i)
  # subset the per chromosome switch types
  # AB switches
  chr_ab <- subset(AB_conditions_assigned_chr, compartment == "AB")
  # BA switches
  chr_ba <- subset(AB_conditions_assigned_chr, compartment == "BA")
  # AA switches
  chr_aa <- subset(AB_conditions_assigned_chr, compartment == "AA")
  # BB switches
  chr_bb <- subset(AB_conditions_assigned_chr, compartment == "BB")

  # plot the chromosome
  kp <- plotKaryotype(chromosomes = i, genome = "hg19", plot.type = 1)
  # add base numbers to the plot
  kpAddBaseNumbers(kp)

  # PR EV plot，对应HMEC EV plot
  # define r0 and r1, the min and max of the plot area. Allows for easy scaling of the plot area
  r0 <- 0.0
  r1 <- 0.35
  # Add the PR EVs segments and scale to the min and max of the EV values
  kpSegments(kp, chr = i, x0 = AB_conditions_assigned_chr$HMECstart, x1 = AB_conditions_assigned_chr$HMECstart, y0 = 0, y1 = AB_conditions_assigned_chr$HMECEV, r0 = r0, r1 = r1, ymax = max(AB_conditions_assigned_chr$HMECEV, na.rm = T), ymin = min(AB_conditions_assigned_chr$HMECEV, na.rm = T))
  # Add the PR EV plot label,对应HMEC
  kpAddLabels(kp, labels = "HMEC", r0 = r0, r1 = r1, data.panel = 1, label.margin = 0.04, cex = 1.4)
  # Add the axis, scaled for min and max of the PR EVs
  kpAxis(kp, r0 = r0, r1 = r1, cex = 0.65, ymax = max(AB_conditions_assigned_chr$HMECEV, na.rm = T), ymin = min(AB_conditions_assigned_chr$HMECEV, na.rm = T), tick.pos = c(max(AB_conditions_assigned_chr$HMECEV, na.rm = T), 0, min(AB_conditions_assigned_chr$HMECEV, na.rm = T)))
  # Add A and B labels to the axis
  kpAddLabels(kp, labels = "B", r0 = r0, r1 = r0 + ((r1 - r0) / 2), cex = 1.0)
  kpAddLabels(kp, labels = "A", r0 = r0 + ((r1 - r0) / 2), r1 = r1, cex = 1.0)

  # CR EV plot，对应BT549 EV plot
  # define r0 and r1, the min and max of the plot area. Allows for easy scaling of the plot area
  r0 <- 0.4
  r1 <- 0.75
  # Add the CR EVs segments and scale to the min and max of the EV values
  kpSegments(kp, chr = i, x0 = AB_conditions_assigned_chr$BT549start, x1 = AB_conditions_assigned_chr$BT549start, y0 = 0, y1 = AB_conditions_assigned_chr$BT549EV, r0 = r0, r1 = r1, ymax = max(AB_conditions_assigned_chr$BT549EV, na.rm = T), ymin = min(AB_conditions_assigned_chr$BT549EV, na.rm = T), col = "#666666")
  # Add the CR EV plot label
  kpAddLabels(kp, labels = "BT549", r0 = r0, r1 = r1, data.panel = 1, label.margin = 0.04, cex = 1.4)
  # Add the axis, scaled for min and max of the CR EVs
  kpAxis(kp, r0 = r0, r1 = r1, cex = 0.65, ymax = max(AB_conditions_assigned_chr$BT549EV, na.rm = T), ymin = min(AB_conditions_assigned_chr$BT549EV, na.rm = T), tick.pos = c(max(AB_conditions_assigned_chr$BT549EV, na.rm = T), 0, min(AB_conditions_assigned_chr$BT549EV, na.rm = T)))
  # Add A and B labels to the axis
  kpAddLabels(kp, labels = "B", r0 = r0, r1 = r0 + ((r1 - r0) / 2), cex = 1.0)
  kpAddLabels(kp, labels = "A", r0 = r0 + ((r1 - r0) / 2), r1 = r1, cex = 1.0)

  # Delta EV plot (CR EV minus PR EV),对应BT549 EV VS HMEC EV
  # define r0 and r1, the min and max of the plot area. Allows for easy scaling of the plot area
  r0 <- 0.8
  r1 <- 1.1
  # Plot the AA delta EV segments
  kpSegments(kp, chr = i, x0 = chr_aa$HMECstart, x1 = chr_aa$HMECstart, y0 = 0, y1 = chr_aa$D.EV, r0 = r0, r1 = r1, ymax = max(AB_conditions_assigned_chr$D.EV, na.rm = T), ymin = min(AB_conditions_assigned_chr$D.EV, na.rm = T), col = "gray")
  # Plot the BB delta EV segments
  kpSegments(kp, chr = i, x0 = chr_bb$HMECstart, x1 = chr_bb$HMECstart, y0 = 0, y1 = chr_bb$D.EV, r0 = r0, r1 = r1, ymax = max(AB_conditions_assigned_chr$D.EV, na.rm = T), ymin = min(AB_conditions_assigned_chr$D.EV, na.rm = T), col = "gray")
  # Plot the AB delta EV segments
  kpSegments(kp, chr = i, x0 = chr_ab$HMECstart, x1 = chr_ab$HMECstart, y0 = 0, y1 = chr_ab$D.EV, r0 = r0, r1 = r1, ymax = max(AB_conditions_assigned_chr$D.EV, na.rm = T), ymin = min(AB_conditions_assigned_chr$D.EV, na.rm = T), col = "blue")
  # Plot the BA delta EV segments
  kpSegments(kp, chr = i, x0 = chr_ba$HMECstart, x1 = chr_ba$HMECstart, y0 = 0, y1 = chr_ba$D.EV, r0 = r0, r1 = r1, ymax = max(AB_conditions_assigned_chr$D.EV, na.rm = T), ymin = min(AB_conditions_assigned_chr$D.EV, na.rm = T), col = "red")
  # Add axis labels
  kpAddLabels(kp, labels = expression(paste(Delta, "EV")), r0 = r0, r1 = r1, data.panel = 1, label.margin = 0.035, cex = 1.4)
  # Add the axis scaled to the min and max of the delta EV values
  kpAxis(kp, r0 = r0, r1 = r1, cex = 0.65, ymax = max(AB_conditions_assigned_chr$D.EV, na.rm = T), ymin = min(AB_conditions_assigned_chr$D.EV, na.rm = T), tick.pos = c(round(max(AB_conditions_assigned_chr$D.EV, na.rm = T), 3), 0, round(min(AB_conditions_assigned_chr$D.EV, na.rm = T), 3)))

  # Add the legend indicating the color of the switches in the delta EV plot
  legend(x = "top", fill = c("gray", "gray", "blue", "red"), legend = c("AA", "BB", "AB", "BA"), y.intersp = 1, ncol = 2, cex = 1)
  
}
#另外就是这些区域，比如说AA,BB的区域配色要不要改一下
#这些区域明显能够提取出来，在每一条染色体上，然后将这些区室变换的区域拿过去做分析，比如说看看序列特征
#其实也有问题，就是这些区域到底做什么分析，主要是区室都是Mb级别的，肯定是跨区域的，就不能做常规的那些内含子外显子等序列分析
#另外清楚这个EV图的意义


#' 
## ----eval=FALSE-----------------------------------------------------------------------
## # Similar to plots in the above code chunk, but plots only the per chromosome change in EVs and groups AA and BB changes into one color. Easier to see regions that contain AB or BA switches
## 
## for (i in chr) {
##   # subset by chromosome
##   AB_conditions_assigned_chr <- subset(AB_conditions_assigned, PRchr == i)
##   # subset the per chromosome switch types
##   # AB switches
##   chr_ab <- subset(AB_conditions_assigned_chr, compartment == "AB")
##   # BA switches
##   chr_ba <- subset(AB_conditions_assigned_chr, compartment == "BA")
##   # AA switches
##   chr_aa <- subset(AB_conditions_assigned_chr, compartment == "AA")
##   # BB switches
##   chr_bb <- subset(AB_conditions_assigned_chr, compartment == "BB")
## 
##   # plot each chromosome
##   kp <- plotKaryotype(chromosomes = i, genome = "hg19", plot.type = 1)
##   # add base numbers
##   kpAddBaseNumbers(kp)
##   # Delta EV plot (CR EV minus PR EV)
##   # define r0 and r1, the min and max of the plot area. Allows for easy scaling of the plot area
##   r0 <- 0.4
##   r1 <- 0.75
##   # Plot the AA delta EV segments
##   kpSegments(kp, chr = i, x0 = chr_aa$PRstart, x1 = chr_aa$PRstart, y0 = 0, y1 = chr_aa$D.EV, r0 = r0, r1 = r1, ymax = max(AB_conditions_assigned_chr$D.EV, na.rm = T), ymin = min(AB_conditions_assigned_chr$D.EV, na.rm = T), col = "#666666")
##   # Plot the BB delta EV segments
##   kpSegments(kp, chr = i, x0 = chr_bb$PRstart, x1 = chr_bb$PRstart, y0 = 0, y1 = chr_bb$D.EV, r0 = r0, r1 = r1, ymax = max(AB_conditions_assigned_chr$D.EV, na.rm = T), ymin = min(AB_conditions_assigned_chr$D.EV, na.rm = T), col = "#666666")
##   # Plot the AB delta EV segments
##   kpSegments(kp, chr = i, x0 = chr_ab$PRstart, x1 = chr_ab$PRstart, y0 = 0, y1 = chr_ab$D.EV, r0 = r0, r1 = r1, ymax = max(AB_conditions_assigned_chr$D.EV, na.rm = T), ymin = min(AB_conditions_assigned_chr$D.EV, na.rm = T), col = "blue")
##   # Plot the BA delta EV segments
##   kpSegments(kp, chr = i, x0 = chr_ba$PRstart, x1 = chr_ba$PRstart, y0 = 0, y1 = chr_ba$D.EV, r0 = r0, r1 = r1, ymax = max(AB_conditions_assigned_chr$D.EV, na.rm = T), ymin = min(AB_conditions_assigned_chr$D.EV, na.rm = T), col = "red")
##   # Add axis labels
##   kpAddLabels(kp, labels = expression(paste(Delta, "EV")), r0 = r0, r1 = r1, data.panel = 1, label.margin = 0.035, cex = 1.4)
##   # Add the axis scaled to the min and max of the delta EV values
##   kpAxis(kp, r0 = r0, r1 = r1, cex = 0.65, ymax = max(AB_conditions_assigned_chr$D.EV, na.rm = T), ymin = min(AB_conditions_assigned_chr$D.EV, na.rm = T), tick.pos = c(round(max(AB_conditions_assigned_chr$D.EV, na.rm = T), 3), 0, round(min(AB_conditions_assigned_chr$D.EV, na.rm = T), 3)))
##   # Add the legend indicating the color of the switches in the delta EV plot
##   legend(x = "topright", fill = c("#666666", "blue", "red"), legend = c("AA and BB", "AB", "BA"), y.intersp = 0.7, ncol = 1, cex = 0.95)
## }

#' 
#' # Bar charts
#' Summary of bar charts for each dataset and differences between datasets on a per chromosomal basis. 
# generates summary bar charts for each data set and the differences between datasets on a per chromosomal basis.
# Note: some regions have NA for EV so the total doesn't always sum to 100%
# !Need to summarize across all chromosomes once chr23 EVs are computed.

# PR summary stacked bar plots,对应HMEC
# subset the summary file for just the PR relevant data
#堆叠柱状图，展示每个染色体中 A 区域和 B 区域的百分比
HMEC_summary <- AB_summary[, c("chr", "HMEC_A", "HMEC_B")]
# reformat the data to the preferred long format for the stacked bar plots. Also remove "chr" from chromosome names
HMEC_summary <- HMEC_summary %>%
  mutate(chr = factor(row_number())) %>%
  gather(variable, value, -chr)
# Order chromosomes
chr_summary <- aggregate(HMEC_summary$value, by = list(HMEC_summary$chr), sum)
HMEC_summary$chr <- factor(HMEC_summary$chr, levels = chr_summary$Group.1[order(chr_summary$x, decreasing = TRUE)])
# define the plot parameters
HMECsum.plot <- ggplot(HMEC_summary, aes(
  x = chr, y = value,
  fill = factor(variable, levels = c("HMEC_B", "HMEC_A"))
)) +
  geom_bar(stat = "identity")
# plot the stacked bar plots
HMECsum.plot + scale_fill_brewer(palette = "Spectral") + xlab("Chromosome") + ylab("Percent Compartment ID") + ggtitle("HMEC AB Compartment switches By Chromosome") + guides(fill = guide_legend(title = "ID Type"))


# CR summary stacked bar plots,对应BT549
# subset the summary file for just the CR relevant data
BT549_summary <- AB_summary[, c("chr", "BT549_A", "BT549_B")]
# reformat the data to the preferred long format for the stacked bar plots. Also remove "chr" from chromosome names
BT549_summary <- BT549_summary %>%
  mutate(chr = factor(row_number())) %>%
  gather(variable, value, -chr)
# Order chromosomes
chr_summary <- aggregate(BT549_summary$value, by = list(BT549_summary$chr), sum)
BT549_summary$chr <- factor(BT549_summary$chr, levels = chr_summary$Group.1[order(chr_summary$x, decreasing = TRUE)])
# define the plot parameters
BT549sum.plot <- ggplot(BT549_summary, aes(
  x = chr, y = value,
  fill = factor(variable, levels = c("BT549_B", "BT549_A"))
)) +
  geom_bar(stat = "identity")
# plot the stacked bar plots
BT549sum.plot + scale_fill_brewer(palette = "Spectral") + xlab("Chromosome") + ylab("Percent Compartment ID") + ggtitle("BT549 AB Compartments switches By Chromosome") + guides(fill = guide_legend(title = "ID Type"))


# compartment switches summary stacked bar plots including percent of regions that switched.
conditions_summary <- AB_summary[, c("chr", "Percent_AA", "Percent_BB", "Percent_AB", "Percent_BA")]
# reformat the data to the preferred long format for the stacked bar plots. Also remove "chr" from chromosome names
conditions_summary <- conditions_summary %>%
  mutate(chr = factor(row_number())) %>%
  gather(variable, value, -chr)
# Sort by largest-to-smallest switches
# Proportions for any switch, to order chromosomes
if (!AB_BA) {
  # Use total AB changes for the genome-wide results
  conditions_summary_AB_BA <- conditions_summary
} else {
  # Use only AB and BA changes for significant results
  conditions_summary_AB_BA <- conditions_summary[conditions_summary$variable %in% c("Percent_AB", "Percent_BA"), ]
}
# Sum them up for any type of switch
conditions_summary_AB_BA <- aggregate(conditions_summary_AB_BA$value, by = list(conditions_summary_AB_BA$chr), FUN = sum)
colnames(conditions_summary_AB_BA) <- c("chr", "percent")
# Order chromosomes
conditions_summary$chr <- factor(conditions_summary$chr, levels = conditions_summary_AB_BA$chr[order(conditions_summary_AB_BA$percent, decreasing = TRUE)])
# define the plot parameters
conditions_AB <- ggplot(conditions_summary, aes(
  x = chr, y = value,
  fill = factor(variable, levels = c("Percent_BA", "Percent_AB", "Percent_BB", "Percent_AA"))
)) +
  geom_bar(stat = "identity")
# plot the stacked bar plots
conditions_AB + scale_fill_brewer(palette = "Spectral") + xlab("Chromosome") + ylab("Percent Compartment ID") + ggtitle("AB Compartment switches By Chromosome") + guides(fill = guide_legend(title = "ID Type"))
#ggsave(paste0("HMEC_vs_BT549_dcHiC_AB_switches_per_chr_", padj_compartment_cutoff, ".svg"), width = 7, height = 2.5)

#' 
## ----fig.height=6---------------------------------------------------------------------
#' Plot the proportion of a selected AB compartment change, 
#' chromosomes sorted from larger to smaller fraction
#' @param conditions conditions to visualize. Any combination of "Percent_AA", "Percent_BB", "Percent_AB", "Percent_BA"
#' @param cols colors. Should be the same number as conditions
plot_condition <- function(conditions = c("Percent_AA"), cols = "#99D594") {
  # Select data
  conditions_summary <- AB_summary[, c("chr", conditions)]
  # Reformat into long form
  condition_summary_long <- melt(conditions_summary, id = "chr")
  condition_summary_long$chr <- factor(condition_summary_long$chr, levels = conditions_summary$chr[order(conditions_summary %>% dplyr::select(starts_with("Percent")) %>% apply(., 1, sum), decreasing = FALSE)])
  # display.brewer.pal(7, "Spectral")
  # brewer.pal(7, "Spectral")
  p <- ggplot(condition_summary_long, aes(x = chr, y = value, fill = factor(variable))) +
    geom_bar(stat = "identity") +
    theme_bw() +
    scale_fill_manual(values = cols) + 
    ylab("Percent of the chromosome length") +
    coord_flip()
  p
}
#定义了一个名为 plot_condition 的函数，用于绘制染色体上特定 AB 区域变化的百分比。函数接受两个参数：conditions 和 cols。
#conditions 参数是一个字符向量，用于指定要可视化的条件，可以是 "Percent_AA"、"Percent_BB"、"Percent_AB" 或 "Percent_BA" 中的任意组合。cols 参数是一个字符向量，用于指定绘图所用的颜色。
#函数首先从 AB_summary 数据框中选择染色体和指定的条件列，并将其转换为长格式。然后，它对染色体进行排序，以便按照染色体的大小排列。最后，它使用 ggplot2 包绘制了堆叠柱状图，其中 x 轴表示染色体，y 轴表示百分比，fill 参数表示区域的标识（A 或 B）
#下面就是这个函数具体的参数输入应用
plot_AA <- plot_condition(conditions = c("Percent_AA"), cols = "#99D594")
plot_BB <- plot_condition(conditions = c("Percent_BB"), cols = "#FC8D59")
plot_AB <- plot_condition(conditions = c("Percent_AB"), cols = "#FEE08B")
plot_BA <- plot_condition(conditions = c("Percent_BA"), cols = "#D53E4F")

grid.arrange(plot_AA, plot_BB, plot_AB, plot_BA, ncol = 2)
#上面绘制的就是AA,AB,BA,BB等转换区域分别对应各自chr的比例



## ----fig.height=3.5, fig.width=5------------------------------------------------------
# All conditions
plot_all <- plot_condition(conditions = c("Percent_AA", "Percent_BB", "Percent_AB", "Percent_BA"), cols = c("#99D594", "#FC8D59", "#FEE08B", "#D53E4F"))
plot_all
#注意和前面的占据区室相互区分（这里是在chr长度比例基础上，前面是在区室内部比例基础上）
#这里是合并在一起占据chr长度比例，不是区室内部比例

#' # Piecharts of compartment switches 
## ----fig.height=6---------------------------------------------------------------------
# Get the number of compartments for each switch 
sum_AA <- sum(AB_summary$Number_AA)
sum_AB <- sum(AB_summary$Number_AB)
sum_BB <- sum(AB_summary$Number_BB)
sum_BA <- sum(AB_summary$Number_BA)
slices <- c(sum_AA, sum_BB, sum_AB, sum_BA)

# transform the counts to percentages for the legend 
pct <- round(slices/sum(slices)*100, digits = 2)
lbls <- c("Percent_AA", "Percent_BB", "Percent_AB", "Percent_BA")
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # add % to labels

#svg(paste0("manuscript/figures/dcHiC_pie_", padj_compartment_cutoff, ".svg"), width = 4.5, height = 4.5)
pie(slices,labels = lbls, col=rainbow(length(lbls)), init.angle = 70,
   main="Compartment changes in HMEC & BT549 conditions") 

#Compartment changes in HMEC & BT549 conditions, excluding non-changed regions
dev.off()


#' # Eigenvector distribution per compartment change
#' #这里就是对EV值，也就是前面使用的区室score进行了一个部分的检验
#' #所以EV值是特征向量值，本质上就是对应区室score值
## -------------------------------------------------------------------------------------
aggregated_ev <- data.frame(Condition = c(rep("HMEC", nrow(AB_conditions_assigned)), rep("BT549", nrow(AB_conditions_assigned))),
                            EV = c(AB_conditions_assigned$HMECEV, AB_conditions_assigned$BT549EV),
                            Switch = c(AB_conditions_assigned$compartment, AB_conditions_assigned$compartment))

ggplot(aggregated_ev, aes(x=Switch, y=EV, fill=Condition)) +
  geom_violin(position=position_dodge(1)) +
  geom_boxplot(width=0.1, position=position_dodge(1))



#' # Genes and Gene Enrichments
#' ## Obtain the ranked genes to be used in GSEA
#' Function obtained from `05_AB_eigenvector.Rmd`
## -------------------------------------------------------------------------------------
# function to Find genes in AB or BA regions by detecting overlap between those regions and gene TSS sites in hg19 from above
#区室切换与gene相互联系，实质上就是寻找区室切换的靶基因，寻找的是区室切换区域重叠的转录起始位点TSS
detect_genes_ranked <- function(rle_dataset) {
  # convert the compartment switch ranges to GRanges object
  rle_dataset.gr <- makeGRangesFromDataFrame(rle_dataset, seqnames.field = "chr", start.field = "start", end.field = "end", keep.extra.columns = TRUE)
  # detect overlaps with hg19 reference genes
  rle_dataset_gene_olap <- findOverlaps(rle_dataset.gr, genomewide.genes)
  # extract the overlaps as gene symbols
  rle_dataset_genes <- rle_dataset_gene_olap %>%
    as.data.frame() %>%
    # group_by(queryHits) %>%
    mutate(genes = genomewide.genes$gene_id[subjectHits]) %>%
    dplyr::select(queryHits, genes) %>%
    distinct()
  # temporary data frame containing the region data
  tmpAB <- rle_dataset %>%
    dplyr::select(chr, start, end, D.EV, switch) %>%
    mutate(id = 1:nrow(rle_dataset))
  # combine gene symbols and regions
  rle_dataset_genes <- left_join(rle_dataset_genes, tmpAB, by = c("queryHits" = "id"))
  # Get gene names
  gene_symbols <- bitr(rle_dataset_genes$genes, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = OrgDb)
  rle_dataset_genes <- left_join(rle_dataset_genes, gene_symbols, by = c("genes" = "ENTREZID"))
  # Aggregate by maximum absolute Eigenvector difference
  rle_dataset_genes <- aggregate(D.EV ~ SYMBOL, rle_dataset_genes, max)
  # Append description
  rle_dataset_genes <- left_join(rle_dataset_genes, gene_annotations, by = c("SYMBOL" = "symbol"))
  # rearrange the columns
  rle_dataset_genes <- rle_dataset_genes[, c("SYMBOL", "D.EV", "description", "biotype")]
  colnames(rle_dataset_genes)[1] <- "genes" # Rename the first column
  # Order by largest to smallest difference
  rle_dataset_genes <- rle_dataset_genes[order(rle_dataset_genes$D.EV, decreasing = TRUE), ]
  # Return the dataset
  rle_dataset_genes 
}
#上面其实是定义了一个函数

#detect_genes_ranked 函数旨在通过检测AB或BA区域与hg19参考基因组中基因TSS位点之间的重叠来查找AB或BA区域中的基因。此函数接受一个包含区室切换范围的数据集 (rle_dataset)，并执行以下步骤：
#将区室切换范围转换为 GRanges 对象。
#检测与 hg19 参考基因的重叠。
#提取重叠作为基因符号。
#结合基因符号和区域。
#使用 org.Hs.eg.db 数据库中的 bitr 函数获取基因名称。
#按最大绝对特征向量差异进行聚合。
#添加基因描述和生物类型。
#重新排列列。
#按最大到最小差异对基因进行排序。
#返回数据集。

AB_conditions_assigned_subset <- AB_conditions_assigned[, c("HMECchr", "HMECstart", "BT549end", "D.EV", "compartment")]
colnames(AB_conditions_assigned_subset) <- c("chr", "start", "end", "D.EV", "switch")

genes_AA_ranked <- detect_genes_ranked(AB_conditions_assigned_subset[AB_conditions_assigned_subset$switch == "AA", ])
genes_BB_ranked <- detect_genes_ranked(AB_conditions_assigned_subset[AB_conditions_assigned_subset$switch == "BB", ])
genes_AB_ranked <- detect_genes_ranked(AB_conditions_assigned_subset[AB_conditions_assigned_subset$switch == "AB", ])
genes_BA_ranked <- detect_genes_ranked(AB_conditions_assigned_subset[AB_conditions_assigned_subset$switch == "BA", ])

# All regions
AB_conditions_assigned_all_subset <- AB_conditions_assigned_all[, c("HMECchr", "HMECstart", "BT549end", "D.EV", "compartment")]
colnames(AB_conditions_assigned_all_subset) <- c("chr", "start", "end", "D.EV", "switch")

genes_all_ranked <- detect_genes_ranked(AB_conditions_assigned_all_subset)

#' 
#' ## Oncogene annotations
#' 
## -------------------------------------------------------------------------------------
# Add onco annotations to the gene lists
# Import PCAWG_C, PCAWG_N, COSMIC oncogene annotation lists将肿瘤基因注释添加到基因列表中。它导入了四个不同来源的肿瘤基因注释列表
#https://github.com/mdozmorov/Cancer_notes/blob/master/README.md 这里有对应解释的各种gene数据的来源，总之都是收集来的cancer gene注释
PCAWG_C <- read.table(file = "/mnt/disk4/haitao/bysj_seu/geo_data/hic/script6/PCAWG_2020_PID_C_87_genes.txt", header = FALSE) #来自 PCAWG 项目的肿瘤基因列表（Cohort）
PCAWG_N <- read.table(file = "/mnt/disk4/haitao/bysj_seu/geo_data/hic/script6/PCAWG_2020_PID_N_93_genes.txt", header = FALSE) #来自 PCAWG 项目的肿瘤基因列表（Normal）
#全基因组泛癌症分析(PCAWG)研究是一项国际合作，旨在识别来自国际癌症基因组联盟的2600多个癌症全基因组的共同突变模式
#The Pan-Cancer Analysis of Whole Genomes (PCAWG) study

COSMIC <- read.table(file = "/mnt/disk4/haitao/bysj_seu/geo_data/hic/script6/COSMIC_genes.txt", header = FALSE)
#Cancer Gene Census (CGC), download COSMIC
BushmanLab <- read.table(file = "/mnt/disk4/haitao/bysj_seu/geo_data/hic/script6/allOnco_May2018.tsv")  #来自 Bushman 实验室的肿瘤基因列表
#切记，注意在读入github上的文件的时候，不能直接读入https的文件，要点击在右上角的raw，然后会转到一个raw相关的http的url，然后直接复制粘贴这个地址即可
#但是要注意进入不进入的了，所以还是梯子的问题
#上面的文件都读入不了，所以只能直接下载到本地再读取


# function add the annotations to gene lists
add_onco <- function(dataset) {
  # DEGs detected in COSMIC
  degs_cosmic <- ifelse(dataset$genes %in% intersect(dataset$genes, COSMIC$V1), "Yes", "")
  degs_bushman <- ifelse(dataset$genes %in% intersect(dataset$genes, BushmanLab$V1), "Yes", "")
  # DEGs detected in PID-C set
  degs_PID_C <- ifelse(dataset$genes %in% intersect(dataset$genes, PCAWG_C$V1), "Yes", "")
  # DEGs detected in PID-N set
  degs_PID_N <- ifelse(dataset$genes %in% intersect(dataset$genes, PCAWG_N$V1), "Yes", "")

  # Attach these annotations
  dataset <- data.frame(dataset, COSMIC = degs_cosmic, BushmanLab = degs_bushman, PID_C = degs_PID_C, PID_D = degs_PID_N)
  # dataset <- dplyr::arrange(dataset, desc(COSMIC), desc(BushmanLab), desc(PID_C), desc(PID_D))
}
#标记数据集中的基因是否在这些列表中。然后，它将这些逻辑向量添加到数据集中作为新的列，并返回修改后的数据集
#了解清楚这些gene list很重要


# add the gene lists，所以是区室转换中的gene+来自gene list 4个来源限制的gene
genes_AA_ranked <- add_onco(genes_AA_ranked)
genes_BB_ranked <- add_onco(genes_BB_ranked)
genes_AB_ranked <- add_onco(genes_AB_ranked)
genes_BA_ranked <- add_onco(genes_BA_ranked)

genes_all_ranked <- add_onco(genes_all_ranked)

#' 
#' # Save the data
#' 
#' Save the gene lists with differential eigenvector values sorted, will be used in GSEA analysis 
#' 
## -------------------------------------------------------------------------------------
# Write out AB and BA gene/compartment data to Xcel spreadsheets
if (rerun) {
  x <- c(list(GenesAA = genes_AA_ranked), list(GenesBB = genes_BB_ranked), list(GenesAB = genes_AB_ranked), list(GenesBA = genes_BA_ranked), list(GenesAll = genes_all_ranked), list(Summary = AB_summary), list(EVs_and_Compartments = AB_conditions_assigned), list(RLE = ABrun_all))
  write_xlsx(x, path = fileNameOut1)
}
#这个文件在/mnt/disk4/haitao/bysj_seu/geo_data/hic/script6/HMEC_vs_BT549_100kb_dchic/DifferentialResult/results/中


#' 
#' 
#' # Save BED files
#' 
#' The following format will use colors, https://www.biostars.org/p/269367/
#' 
#' ```
#' track itemRgb=On								
#' chr1	1000	2000	AB	0	.	1000	2000	167,203,104
#' chr1	1500	2500	BB	0	.	1000	2000	242,51,16
#' ```
#' 
#' - "chr", "start", "end" - coordinates
#' - "compartment" - type of AB compartment change
#' - "padj - difference in compartments
#' - "." - strand
#' - "thickStart" and "thickEnd" 
#' - RGB colors
#'     - AA - red, 255,0,0
#'     - AB - orange, 255,165,0
#'     - BA - light purple, 158,121,240
#'     - BB - blue, 0,0,255
#' 
## -------------------------------------------------------------------------------------
x <- data.frame(chr        = AB_conditions_assigned$HMECchr,
                start      = AB_conditions_assigned$HMECstart,
                end        = AB_conditions_assigned$HMECend,
                name       = AB_conditions_assigned$compartment,
                score      = -log10(AB_conditions_assigned$padj),
                strand     = ".",
                thickStart = AB_conditions_assigned$BT549start,
                thickEnd   = AB_conditions_assigned$BT549end,
                color = ifelse(AB_conditions_assigned$compartment == "AA", "255,0,0",
                               ifelse(AB_conditions_assigned$compartment == "AB", "255,165,0",
                                      ifelse(AB_conditions_assigned$compartment == "BA", "158,121,240",
                                             ifelse(AB_conditions_assigned$compartment == "BB", "0,0,255", ""))))
)
# Keep full numbers
x$start <- format(x$start, scientific = FALSE, trim = TRUE, justify = "none")
x$end   <- format(x$end, scientific = FALSE, trim = TRUE, justify = "none")
x$thickStart <- format(x$thickStart, scientific = FALSE, trim = TRUE, justify = "none")
x$thickEnd   <- format(x$thickEnd, scientific = FALSE, trim = TRUE, justify = "none")
# Append header
x <- rbind(c("track itemRgb=On", rep("", ncol(x) - 1)), x)
# Save the BED file
if (rerun) {
  fwrite(round_df(x), file = fileNameOut2, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}
#保存在/mnt/disk4/haitao/bysj_seu/geo_data/hic/script6/HMEC_vs_BT549_100kb_dchic/DifferentialResult/results/

#' 
#' ## Summary and EV file
#' 
#' **Summary:** Summary of per chromosome statistics. "chr"- chromosome number. "PT_A" - Percent of chromosome in A compartment in primary tumor dataset. "PT_B" - Percent of chromosome in A compartment in primary tumor dataset. "CR_A" - Percent of chromosome in A compartment in CR dataset. "CR_B" - Percent of chromosome in A compartment in CR dataset. "Percent_AA", "Percent_BB", "Percent_AB", and "Percent_BA" - percent of each chromosome that is in one of the four compartment switches between the two datasets. "AB_min", "AB_mean", "AB_max" - Min mean and max length of AB compartment per chromosome. "BA_min", "BA_mean", "BA_max" - Min mean and max length of BA compartment per chromosome. 
#' 
#' **EVs_and_Compartments** A summary of regions on each chromosome and their associated EV and compartment. For both PR and CR datasets. "PRchr" - chromosome in PR data. "PRstart" - start coordinates of the region in PR dataset. "PRend" end coordinate of the region in the PR data. "PREV" EV for the region in PR data. "compartment" Either A or B, the assigned compartment for the region. "CRchr" - chromosome in CR data. "CRstart" - start coordinates of the region in CR dataset. "CRend" end coordinate of the region in the CR data. "CREV" EV for the region in CR data. "compartment" Either A or B, the assigned compartment for the region. "D.EV" the difference in the EV between datasets. Determined by CREV - PREV. "switch" whether the region switched from A to A (AA), A to B (AB), A to A (AA), or B to B (BB).
#' 
#' **RLE** - Summarizes the run length of regions and their associated switches. "lengths" - length of the run. "comparment" what compartment the region belongs to. "start" start coordinate of the region. "end" coordinate of the region. "chr" chromosome that the region belongs to.  
#' 
#' # EDA of eigenvectors

#摘要（Summary）：每个染色体的统计摘要。"chr" - 染色体编号。"PT_A" - 原发肿瘤数据集中 A 区室的百分比。"PT_B" - 原发肿瘤数据集中 B 区室的百分比。"CR_A" - CR 数据集中 A 区室的百分比。"CR_B" - CR 数据集中 B 区室的百分比。"Percent_AA"、"Percent_BB"、"Percent_AB" 和 "Percent_BA" - 每个染色体中四种区室切换的百分比。"AB_min"、"AB_mean" 和 "AB_max" - 每个染色体中 AB 区室的最小、平均和最大长度。"BA_min"、"BA_mean" 和 "BA_max" - 每个染色体中 BA 区室的最小、平均和最大长度。
#EVs_and_Compartments：每个染色体上区域的摘要及其相关的 EV 和区室。对于 PR 和 CR 数据集。"PRchr" - PR 数据中的染色体。"PRstart" - PR 数据集中区域的起始坐标。"PRend" - PR 数据集中区域的结束坐标。"PREV" - PR 数据集中区域的 EV。"compartment" - 区域的分配区室（A 或 B）。"CRchr" - CR 数据中的染色体。"CRstart" - CR 数据集中区域的起始坐标。"CRend" - CR 数据集中区域的结束坐标。"CREV" - CR 数据集中区域的 EV。"compartment" - 区域的分配区室（A 或 B）。"D.EV" - 数据集之间 EV 的差异。由 CREV - PREV 确定。"switch" - 区域是否从 A 到 A（AA）、A 到 B（AB）、A 到 A（AA）或 B 到 B（BB）。
#RLE：区域的运行长度及其相关的切换摘要。"lengths" - 区域的长度。"compartment" - 区域所属的区室。"start" - 区域的起始坐标。"end" - 区域的结束坐标。"chr" - 区域所属的染色体。





#算了，下面的步骤就不做了，原始dchic执行命令不清楚，不清楚那个tsv文件如何对应上
#！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
#' Independent code, for exploratory analysis of eigenvector concordance.
#' 
## ----fig.height=2, fig.width=2--------------------------------------------------------
library(data.table)
library(pheatmap)
# Full dcHiC results, root folder
dir_data <- "/mnt/disk4/haitao/bysj_seu/geo_data/hic/script6/HMEC_vs_BT549_100kb_dchic"
#/Users/mdozmorov/Documents/Data/GoogleDrive/HiC_files/results/AB_compartments/dcHiC_2021-09-03/CR_PR_Dataset_250Kb类比
#
# bedGraph folders
dir_bedGraph_HMEC <- file.path(dir_data, "HMEC_rep2_100000_pca/intra_pca/HMEC_rep2_100000_mat")
#PR_250Kb_pca/intra_pca/PR_250Kb_mat类比
#实际上就是"/mnt/disk4/haitao/bysj_seu/geo_data/hic/script6/HMEC_vs_BT549_100kb_dchic/HMEC_rep2_100000_pca/intra_pca/HMEC_rep2_100000_mat/"

#注意我这里实际上只有2个rep的数据，所以只能还是取rep2的数据
dir_bedGraph_BT549 <- file.path(dir_data, "BT549_rep2_100000_pca/intra_pca/BT549_rep2_100000_mat")
#CR_250Kb_pca/intra_pca/CR_250Kb_mat类比
#实际上就是"/mnt/disk4/haitao/bysj_seu/geo_data/hic/script6/HMEC_vs_BT549_100kb_dchic/BT549_rep2_100000_pca/intra_pca/BT549_rep2_100000_mat/"


# Selected PCs
pc_selected <- fread(file.path(dir_data, "chr_pc_selected.tsv"))
#主要会是这一步读取的文件不知道在哪里

pc_selected <- as.data.frame(pc_selected)
# Discordant selection
print("Discordant PC selection")
kable(pc_selected[pc_selected$HMEC != pc_selected$BT549, ])

# bedGraph correlation for each chromosome
chromosomes <- paste0("chr", c(1:22, "X"))
pcs <- c("PC1", "PC2", "PC3")
for (chr in chromosomes) {
  print(paste("Selected:", pc_selected[pc_selected$chr == chr, ]))
  mtx <- matrix(0, nrow = length(pcs), ncol = length(pcs))
  for (i in 1:length(pcs)) {
    for (j in 1:length(pcs))
      if(i >= j) {
        pc_PR <- fread(file.path(dir_bedGraph_HMEC, paste0(chr, ".", pcs[i], ".bedGraph")))
        pc_CR <- fread(file.path(dir_bedGraph_BT549, paste0(chr, ".", pcs[j], ".bedGraph")))
        if (nrow(pc_HMEC) != nrow(pc_BT549)) {
          print(paste("Length mismatch:", chr, "HMEC", pcs[i], "BT549", pcs[j]))
        } else {
          mtx[i, j] <- mtx[j, i] <- cor(pc_HMEC$V4, pc_BT549$V4)
        }
      }
  }
  colnames(mtx) <- rownames(mtx) <- paste(chr, pcs, sep = ".")
  pheatmap(mtx, display_numbers = TRUE, cluster_cols = FALSE, cluster_rows = FALSE, legend = FALSE, treeheight_col = 0, treeheight_row = 0)
}

#' 
