setwd("/mnt/disk4/haitao/bysj_seu/geo_data/hic/script5/cool_scc/")  
library(tidyverse)
library(plyr)
library(stringi)

mydir = "/mnt/disk4/haitao/bysj_seu/geo_data/hic/script5/cool_scc/20kb/" #这里设置含有hicrep的txt文件的文件夹
myfiles = list.files(path=mydir, pattern = "*.txt", full.names=TRUE)
dat_csv = lapply(myfiles, read_csv) #read all txt files



df <- data.frame(dat_csv) #dataframe
df <- df[-1,] #remove first row: #
rownames(df) <- c(1:24) #rename rownames,刚好就是chr1-chr24（x和y也算进去了）
#dim(df)查看维度
#summary(df)

cn <- c() #contains all the pairwise comparison name，包含所有配对的比对文件夹
n_col <- colnames(df) #names of column，下面就是从列名中提取子字符串，按照自己实际来处理
#注意下列代码中的rep名称的修改，需要进入到hicrep的脚本输出中手动修改第一行的名字，直到修改为rep.txt格式，以及在rep之间使用双横线
for (i in 1:length(n_col)){ 
  name <- n_col[i]
  cn <- append(cn, sub("\\.txt","",str_sub(name, 4)))  #.*匹配任意字符，\\.需要转义
  #这里是提取从倒数第15个字符到倒数第5个字符的子字符串（包括这两个位置的字符），主要是最长的是MB231和HCC70,得提取到16，两道工序
}
colnames(df) <- cn


df_chr1 <- df[1,] #contains the first row of df (chromosome 1)，即染色体1

col_1 <- c() #contains the first sample (pair-wise)，比对一边的样本
col_2 <- c() #contains the 2nd sample (pair-wise)，比对另外一边的样本
scc <- c() #contains scc value for the pair-wise comparison，该比对的scc值
for (i in 1:length(df_chr1)){
  str <- colnames(df_chr1)[i]
  two_col <- str_split_fixed(str, "__", 2)
  col_1 <- append(col_1, two_col[1])
  col_2 <- append(col_2, two_col[2])
  scc <- append(scc, df_chr1[i])
}
#本意是1,2,3比最后是1vs2,2vs3,3vs1，这样每一边都取一个值得到123不重复，
#如果hicrep中是j>i的话，那比对的每一边都有可能不完整且重复，所以需要下面df翻转之后再合并，但是自身比对的数据就没法获取，得自己添加
#如果hicrep中是j>=i的话，每一边都完整但会有重复，如果这样的话其实就没必要df翻转再合并了
#如何让每一边都只留下一个rep名字

col_1 <- data.frame(col_1)
col_2 <- data.frame(col_2)
scc <- data.frame(scc)  #这里实际上scc就是提取了df的第一行而已
scc <- data.frame(t(scc))  #变成了一列
colnames(scc) <- c("scc") #上面转置过后列名实际上变成了t.scc.,这里只是把列名改回去scc


df_chr1_hm <- data.frame()
df_chr1_hm <- cbind(col_1, col_2, scc) #dataframe with all the values，横向合并

#无论是41还是40都建议执行下列操作
df_chr1_rev <- data.frame()
df_chr1_rev <- cbind(col_2, col_1, scc) #dataframe with all the values, in reverse，横向合并，但是与上面相反
colnames(df_chr1_rev) <- c("col_1", "col_2", "scc") #to add columns into correct columns in rbind，这一步故意反着来，修改了列名
df_chr1_hm <- rbind(df_chr1_hm, df_chr1_rev) #all values，纵向合并，就解决了原始不全会重复的问题，因为这样就双方都有

#下面1两行其实就是修改双方矩阵比对中的样本名字，没必要，主要是自己命名得挺好的，之后如果名字不好出现问题就选这一步
#df_chr1_hm$col_1 <- str_replace_all(df_chr1_hm$col_1, c("105259" = "PR1_105259", "102770" = "PR2_102770", "102933" = "PR3_102933", "105242" = "CR1_105242", "105244" = "CR2_105244", "105246" = "CR3_105246", "100887" = "LM1_100887", "100888" =  "LM2_100888", "100889" = "LM3_100889"))
#df_chr1_hm$col_2 <- str_replace_all(df_chr1_hm$col_2, c("105259" = "PR1_105259", "102770" = "PR2_102770", "102933" = "PR3_102933", "105242" = "CR1_105242", "105244" = "CR2_105244", "105246" = "CR3_105246", "100887" = "LM1_100887", "100888" =  "LM2_100888", "100889" = "LM3_100889"))
#df_chr1_hm$col_1 <- rbind(df_chr1_hm$col_1, df_chr1_hm$col_2)

df_chr1_hm$scc <- as.numeric(as.character(df_chr1_hm$scc)) #to make scc continuous (not discrete); helps with making heatmap
#class查看之后scc这列是字符型，如果 scc 列是字符型的，那么它会被视为离散的因子（factor），而不是连续的数值


# make heatmap for chr1，仅仅是chr1
#如果是40，没有自身比对的话需要自己添加对角线，如果是41的话就不用
plt_chr1 <- ggplot(df_chr1_hm, aes(x=col_1, y=col_2, fill=scc)) + geom_tile() + geom_text(aes(label = round(scc, 1))) + coord_fixed() + scale_fill_gradient(low = "blue", high = "red") + theme(axis.text.x = element_text(angle = 90), axis.title.x=element_blank(), axis.title.y=element_blank()) + ggtitle("Chromosome 1")
ggsave("scc_heatmap_chr1.jpg", plot=plt_chr1)


# average,上面是染色体1，下面要全局均值,因为代码的一些执行漏洞，还是建议chr1——》全局都执行下来

df_num <- lapply(df, as.numeric)
df_num <- data.frame(df_num)
df_mean <- data.frame(colMeans(df_num)) #每一列的均值，实际上就是每一对比较的所有chr的均值,此时df_mean是列向量


#下面应该是和最开始提取的样本名字一样的操作，但是上面chr1已经操作过了，就不用再做了
#而且其实是可以像上面那样将所有染色体上面的数据都显示一次
#rn <- c() #contains all the pairwise comparison name，不止局限于chr1
#n_row <- rownames(df_mean) #names of column
#for (i in 1:length(n_row)){ 
#  name <- n_row[i]
#  rn <- append(rn, str_sub(name, 2, -1))
#}
#rownames(df_mean) <- rn


#下面取名字对矩阵的操作完全就是一样，和上面chr1一样操作
col_1 <- c() #contains the first sample (pair-wise),同上chr1时候的第一个样本rep
col_2 <- c() #contains the 2nd sample (pair-wise)
scc <- c() #contains scc value for the pair-wise comparison
for (i in 1:nrow(df_mean)){
  str <- rownames(df_mean)[i]
  two_col <- str_split_fixed(str, "__", 2)
  col_1 <- append(col_1, two_col[1])
  col_2 <- append(col_2, two_col[2])
  scc <- append(scc, df_mean[i,])
}
#此处应该取得df_mean每行的值,此时scc和col一样是列向量,但是非数据框

col_1 <- data.frame(col_1) 
col_2 <- data.frame(col_2)
scc <- data.frame(scc)  #全都正式变成列向量的数据框

df_mean_hm <- data.frame()
df_mean_hm <- cbind(col_1, col_2, scc) #dataframe with all the values，按照列合并

#下面同理翻转再合并
df_mean_rev <- data.frame()
df_mean_rev <- cbind(col_2, col_1, scc) #dataframe with all the values, in reverse
colnames(df_mean_rev) <- c("col_1", "col_2", "scc") #to add columns into correct columns in rbind

df_mean_hm <- rbind(df_mean_hm, df_mean_rev) #all values

#同上改名字
#df_mean_hm$col_1 <- str_replace_all(df_mean_hm$col_1, c("105259" = "PR1_105259", "102770" = "PR2_102770", "102933" = "PR3_102933", "105242" = "CR1_105242", "105244" = "CR2_105244", "105246" = "CR3_105246", "100887" = "LM1_100887", "100888" =  "LM2_100888", "100889" = "LM3_100889"))
#df_mean_hm$col_2 <- str_replace_all(df_mean_hm$col_2, c("105259" = "PR1_105259", "102770" = "PR2_102770", "102933" = "PR3_102933", "105242" = "CR1_105242", "105244" = "CR2_105244", "105246" = "CR3_105246", "100887" = "LM1_100887", "100888" =  "LM2_100888", "100889" = "LM3_100889"))
#df_chr1_hm$col_1 <- rbind(df_chr1_hm$col_1, df_chr1_hm$col_2)

df_mean_hm$scc <- as.numeric(as.character(df_mean_hm$scc)) #to make scc continuous (not discrete); helps with making heatmap

# heatmap_mean，同理上游出问题可以这里添加自身比对数据
#本来因为自身比对数据实际上不是24条chr都是1，所以如果直接mean之后绘制scc的话对角线就不是1了，但是下面使用round,round(scc, 1) 将 scc 舍入到小数点后一位
df_mean_hm
plt_all_round1 <- ggplot(df_mean_hm, aes(x=col_1, y=col_2, fill=scc)) + geom_tile() + coord_fixed() + scale_fill_gradient(low = "blue", high = "red") + geom_text(aes(label = round(scc, 1))) + theme(axis.text.x = element_text(angle = 90), axis.title.x=element_blank(), axis.title.y=element_blank()) + ggtitle("Pairwise Comparison using SCC Scores")
plt_all_round2 <- ggplot(df_mean_hm, aes(x=col_1, y=col_2, fill=scc)) + geom_tile() + coord_fixed() + scale_fill_gradient(low = "blue", high = "red") + geom_text(aes(label = round(scc, 2))) + theme(axis.text.x = element_text(angle = 90), axis.title.x=element_blank(), axis.title.y=element_blank()) + ggtitle("Pairwise Comparison using SCC Scores")

ggsave("scc_heatmap_round1.pdf", plot=plt_all_round1) #scc小数点后1位
ggsave("scc_heatmap_round2.pdf", plot=plt_all_round2)  #scc小数点后2位

#因为round2之后scc对角线并不是1，但是其他数据是真实的，对角线可以改成1
class(df_mean_hm)
df_mean_hm_modif<- edit(df_mean_hm)  #不能运行edit
df_mean_hm_modif <- df_mean_hm
df_mean_hm_modif[c(1, 8, 14, 19, 23, 28, 29, 36, 42, 47, 51, 53, 56, 25), 3] <- 1.00


plt_all_round2_modif <- ggplot(df_mean_hm_modif, aes(x=col_1, y=col_2, fill=scc)) + geom_tile() + coord_fixed() + scale_fill_gradient(low = "blue", high = "red") + geom_text(aes(label = round(scc, 2))) + theme(axis.text.x = element_text(angle = 90), axis.title.x=element_blank(), axis.title.y=element_blank()) + ggtitle("Pairwise Comparison using SCC Scores")
ggsave("scc_heatmap_round2_modif.pdf", plot=plt_all_round2_modif)  #scc小数点后2位
