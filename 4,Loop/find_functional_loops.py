# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 10:05:30 2020

@author: Kara
"""

import sys

#reads a bed file of promoters into a dictionary.
def readbed(filename):   #读取bed文件中的3列，用于读取后面的P以及E的bed文件
    i=0
    dictname = {}
    with open(filename, encoding="utf8", errors='ignore') as bedfile:
        for line in bedfile:
            sline = line.rstrip().split()
            chrom = sline[0]
            start=int(sline[1])
            end=int(sline[2])
            dictname[i]=[chrom,start,end]
            i+=1
    return dictname
    
    
    
        
def find_cre(chromosome,position,cre_dict,threshold,side,bidirectional):#mess around with minval and maxval根据给定的染色体、位置和搜索条件，在指定的调控元件字典中查找匹配的调控元件，并返回是否找到以及找到的调控元件信息，重要的是threshold这一个参数，是用于设置查找范围的距离参数，主要是左右两侧，大概和E的长度有关
    foundstuff={}
    if side=="right":
        maxval=position+threshold
        if bidirectional ==True:
            minval=position - threshold
        else:
            minval=position
    else:
        if bidirectional == True:
            maxval=position + threshold
        else: 
            maxval=position
        minval=position-threshold
    foundone=False
    i=0
    for key in cre_dict:
        cre=cre_dict[key]
        cre_chrom=cre[0]
        if cre_chrom==chromosome:
            cre_start=cre[1]
            if cre_start >= minval:
                cre_end = cre_dict[key][2]
                #print("%s\t%s\t%s"%(cre_chrom,cre_start,cre_end))
                if cre_start<=maxval:
                    foundone=True
                    foundstuff[i]=[cre_chrom,cre_start,cre_end]
                    i+=1
    return foundone,foundstuff
                    
#在基因组学和调控区域分析中，阈值参数 threshold 所代表的意义是对于给定的基因组位置，要查找在该位置附近一定范围内的调控元件（启动子或增强子）。具体来说：
#threshold 表示了搜索的范围，以基因组上的位置为中心，向左右两侧延伸 threshold 的距离，形成一个搜索窗口。
#如果在这个搜索窗口内找到了符合条件的调控元件（位置在搜索窗口内），则会被记录下来，否则将被忽略


#python find_functional_loops.py "/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/SE/enhancer/sortbed/HMEC_SE_sort.bed"   "/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result/TNBC_preprocessing_any/loops_HMEC_all_10000_summit.bed"  HMEC_SEPC.txt 
#################MAIN SCRIPT###################################################

       
promoter_file="/mnt/disk4/haitao/bysj_seu/ref_genome/new_hg19/hg19_promoter_sort.bed"
promoters = readbed(promoter_file)
enhancer_file = sys.argv[0]   #可以另外写一个bash脚本在for循环中输入第一个参数
#"/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/SE/enhancer/sortbed/HMEC_SE_sort.bed"
loop_loc_file = sys.argv[1] #查看下面的操作，实际上这里获取的bed文件只有3列，看起来是anchor文件，但是看最下面左E右P左P右E的脚本处理，这里应该提供的是一个gr对象，确实是左右anchor范围的loop的bed文件
#threshold = int(sys.argv[2])
EPCoutput_file =sys.argv[2]    #"human_all_enhancer_promoter_candidates_10kb.txt"
bidirectional = True   #所以一开始就是设置为双向寻找
#looplist = []#list of loopnames for loops of interest
#with open(loopfile) as f:
#    for line in f:
#        sline = line.rstrip().split()
#        looplist.append(sline[0])
enhancers=readbed(enhancer_file)
loopdict={}
with open(loop_loc_file, encoding="utf8", errors='ignore') as f:
    for line in f:
        sline = line.rstrip().split()
        loopname = sline[3]
        ch = sline[0]
        start = int(sline[1])
        end = int(sline[2])
        loopdict[loopname] = [ch,start,end,loopname]

candidates_dict={}
cres_dict={}
looplengths=[]
for key in loopdict:
    right_promoters={}
    right_enhancers={}
    left_promoters={}
    left_enhancers={}
    left_anchor=int(loopdict[key][1])
    right_anchor=int(loopdict[key][2])
    chromosome=loopdict[key][0]
    loop_length=right_anchor-left_anchor
    looplengths.append(loop_length)
        #check for promoters within 10kb of right anchor，P区域使用2kb的数据就设置为7kb左右
    right_promoter,right_promoters=find_cre(chromosome,right_anchor,promoters,7000,"right",bidirectional)
    if right_promoter==True:
        #if found, check for enhancers within 10kb of left anchor,还是使用12.5kb的数据，那就设置为17.5kb大小
         left_enhancer,left_enhancers=find_cre(chromosome,left_anchor,enhancers,17500,"left",bidirectional)
         if left_enhancer==True:
            #if found, add loop to dictionary of candidates
            candidates_dict[key]=loopdict[key]
    #check for promoters within 10kb of left anchor
    left_promoter,left_promoters=find_cre(chromosome,left_anchor,promoters,7000,"left",bidirectional)
    if left_promoter==True:
        #if found, check for enhancers within 10kb of right anchor
        right_enhancer,right_enhancers=find_cre(chromosome,right_anchor,enhancers,17500,"right",bidirectional)
        if right_enhancer==True:
            #if found, add loop to dictionary of candidates
            if key not in candidates_dict.keys():
                candidates_dict[key]=loopdict[key]
    cres_dict[key]=[right_promoters,left_enhancers,left_promoters,right_enhancers]
                    

with open(EPCoutput_file,"w") as nf:
    nf.write("loop_name\tchromosome\tstart\tend\t\n")
    for key in candidates_dict:
        nf.write("%s\t%s\t%s\t%s\t\n"%(key,candidates_dict[key][0],candidates_dict[key][1],candidates_dict[key][2]))

































