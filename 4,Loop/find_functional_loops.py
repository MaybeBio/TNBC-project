# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 10:05:30 2020

@author: Kara
"""

import sys

#reads a bed file of promoters into a dictionary.
def readbed(filename):   #��ȡbed�ļ��е�3�У����ڶ�ȡ�����P�Լ�E��bed�ļ�
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
    
    
    
        
def find_cre(chromosome,position,cre_dict,threshold,side,bidirectional):#mess around with minval and maxval���ݸ�����Ⱦɫ�塢λ�ú�������������ָ���ĵ���Ԫ���ֵ��в���ƥ��ĵ���Ԫ�����������Ƿ��ҵ��Լ��ҵ��ĵ���Ԫ����Ϣ����Ҫ����threshold��һ�����������������ò��ҷ�Χ�ľ����������Ҫ���������࣬��ź�E�ĳ����й�
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
                    
#�ڻ�����ѧ�͵�����������У���ֵ���� threshold ������������Ƕ��ڸ����Ļ�����λ�ã�Ҫ�����ڸ�λ�ø���һ����Χ�ڵĵ���Ԫ���������ӻ���ǿ�ӣ���������˵��
#threshold ��ʾ�������ķ�Χ���Ի������ϵ�λ��Ϊ���ģ��������������� threshold �ľ��룬�γ�һ���������ڡ�
#�������������������ҵ��˷��������ĵ���Ԫ����λ�������������ڣ�����ᱻ��¼���������򽫱�����


#python find_functional_loops.py "/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/SE/enhancer/sortbed/HMEC_SE_sort.bed"   "/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result/TNBC_preprocessing_any/loops_HMEC_all_10000_summit.bed"  HMEC_SEPC.txt 
#################MAIN SCRIPT###################################################

       
promoter_file="/mnt/disk4/haitao/bysj_seu/ref_genome/new_hg19/hg19_promoter_sort.bed"
promoters = readbed(promoter_file)
enhancer_file = sys.argv[0]   #��������дһ��bash�ű���forѭ���������һ������
#"/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/SE/enhancer/sortbed/HMEC_SE_sort.bed"
loop_loc_file = sys.argv[1] #�鿴����Ĳ�����ʵ���������ȡ��bed�ļ�ֻ��3�У���������anchor�ļ������ǿ���������E��P��P��E�Ľű���������Ӧ���ṩ����һ��gr����ȷʵ������anchor��Χ��loop��bed�ļ�
#threshold = int(sys.argv[2])
EPCoutput_file =sys.argv[2]    #"human_all_enhancer_promoter_candidates_10kb.txt"
bidirectional = True   #����һ��ʼ��������Ϊ˫��Ѱ��
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
        #check for promoters within 10kb of right anchor��P����ʹ��2kb�����ݾ�����Ϊ7kb����
    right_promoter,right_promoters=find_cre(chromosome,right_anchor,promoters,7000,"right",bidirectional)
    if right_promoter==True:
        #if found, check for enhancers within 10kb of left anchor,����ʹ��12.5kb�����ݣ��Ǿ�����Ϊ17.5kb��С
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

































