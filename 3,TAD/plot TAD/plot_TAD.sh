#!/bin/bash

#mamba activate HiCExplorer

region=$2    #ָ�����Ƶ�����������chrX:6800000-8500000����ʽ��bed�ļ��������Լ���һЩ����λ��Ҳ����
ini=$1   #ָ�����Ƶ�ini�����ļ�
output_prefix=$3  #ָ�����ͼƬ��ǰ׺

#common�Ƚϣ���ȫ��TAD trackͼ���Ƚϣ����漰������TAD��������show��ʹ��TADs-common.ini�����ڲ����TAD�ȽϿ�����������һ��ini�ļ�
#������ini�н����е�TAD��ͳһ������Ҳ����д��ѭ����Ȼ�����ò�ͬ��ini����ʱ����ѭ����ͳһ����ͬһ������

hicPlotTADs --tracks  ${ini} --region ${region} -o ${output_prefix}.pdf