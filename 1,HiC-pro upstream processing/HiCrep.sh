#!/bin/bashs

export HDF5_USE_FILE_LOCKING='FALSE'

DIRIN=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script5    #mcool�ļ��Ĺ���Ŀ¼
DIROUT=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script5/scc  #ָ������ļ��ļ���

array=(BT549_rep1 BT549_rep2 HCC70_rep1 HCC70_rep2 MB231 HMEC_rep1 HMEC_rep2)    #������һϵ���ַ��������飬ÿ���ַ�������һ�� .mcool �ļ����ļ����������� .mcool ��չ������ʱ���޸ĳɶ�Ӧ�ļ�����

total=$((${#array[@]}-1)) # Minus 1 to match array indices��������С��
# total=
for i in $(seq 0 $total); do
    for j in $(seq 0 $total); do
        if [ $j -ge $i ]    
        then
            FILEIN1=$DIRIN"/"${array[$i]}".mcool"
            FILEIN2=$DIRIN"/"${array[$j]}".mcool"
            FILEOUT=$DIROUT"/"${array[$i]}"_"${array[$j]}".txt"
            hicrep  --binSize 10000 --h 10 --dBPMax 5000000 $FILEIN1 $FILEIN2 $FILEOUT
        fi
    done
done

#ԭ���Ǵ���gt�����Ǻ�������û���Լ��ȶԵ�scc�����Ըĳ���ge���ڵ���
