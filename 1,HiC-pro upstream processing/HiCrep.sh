#!/bin/bashs

export HDF5_USE_FILE_LOCKING='FALSE'

DIRIN=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script5    #mcool文件的工作目录
DIROUT=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script5/scc  #指定输出文件文件夹

array=(BT549_rep1 BT549_rep2 HCC70_rep1 HCC70_rep2 MB231 HMEC_rep1 HMEC_rep2)    #包含了一系列字符串的数组，每个字符串都是一个 .mcool 文件的文件名，不包含 .mcool 扩展名，到时候修改成对应文件名字

total=$((${#array[@]}-1)) # Minus 1 to match array indices，数组上小标
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

#原来是大于gt，但是后来发现没有自己比对的scc，所以改成了ge大于等于
