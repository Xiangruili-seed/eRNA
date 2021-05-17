dir=$1
bed=$2
annotation=$3
all=$4

###########################get rna-seq signal
	nohup bwtool agg 300:300 ${bed} ${dir}/${annotation}/data/filtered.bigWig ${dir}/${annotation}/data/fpos.txt &
	nohup bwtool agg 300:300 ${bed} ${dir}/${annotation}/data/filtered_reverse.bigWig ${dir}/${annotation}/data/fpos_plus.txt &
	nohup bwtool agg 300:300 ${bed} ${dir}/${annotation}/data/filtered_forward.bigWig  ${dir}/${annotation}/data/fpos_minus.txt &
wait;
	nohup Rscript ~/eRNA/script/RNA.R ${dir}/${annotation}/data/ &
#################chip signal

ls -R ${dir}/chip/*.bigWig > ${dir}/chip/annotation.txt
n=$(cat ${dir}/chip/annotation.txt|wc -l)



cat ${all}|awk '{print $1"\t"$2"\t"$3}' >${all}.bed
for k in $(seq $n)\t
do
    peakFile=$(cat ${dir}/chip/annotation.txt| awk -F "\t" '{if (NR == '$k') print $1}')
    bwtool agg 500:500 ${bed} $peakFile $peakFile.pos.txt
    bwtool agg 500:500 ${all}.bed $peakFile $peakFile.txt

done
wait;


for k in $(seq $n)
do
    peakFile=$(cat ${dir}/chip/annotation.txt| awk -F "\t" '{if (NR == '$k') print $1}')
    nohup Rscript ~/eRNA/script/chip_t.R $peakFile pos &


done

