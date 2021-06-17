########################################signal_eRNA.sh
dir=$1
type=$2

source activate r-env
nohup Rscript ~/eRNA/script/cut.R ${dir}/data/f_${type}_reads.txt &

wait;
#########cut dot

	d=$(cat ${dir}/data/f_${type}_reads.txt_dot.txt|awk '{print $0}')
	awk '{if($7 >n)  print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$4}' n=$d ${dir}/data/f_${type}_reads.txt >${dir}/data/f_${type}_reads.bed

wait;


mkdir -p ${dir}/data/reads_${type}

n1=$(cat ${dir}/data/f_barcode.txt|wc -l)
n=$(cat ${dir}/data/f_barcode.txt ${dir}/data/f_barcode.txt|wc -l)

rm ${dir}/data/reads.sh
echo "source activate r-env">>${dir}/data/reads.sh
for q in $(seq $n1)
do
	echo "nohup ~/anaconda2/envs/r-env/bin/samtools index -@ 48 ${dir}/data/splits/CB_$q.bam &" >>${dir}/data/reads.sh
done

rm ${dir}/data/reads_${type}/reads.sh

echo "source activate r-env">>${dir}/data/reads_${type}/reads.sh

for q in $(seq $n1)
do
	echo "nohup bedtools multicov -bams ${dir}/data/splits/CB_$q.bam -bed ${dir}/data/f_${type}_reads.bed > ${dir}/data/reads_${type}/reads_$q.txt &" >>${dir}/data/reads_${type}/reads.sh

done

wait;


for i in $(seq 1 5 $n)
do
	sed -i ''"$i"'i wait;' ${dir}/data/reads.sh
done


for i in $(seq 1 5 $n)
do
	sed -i ''"$i"'i wait;' ${dir}/data/reads_${type}/reads.sh
done



wait;
echo "wait;">>${dir}/data/reads_${type}/reads.sh
echo "wait;">>${dir}/data/reads.sh

ulimit -n 20000
nohup bash ${dir}/data/reads.sh &
wait;
nohup bash ${dir}/data/reads_${type}/reads.sh &
wait;

