########################################signal_eRNA.sh
dir=$1
source activate r-env
nohup Rscript ~/eRNA/script/cut.R ${dir}/data/f_dELS_reads.txt &

wait;
#########cut dot

	d=$(cat ${dir}/data/f_dELS_reads.txt_dot.txt|awk '{print $0}')
	awk '{if($7 >n)  print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$6}' n=$d ${dir}/data/f_dELS_reads.txt >${dir}/data/f_dELS_reads.bed

wait;


mkdir -p ${dir}/data/reads

n1=$(cat ${dir}/data/f_barcode.txt|wc -l)
rm ${dir}/data/reads.sh
echo "source activate r-env">>${dir}/data/reads.sh
for q in $(seq $n1)
do
	echo "nohup ~/anaconda2/envs/r-env/bin/samtools index -@ 48 ${dir}/data/splits/CB_$q.bam &" >>${dir}/data/reads.sh
done
wait
for i in $(seq 1 30 $n1)
do
	sed -i ''"$i"'i wait;' ${dir}/data/reads.sh
done
wait

rm ${dir}/data/reads/reads.sh

echo "source activate r-env">>${dir}/data/reads/reads.sh

for q in $(seq $n1)
do
	echo "nohup bedtools multicov -bams ${dir}/data/splits/CB_$q.bam -bed ${dir}/data/f_dELS_reads.bed > ${dir}/data/reads/reads_$q.txt &" >>${dir}/data/reads/reads.sh

done

wait;
#################cut off>0

for i in $(seq 1 20 $n1)
do
	sed -i ''"$i"'i wait;' ${dir}/data/reads/reads.sh
done
wait;
ulimit -n 20000
nohup bash ${dir}/data/reads.sh &
wait;
nohup bash ${dir}/data/reads/reads.sh &
wait;

