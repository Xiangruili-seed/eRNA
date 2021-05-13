########################################signal_eRNA.sh
dir=$1
nohup Rscript ~/eRNA/script/cut.R ${dir}/data/f_dELS_reads.txt &

wait;
#########cut dot

	d=$(cat ${dir}/data/f_dELS_reads.txt_dot.txt|awk '{print $0}')
	awk '{if($7 >n)  print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$6}' n=$d ${dir}/data/f_dELS_reads.txt >${dir}/data/f_dELS_reads.bed

wait;


mkdir -p ${dir}/data/reads

n=$(cat ${dir}/data/f_barcode.txt|wc -l)

for q in $(seq $n)
do
	~/anaconda2/envs/r-env/bin/samtools index -@ 48 ${dir}/data/splits/CB_$q.bam
done
wait;


for q in $(seq $n)
do
	echo "nohup bedtools multicov -bams ${dir}/data/splits/CB_$q.bam -bed ${dir}/data/f_dELS_reads.bed > ${dir}/data/reads/reads_$q.txt &" >>${dir}/data/reads/reads.sh

done
wait;
for i in $(seq 1 20 $n)
do
	sed -i ''"$i"'i wait;' ${dir}/data/reads/reads.sh
done
wait;
bash ${dir}/data/reads/reads.sh
wait;
#################cut off>0
cd ${dir}/data/reads/
for q in $(seq $n)
do
	cat ${dir}/data/reads/reads_$q.txt|awk '{if($5>0) print NR" ""'$q'"" "$5}'>${dir}/data/reads/matrix_$q.txt

done
