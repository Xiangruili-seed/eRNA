dir=$1
bed=$2
BAM_FILE=$3

nohup Rscript ~/eRNA/script/step1_filter.R ${dir} &
wait;
# Save the header lines20.12.27
source activate r-env
nohup ~/anaconda2/envs/r-env/bin/samtools view -@ 16 -H ${BAM_FILE} > ${dir}/data/SAM_header &
##############去掉pca的outliers
######### run.sh
cat ${dir}/data/cluster.txt|awk '{print "CB:Z:"$1}'>${dir}//data/f_barcode.txt

sed -i '1d' ${dir}/data/f_barcode.txt

nohup ~/anaconda2/envs/r-env/bin/samtools view -@ 48 ${dir}/data/*.bam | LC_ALL=C grep -F -f ${dir}/data/f_barcode.txt|grep 'UB:Z' > ${dir}/data/f_filtered_SAM_body &

wait;
nohup cat ${dir}/data/SAM_header ${dir}/data/f_filtered_SAM_body > ${dir}/data/f_filtered.sam &
wait;
# Convert filtered.sam to BAM format
nohup ~/anaconda2/envs/r-env/bin/samtools view -u -@ 48 -b -S ${dir}/data/f_filtered.sam -h -o ${dir}/data/f_filtered_1.bam &
wait;
nohup samtools view -H ${dir}/data/f_filtered_1.bam | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - ${dir}/data/f_filtered_1.bam > ${dir}/data/f_filtered.bam &

wait;
nohup ~/anaconda2/envs/r-env/bin/samtools index -@ 48 ${dir}/data/f_filtered.bam &
wait;
#nohup ~/anaconda2/envs/r-env/bin/samtools sort -@ 48 -t CB f_filtered.bam -o sorted_tags.bam -T /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/tmp &
ulimit -n 20000

nohup bedtools multicov -bams ${dir}/data/f_filtered.bam -bed ${bed} -s >${dir}/data/f_minus_dELS_reads.txt &
nohup bedtools multicov -bams ${dir}/data/f_filtered.bam -bed ${bed} -S >${dir}/data/f_plus_dELS_reads.txt &
nohup bedtools multicov -bams ${dir}/data/f_filtered.bam -bed ${bed} >${dir}/data/f_dELS_reads.txt &
mkdir -p ${dir}/data//tmp
nohup ~/anaconda2/envs/r-env/bin/samtools sort -@ 48 -t CB ${dir}/data/f_filtered.bam -o ${dir}/data/sorted_tags.bam -T ${dir}/data//tmp &
wait;
mkdir -p ${dir}/data/splits
cd ${dir}/data/splits/
################split bam
nohup ~/anaconda2/envs/r-env/bin/python ~/eRNA/script/split.py ${dir}/data/sorted_tags.bam ${dir}/data/splits/ ${dir}/data/f_barcode.txt &
nohup bamCoverage -bs 10 -b ${dir}/data/f_filtered.bam --filterRNAstrand forward -o ${dir}/data/filtered_forward.bigWig & 负
nohup bamCoverage -bs 10 -b ${dir}/data/f_filtered.bam -o ${dir}/data/filtered.bigWig &
nohup bamCoverage -bs 10 -b ${dir}/data/f_filtered.bam --filterRNAstrand reverse -o ${dir}/data/filtered_reverse.bigWig &正
