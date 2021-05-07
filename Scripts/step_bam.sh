###########################################################step to prepare bam.sh
dir=$1
bed=$2
nohup Rscript ~/eRNA/script/step1_filter.R ${dir} &
wait;
export BAM_FILE= ${dir}/data/*.bam
# Save the header lines
nohup ~/anaconda2/envs/r-env/bin/samtools view -@ 16 -H $BAM_FILE > ${dir}/data/SAM_header &
##############去掉pca的outliers
######### run.sh
cat ${dir}/data/cluster.txt|awk '{print "CB:Z:"$1}'>${dir}//data/f_barcode.txt

sed -i '1d' ${dir}/data/f_barcode.txt

nohup ~/anaconda2/envs/r-env/bin/samtools view -@ 48 ${dir}/data/*.bam | LC_ALL=C grep -F -f ${dir}/data/f_barcode.txt|grep 'UB:Z' > ${dir}/data/f_filtered_SAM_body &

wait;
nohup cat ${dir}/data/SAM_header ${dir}/data/f_filtered_SAM_body > ${dir}/data/f_filtered.sam &
wait;
# Convert filtered.sam to BAM format
nohup ~/anaconda2/envs/r-env/bin/samtools view -u -@ 48 -b -S ${dir}/data/f_filtered.sam -h -o ${dir}/data/f_filtered.bam &
wait;
nohup ~/anaconda2/envs/r-env/bin/samtools index -@ 48 ${dir}/data/f_filtered.bam &
wait;
#nohup ~/anaconda2/envs/r-env/bin/samtools sort -@ 48 -t CB f_filtered.bam -o sorted_tags.bam -T /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/tmp &
ulimit -n 20000

nohup ~/anaconda2/envs/r-env/bin/bedtools multicov -bams ${dir}/data/f_filtered.bam -bed ${bed} -s >${dir}/data/f_minus_dELS_reads.txt &
nohup ~/anaconda2/envs/r-env/bin/bedtools multicov -bams ${dir}/data/f_filtered.bam -bed ${bed} -S >${dir}/data/f_plus_dELS_reads.txt &
nohup ~/anaconda2/envs/r-env/bin/bedtools multicov -bams ${dir}/data/f_filtered.bam -bed ${bed} >${dir}/data/f_dELS_reads.txt &
mkdir -p ${dir}/data//tmp
nohup ~/anaconda2/envs/r-env/bin/samtools sort -@ 48 -t CB ${dir}/data/f_filtered.bam -o ${dir}/data/sorted_tags.bam -T ${dir}/data//tmp &
wait;
mkdir -p ${dir}/data/split
cd ${dir}/data/splits/
################split bam
nohup ~/anaconda2/envs/r-env/bin/python ~/eRNA/script/split.py ${dir}/data/sorted_tags.bam ${dir}/data/splits/ ${dir}/data/f_barcode.txt &
