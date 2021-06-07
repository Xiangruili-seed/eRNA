dir=$1
type=$2


#################erna
n1=$(cat ${dir}/data/f_barcode.txt|wc -l)

  #########eRNA
#dir=~/eRNA/12878/outs/filtered_feature_bc_matrix 
for q in $(seq $n1)
do
	cat ${dir}/data/reads_${type}/reads_$q.txt|awk '{if($7>0) print NR"\t""'$q'""\t"$7}'>${dir}/data/reads_${type}/matrix_$q.txt

done

mkdir -p ${dir}/pos_${type}
mkdir -p ${dir}/pos_add_${type}
gzip -d -c ${dir}/matrix.mtx.gz|head -2 > ${dir}/pos_${type}/matrix.mtx
cell=$(cat ${dir}/data/f_barcode.txt|wc -l)
eRNA=$(cat ${dir}/data/f_${type}_reads.bed|wc -l)
cat ${dir}/data/reads_${type}/matrix_* >> ${dir}/pos_${type}/mtx
line=$(cat ${dir}/pos_${type}/mtx|wc -l)
echo -e $eRNA"\t"$cell"\t"$line >>${dir}/pos_${type}/matrix.mtx
cat ${dir}/pos_${type}/mtx>>${dir}/pos_${type}/matrix.mtx
cp ${dir}/data//f_barcode.txt ${dir}//pos_${type}/barcodes.tsv
cat ${dir}/data/f_${type}_reads.bed|awk '{print $4"\t"$1"_"$4"\tGene"}' >${dir}//pos_${type}/genes.tsv
#############eRNA_GENES
features=$(gzip -d -c ${dir}/features.tsv.gz|wc -l)
gzip -d -c ${dir}/barcodes.tsv.gz|awk '{print NR" CB:Z:"$0}'>${dir}/pos_add_${type}/all_barcode.txt
cat ${dir}/data//f_barcode.txt|awk '{print NR"\t"$0}'>${dir}/pos_add_${type}//f_barcode.txt
join -1 2 -2 2 -o 1.1 2.1 ${dir}/pos_add_${type}/all_barcode.txt ${dir}/pos_add_${type}/f_barcode.txt|sort -k 2 >${dir}/pos_add_${type}/id.txt
tail -n +4 ${dir}/pos_${type}/matrix.mtx|sort -k 2 >${dir}/pos_add_${type}/m.txt
join -1 2 -2 2 -o 1.1 2.1 1.3 ${dir}/pos_add_${type}/m.txt ${dir}/pos_add_${type}/id.txt|awk '{print $1+features"\t"$2"\t"$3}' features=$features >${dir}/pos_ad_${type}d/t.txt
gzip -d -c ${dir}/matrix.mtx.gz|head -2 > ${dir}/pos_add_${type}/matrix.mtx
gzip -d -c ${dir}/matrix.mtx.gz|sed -n '3p'|awk '{print $1+eRNA"\t"$2"\t"$3+line}' eRNA=$eRNA line=$line>>${dir}/pos_add_${type}/matrix.mtx
gzip -d -c ${dir}/matrix.mtx.gz|tail -n +4 |awk '{print $1"\t"$2"\t"$3}'>${dir}/pos_add_${type}/a.txt
cat ${dir}/pos_add_${type}/a.txt ${dir}/pos_add_${type}/t.txt>>${dir}/pos_add_${type}/matrix.mtx
cp ${dir}/barcodes.tsv.gz ${dir}/pos_add_${type}/barcodes.tsv.gz
gzip -d -c ${dir}/features.tsv.gz|awk '{print $1"\t"$2"\t"$3}' >${dir}/pos_add_${type}/gene.tsv
cat ${dir}/pos_add_${type}/gene.tsv ${dir}/pos_${type}/genes.tsv>${dir}/pos_add_${type}/genes.tsv

##############gzip
gzip ${dir}/pos_add_${type}/genes.tsv
mv ${dir}/pos_add_${type}/genes.tsv.gz ${dir}/pos_add_${type}/features.tsv.gz

gzip ${dir}/pos_${type}/genes.tsv
mv ${dir}/pos_${type}/genes.tsv.gz ${dir}/pos_${type}/features.tsv.gz
gzip ${dir}/pos_add_${type}/matrix.mtx
gzip ${dir}/pos_${type}/barcodes.tsv
gzip ${dir}/pos_${type}/matrix.mtx

