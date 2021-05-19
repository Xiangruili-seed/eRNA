dir=$1


#################erna
n1=$(cat ${dir}/data/f_barcode.txt|wc -l)

  #########eRNA
#dir=~/eRNA/12878/outs/filtered_feature_bc_matrix 
for q in $(seq $n1)
do
	cat ${dir}/data/reads/reads_$q.txt|awk '{if($7>0) print NR"\t""'$q'""\t"$7}'>${dir}/data/reads/matrix_$q.txt

done

mkdir -p ${dir}/pos
mkdir -p ${dir}/pos_add
gzip -d -c ${dir}/matrix.mtx.gz|head -2 > ${dir}/pos/matrix.mtx
cell=$(cat ${dir}/data/f_barcode.txt|wc -l)
eRNA=$(cat ${dir}/data/f_dELS_reads.bed|wc -l)
cat ${dir}/data/reads/matrix_* >> ${dir}/pos/mtx
line=$(cat ${dir}/pos/mtx|wc -l)
echo -e $eRNA"\t"$cell"\t"$line >>${dir}/pos/matrix.mtx
cat ${dir}/pos/mtx>>${dir}/pos/matrix.mtx
cp ${dir}/data//f_barcode.txt ${dir}//pos/barcodes.tsv
cat ${dir}/data/f_dELS_reads.bed|awk '{print $4"\t"$1"_"$4"\tGene"}' >${dir}//pos/genes.tsv
#############eRNA_GENES
features=$(gzip -d -c ${dir}/features.tsv.gz|wc -l)
gzip -d -c ${dir}/barcodes.tsv.gz|awk '{print NR" CB:Z:"$0}'>${dir}/pos_add/all_barcode.txt
cat ${dir}/data//f_barcode.txt|awk '{print NR"\t"$0}'>${dir}/pos_add//f_barcode.txt
join -1 2 -2 2 -o 1.1 2.1 ${dir}/pos_add/all_barcode.txt ${dir}/pos_add/f_barcode.txt|sort -k 2 >${dir}/pos_add/id.txt
tail -n +4 ${dir}/pos/matrix.mtx|sort -k 2 >${dir}/pos_add/m.txt
join -1 2 -2 2 -o 1.1 2.1 1.3 ${dir}/pos_add/m.txt ${dir}/pos_add/id.txt|awk '{print $1+features"\t"$2"\t"$3}' features=$features >${dir}/pos_add/t.txt
gzip -d -c ${dir}/matrix.mtx.gz|head -2 > ${dir}/pos_add/matrix.mtx
gzip -d -c ${dir}/matrix.mtx.gz|sed -n '3p'|awk '{print $1+eRNA"\t"$2"\t"$3+line}' eRNA=$eRNA line=$line>>${dir}/pos_add/matrix.mtx
gzip -d -c ${dir}/matrix.mtx.gz|tail -n +4 |awk '{print $1"\t"$2"\t"$3}'>${dir}/pos_add/a.txt
cat ${dir}/pos_add/a.txt ${dir}/pos_add/t.txt>>${dir}/pos_add/matrix.mtx
cp ${dir}/barcodes.tsv.gz ${dir}/pos_add/barcodes.tsv.gz
gzip -d -c ${dir}/features.tsv.gz|awk '{print $1"\t"$2"\t"$3}' >${dir}/pos_add/gene.tsv
cat ${dir}/pos_add/gene.tsv ${dir}/pos/genes.tsv>${dir}/pos_add/genes.tsv

##############gzip
gzip ${dir}/pos_add/genes.tsv
mv ${dir}/pos_add/genes.tsv.gz ${dir}/pos_add/features.tsv.gz

gzip ${dir}/pos/genes.tsv
mv ${dir}/pos/genes.tsv.gz ${dir}/pos/features.tsv.gz
gzip ${dir}/pos_add/matrix.mtx
gzip ${dir}/pos/barcodes.tsv
gzip ${dir}/pos/matrix.mtx

