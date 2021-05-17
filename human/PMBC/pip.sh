cat /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/filtered_feature_bc_matrix/barcodes.tsv|awk '{ print "CB:Z:"$0}'>barcode.txt
cp /data/projects/encode/Registry/V2/GRCh38/GRCh38-ccREs.dELS.bed GRCh38-ccREs.dELS.bed 

cat /data/tusers/lixiangr/eRNA/refer/human/hg38/gencode.v35.annotation.gtf|awk -F '[\t|;]' '$3=="exon"{print $1"\t"$4"\t"$5"\t"$10}'>exon.bed
bedtools intersect -a GRCh38-ccREs.dELS.bed -b exon.bed -v>filtered_dELS.bed
cat filtered_dELS.bed|awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t-"}'>GRCh38-ccREs.dELS_minus.bed
##########################
export BAM_FILE='/data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/SC3_v3_NextGem_SI_PBMC_10K_possorted_genome_bam.bam'
# Save the header lines
nohup ~/anaconda2/envs/r-env/bin/samtools view -@ 16 -H $BAM_FILE > SAM_header &



##############去掉pca的outliers
######### run.sh
cat ~/eRNA/pbmcs/cells/reads/cluster.txt|awk '{print "CB:Z:"$1}'>/data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/f_barcode.txt

sed -i '1d' /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/f_barcode.txt

nohup ~/anaconda2/envs/r-env/bin/samtools view -@ 48 /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/SC3_v3_NextGem_SI_PBMC_10K_possorted_genome_bam.bam | LC_ALL=C grep -F -f /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/f_barcode.txt|grep 'UB:Z' > f_filtered_SAM_body &

wait;
nohup cat SAM_header f_filtered_SAM_body > f_filtered.sam &
wait;
# Convert filtered.sam to BAM format
nohup ~/anaconda2/envs/r-env/bin/samtools view -u -@ 48 -b -S f_filtered.sam -h -o f_filtered.bam &
wait;
nohup ~/anaconda2/envs/r-env/bin/samtools index -@ 48 f_filtered.bam &
wait;
#nohup ~/anaconda2/envs/r-env/bin/samtools sort -@ 48 -t CB f_filtered.bam -o sorted_tags.bam -T /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/tmp &
ulimit -n 20000

nohup bedtools multicov -bams f_filtered.bam -bed GRCh38-ccREs.dELS_minus.bed -s >f_minus_dELS_reads.txt &
nohup bedtools multicov -bams f_filtered.bam -bed GRCh38-ccREs.dELS_minus.bed -S >f_plus_dELS_reads.txt &
nohup bedtools multicov -bams f_filtered.bam -bed GRCh38-ccREs.dELS_minus.bed >f_dELS_reads.txt &

nohup ~/anaconda2/envs/r-env/bin/samtools sort -@ 48 -t CB /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/f_filtered.bam -o ~/eRNA/pbmcs/sorted_tags.bam -T ~/eRNA/pbmcs//tmp &
wait;
cd ~/eRNA/pbmcs/
################split bam
nohup python split.py &
wait;
####################并行 index
ulimit -n 20000
cd ~/eRNA/pbmcs/cell
for q in $(seq 10135)
do
	~/anaconda2/envs/r-env/bin/samtools index -@ 48 CB_$q.bam
done
wait;
for q in $(seq 10135)
do
	samtools view CB_$q.bam|head -1|awk '{for (i=1;i<=NF;i++){if ($i ~/CB/) {print $i}}}'>>~/eRNA/pbmcs/cell/barcode.txt
done

#################################################################################get reads count

cat pos_minus_*.bed|awk '{print $1"\t"$2"\t"$3"\t"$4}'|sort -k 4 -u >pos_new.bed

for q in $(seq 10135)
do
	echo "nohup bedtools multicov -bams ~/eRNA/pbmcs/cell/CB_$q.bam -bed /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/pos_new.bed > ~/eRNA/pbmcs/cell/reads/reads_$q.txt &" >>~/eRNA/pbmcs/cell/reads/reads.sh

done

for i in $(seq 1 50 10135)
do
	sed -i ''"$i"'i wait;' ~/eRNA/pbmcs/cell/reads/reads.sh
done

nohup bash ~/eRNA/pbmcs/cell/reads/reads.sh &

#################cut off>0
cd ~/eRNA/pbmcs/cell/reads/
for q in $(seq 10135)
do
	cat reads_$q.txt|awk '{if($5>0) print NR" ""'$q'"" "$5}'>matrix_$q.txt

done

#################erna
vi ~/eRNA/pbmcs/cell/reads/pos/matrix.mtx
%%MatrixMarket matrix coordinate integer general
%metadata_json: {"software_version": "cellranger-4.0.0", "format_version": 2}
690 10135 587486


cat ~/eRNA/pbmcs/cell/reads/matrix_* >> ~/eRNA/pbmcs/cell/reads/pos/matrix.mtx
cp ~/eRNA/pbmcs/cell/barcode.txt ~/eRNA/pbmcs/cell/reads/pos/barcodes.tsv

cat /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/pos_new.bed|awk '{print $4"\t"$1"_"$4}' >~/eRNA/pbmcs/cell/reads/pos/genes.tsv




################################gene+erna
cat /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/barcode.txt|awk '{print NR"\t"$0}'>~/eRNA/pbmcs/all_barcode.txt

cat /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/f_barcode.txt|awk '{print NR"\t"$0}'>~/eRNA/pbmcs/f_barcode.txt
join -1 2 -2 2 -o 1.1 2.1 ~/eRNA/pbmcs/all_barcode.txt ~/eRNA/pbmcs/f_barcode.txt|sort -k 2 >~/eRNA/pbmcs/cell/reads/pos_add/id.txt


tail -n +4 ~/eRNA/pbmcs/cell/reads/pos/matrix.mtx|sort -k 2 >~/eRNA/pbmcs/cell/reads/pos_add/m.txt


join -1 2 -2 2 -o 1.1 2.1 1.3 ~/eRNA/pbmcs/cell/reads/pos_add/m.txt ~/eRNA/pbmcs/cell/reads/pos_add/id.txt|awk '{print $1+36601"\t"$2"\t"$3}'>~/eRNA/pbmcs/cell/reads/pos_add/t.txt

vi ~/eRNA/pbmcs/cell/reads/pos_add/matrix.mtx

%%MatrixMarket matrix coordinate integer general
%metadata_json: {"software_version": "cellranger-4.0.0", "format_version": 2}
37291 10985 29052410
tail -n +4 /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/filtered_feature_bc_matrix/hg38/matrix.mtx>~/eRNA/pbmcs/cell/reads/pos_add/a.txt
cat ~/eRNA/pbmcs/cell/reads/pos_add/a.txt ~/eRNA/pbmcs/cell/reads/pos_add/t.txt>>~/eRNA/pbmcs/cell/reads/pos_add/matrix.mtx

cp /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/barcode.txt barcodes.tsv

cat /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/filtered_feature_bc_matrix/hg38/genes.tsv ~/eRNA/pbmcs/cell/reads/pos/genes.tsv>~/eRNA/pbmcs/cell/reads/pos_add/genes.tsv
######
Rscript PMBC_pos.R

cat ~/eRNA/pbmcs/pos/top_0.05.txt|grep EH38D|sort -k 8 -u>~/eRNA/pbmcs/pos/top_0.05_erna.txt
cat ~/eRNA/pbmcs/pos/top_10.txt|grep EH38D|sort -k 8 -u>~/eRNA/pbmcs/pos/top_10_erna.txt
cat ~/eRNA/pbmcs/pos/top_100.txt|grep EH38D|sort -k 8 -u>~/eRNA/pbmcs/pos/top_100_erna.txt
cd ~/eRNA/refdata-gex-GRCh38-2020-A/genes
cat genes.gtf|awk '{if($3 =="gene") print $1"\t"$4"\t"$5"\t"$16}'|grep 'chr'>gene.bed

cat /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/GRCh38-ccREs.dELS.bed|sort -k 4 >ccre.bed
cat ~/eRNA/pbmcs/pos/top_0.05.txt|grep EH38D|sort -k 8 -k3nr|sort -k 8 -u|sort -k7n -k3nr|sort -k7n -u>~/eRNA/pbmcs/pos/top_0.05_erna_top1.txt
cat ~/eRNA/pbmcs/pos/top_100.txt|grep EH38D|sort -k 8 -k3nr|sort -k 8 -u|sort -k7n -k3nr|sort -k7n -u>~/eRNA/pbmcs/pos/top_100_erna_top1.txt


cat ~/eRNA/pbmcs/pos/top_10.txt|grep EH38D|sort -k 8 -k3nr|sort -k 8 -u|sort -k7n -k3nr|sort -k7n -u>~/eRNA/pbmcs/pos/top_10_erna_top1.txt
cat ~/eRNA/pbmcs/pos/top_0.05_erna.txt|awk '{print $8}'|awk -F "-" '{print $2}'|sort >select.bed

join -1 1 -2 4 select.bed ccre.bed|awk '{print $2"\t"$3-100000"\t"$4+100000"\t"$1"\t"$6}'> select_ccre.bed
 bedtools intersect -a select_ccre.bed -b gene.bed -wo >select_gene.bed
 cat select_gene.bed|awk '{print $1"-"$4"\t"$9}'| sed 's/\"//g'| sed 's/\;//g'|sort -k 2 >gene_ccre.txt
sort -k 8 ~/eRNA/pbmcs/pos/top_0.05.txt>top_0.05.txt
join -1 2 -2 8 -o 1.1 1.2 gene_ccre.txt top_0.05.txt|sort -u >s_gene_id.txt

cat ~/eRNA/pbmcs/pos/top_0.05_cor_cut.txt|awk -F "-" '{print $2}'|sort >select_cut.bed
join -1 1 -2 4 -o 2.1 2.2 2.3 2.4 select_cut.bed ccre.bed> select_ccre_cut.bed
###############
head ~/eRNA/pbmcs/cell/reads/pos//cluster.txt
head ~/eRNA/pbmcs/cells/reads/cell_type.txt
head ~/eRNA/pbmcs/cells/reads/cluster.txt
cat ~/eRNA/pbmcs/cells/reads/cluster.txt|awk '{ print "CB:Z:"$0}'|sort -k 1 >~/eRNA/pbmcs/cell/reads/pos/id

cat cluster.txt|sort -k 1 >~/eRNA/pbmcs/cell/reads/pos/id_sort
cat ~/eRNA/pbmcs/cells/reads/cell_type.txt|sort -k 3 >~/eRNA/pbmcs/cell/reads/pos/id_c_sort
join -1 1 -2 1 -o 1.1 1.2 2.2 ~/eRNA/pbmcs/cell/reads/pos/id_sort ~/eRNA/pbmcs/cell/reads/pos/id|sort -k 3 >s
join -1 3 -2 3 -o 1.1 1.2 1.3 2.2 s ~/eRNA/pbmcs/cell/reads/pos/id_c_sort >o

#################################################################################################










##############

conda install -c conda-forge r-corrplot 

conda create --name r4-base
source activate r4-base
conda install -c conda-forge r-base


conda install -c bioconda bioconductor-biocfilecache 
conda install -c bioconda bioconductor-nebulosa 

conda install -c conda-forge r-patchwork 

conda install -c conda-forge r-seurat

conda install -c conda-forge r-dplyr 



#/data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/filtered_feature_bc_matrix/.RData
#########
~/eRNA/pbmcs/cells/reads/cell_type.txt
~/eRNA/pbmcs/cells/reads/cluster.txt
load("/data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/filtered_feature_bc_matrix/MAESTRO/data/human.immune.CIBERSORT.RData")


cat ~/eRNA/pbmcs/cells/reads/cluster.txt|awk '{if($2==0) print "CB:Z:"$1}'>/data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/CD8Tcells.txt
cat ~/eRNA/pbmcs/cells/reads/cluster.txt|awk '{if ($2==1 ||$2==2||$2==7||$2==8) print "CB:Z:"$1}'>/data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/Monocytes.txt
cat ~/eRNA/pbmcs/cells/reads/cluster.txt|awk '{if($2==3) print "CB:Z:"$1}'>/data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/RestMemCD4Tcells.txt
cat ~/eRNA/pbmcs/cells/reads/cluster.txt|awk '{if ($2==4 ||$2==12) print "CB:Z:"$1}'>/data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/NaiveBcells.txt
cat ~/eRNA/pbmcs/cells/reads/cluster.txt|awk '{if($2==5) print "CB:Z:"$1}'>/data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/CD8Tex.txt
cat ~/eRNA/pbmcs/cells/reads/cluster.txt|awk '{if($2==6) print "CB:Z:"$1}'>/data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/ActNK.txt
cat ~/eRNA/pbmcs/cells/reads/cluster.txt|awk '{if ($2==9 ||$2==10) print "CB:Z:"$1}'>/data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/RestpDCs.txt
cat ~/eRNA/pbmcs/cells/reads/cluster.txt|awk '{if($2==11) print "CB:Z:"$1}'>/data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/MacrophagesM0.txt


###################cluster.sh
data=$1
nohup ~/anaconda2/envs/r-env/bin/samtools view -@ 48 /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/SC3_v3_NextGem_SI_PBMC_10K_possorted_genome_bam.bam | LC_ALL=C grep -F -f /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/$data.txt|grep 'UB:Z' > f_filtered_$data &
wait;
nohup cat SAM_header f_filtered_$data > f_filtered_$data.sam &
wait;
# Convert filtered.sam to BAM format
nohup ~/anaconda2/envs/r-env/bin/samtools view -u -@ 48 -b -S f_filtered_$data.sam -h -o f_filtered_$data.bam &
wait;
nohup ~/anaconda2/envs/r-env/bin/samtools index -@ 48 f_filtered_$data.bam &
wait;
#计算ccre reads
nohup bedtools multicov -bams f_filtered_$data.bam -bed GRCh38-ccREs.dELS_minus.bed -s >minus_dELS_reads_$data.txt &
nohup bedtools multicov -bams f_filtered_$data.bam -bed GRCh38-ccREs.dELS_minus.bed -S >plus_dELS_reads_$data.txt &
nohup bedtools multicov -bams f_filtered_$data.bam -bed GRCh38-ccREs.dELS_minus.bed >dELS_reads_$data.txt &
############
for q in {'CD8Tcells','Monocytes','RestMemCD4Tcells','NaiveBcells','CD8Tex','ActNK','RestpDCs','MacrophagesM0'}
do
	echo "nohup bash cluster.sh ${q} &"
done




#################
F5.hg38.enhancers.bed
bedtools intersect -a filtered_dELS.bed -b F5.hg38.enhancers.bed >pos_dELS.bed

bedtools intersect -a all.bed -b F5.hg38.enhancers.bed |wc -l



http://genome-asia.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr5%3A172766651%2D172766763&hgsid=746050145_nAMBwB6BDDKQRjXfYEeAEgJLaVAq



#nohup bedtools multicov -bams f_filtered.bam -bed ccREs.dELS_minus.bed >dELS_reads_over.txt &

nohup bedtools multicov -bams f_filtered.bam -bed pos_dELS.bed >dELS_reads_pos.txt &


nohup Rscript /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/all_cut.R pos &
#########cut dot
for data in 'pos'
do

	d=$(cat dELS_reads_dot_$data.txt|awk '{print $0}')
	awk '{if($7 >n)  print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$6}' n=$d dELS_reads_$data.txt >$data.bed

done

#####################just >0
for q in {'CD8Tcells','Monocytes','RestMemCD4Tcells','NaiveBcells','CD8Tex','ActNK','RestpDCs','MacrophagesM0'}
do
	nohup bedtools intersect -a plus_dELS_reads_$q.txt -b pos.bed |awk '{if($7>0) print $1"\t"$2"\t"$3"\t"$4"\t0\t-\t"$7}'>fpos_plus_dELS_$q.bed &
	nohup bedtools intersect -a minus_dELS_reads_$q.txt -b pos.bed |awk '{if($7>0) print $1"\t"$2"\t"$3"\t"$4"\t0\t-\t"$7}'>fpos_minus_dELS_$q.bed &
done

wait;
####################overlap plus_minus
for q in {'CD8Tcells','Monocytes','RestMemCD4Tcells','NaiveBcells','CD8Tex','ActNK','RestpDCs','MacrophagesM0'}
do
	nohup bedtools intersect -a fpos_plus_dELS_$q.bed -b fpos_minus_dELS_$q.bed|awk '{print $1"\t"$2"\t"$3"\t"$4"\t0\t.\t"$7}'>pos_$q.bed &
	nohup bedtools intersect -a fpos_minus_dELS_$q.bed -b fpos_plus_dELS_$q.bed|awk '{print $1"\t"$2"\t"$3"\t"$4"\t0\t.\t"$7}'>pos_minus_$q.bed &


done
wait;
####################plot for rna seq
nohup bamCoverage -bs 10 -b f_filtered.bam --filterRNAstrand forward -o filtered_forward.bigWig & 负
nohup bamCoverage -bs 10 -b f_filtered.bam -o filtered.bigWig &


nohup bamCoverage -bs 10 -b f_filtered.bam --filterRNAstrand reverse -o filtered_reverse.bigWig &正
for q in {'CD8Tcells','Monocytes','RestMemCD4Tcells','NaiveBcells','CD8Tex','ActNK','RestpDCs','MacrophagesM0'}
do
	nohup bamCoverage -bs 10 -b f_filtered_$q.bam --filterRNAstrand forward -o filtered_forward_$q.bigWig &
	nohup bamCoverage -bs 10 -b f_filtered_$q.bam --filterRNAstrand reverse -o filtered_reverse_$q.bigWig &
	nohup bamCoverage -bs 10 -b f_filtered_$q.bam -o filtered_$q.bigWig &


done







for q in {'CD8Tcells','Monocytes','RestMemCD4Tcells','NaiveBcells','CD8Tex','ActNK','RestpDCs','MacrophagesM0'}
do

	nohup bwtool agg 300:300 pos_$q.bed filtered_$q.bigWig fpos_$q.txt &
	nohup bwtool agg 300:300 pos_$q.bed filtered_reverse_$q.bigWig fpos_plus_$q.txt &
	nohup bwtool agg 300:300 pos_$q.bed filtered_forward_$q.bigWig  fpos_minus_$q.txt &
done


wait;



	nohup bwtool agg 300:300 pos.bed filtered.bigWig fpos.txt &
	nohup bwtool agg 300:300 pos.bed filtered_reverse.bigWig fpos_plus.txt &
	nohup bwtool agg 300:300 pos.bed filtered_forward.bigWig  fpos_minus.txt &






for q in {'CD8Tcells','Monocytes','RestMemCD4Tcells','NaiveBcells','CD8Tex','ActNK','RestpDCs','MacrophagesM0'}
do
	nohup Rscript fpos.R $q &

done


#######################plot for chip
for q in {'CD8Tcells','Monocytes','RestMemCD4Tcells','NaiveBcells','CD8Tex','ActNK','RestpDCs','MacrophagesM0'}
do
	ls -R /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/chip/$q/*.bigWig > /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/chip/$q/$q.txt
done



for q in {'CD8Tcells','Monocytes','RestMemCD4Tcells','NaiveBcells','CD8Tex','ActNK','RestpDCs','MacrophagesM0'}
do
    n=$(cat /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/chip/$q/$q.txt|wc -l)

    for k in $(seq $n)
    do
        peakFile=$(cat /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/chip/$q/$q.txt | awk -F "\t" '{if (NR == '$k') print $1}')
        nohup bwtool agg 500:500 /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/pos_$q.bed $peakFile $peakFile.pos.txt &
        nohup bwtool agg 500:500 /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/GRCh38-ccREs.dELS.bed $peakFile $peakFile.txt &

    done
done
wait;
cd /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/chip


for q in {'CD8Tcells','Monocytes','RestMemCD4Tcells','NaiveBcells','CD8Tex','ActNK','RestpDCs','MacrophagesM0'}
do
    n=$(cat /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/chip/$q/$q.txt|wc -l)

    for k in $(seq $n)
    do
        peakFile=$(cat /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/chip/$q/$q.txt | awk -F "\t" '{if (NR == '$k') print $1}')
        nohup Rscript chip_t.R $peakFile pos &

    done
done



for q in {'CD8Tcells','Monocytes','RestMemCD4Tcells','NaiveBcells','CD8Tex','ActNK','RestpDCs','MacrophagesM0'}
do

    wc -l /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/pos_$q.bed 
   
done



