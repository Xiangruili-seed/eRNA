#####example for gm12878_gm185021

#########step 1: use seruat to filter data
################dir/barcodes.tsv dir/genes.tsv dir/matrix.mtx dir/data/bam
source activate r4-base
#########get fantom enhancer
nohup Rscript ~/eRNA/12878/step1_filter.R ~/eRNA/185021_GM12878/outs/filtered_feature_bc_matrix GM185021_GM12878 &
nohup Rscript ~/eRNA/12878/step1_filter.R ~/eRNA/12878/outs/filtered_feature_bc_matrix GM12878 &
nohup Rscript ~/eRNA/12878/step1_filter.R ~/eRNA/185021/outs/filtered_feature_bc_matrix GM185021 &


wait;
##############################run split bam by barcodes and get signal for positive enhancer 
#nohup bash ~/eRNA/12878/step_bam.sh ~/eRNA/185021_GM12878/outs/filtered_feature_bc_matrix /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/pos_dELS.bed ~/eRNA/185021_GM12878/outs/filtered_feature_bc_matrix/data/possorted_genome_bam.bam  &
nohup bash ~/eRNA/12878/step_bam.sh ~/eRNA/12878/outs/filtered_feature_bc_matrix /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/pos_dELS.bed ~/eRNA/12878/outs/filtered_feature_bc_matrix/data/possorted_genome_bam.bam &
nohup bash ~/eRNA/12878/step_bam.sh ~/eRNA/185021/outs/filtered_feature_bc_matrix /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/pos_dELS.bed ~/eRNA/185021/outs/filtered_feature_bc_matrix/data/possorted_genome_bam.bam &
wait;


source activate r-env

#nohup bash ~/eRNA/12878/signal_eRNA.sh ~/eRNA/185021_GM12878/outs/filtered_feature_bc_matrix &
nohup bash ~/eRNA/12878/signal_eRNA.sh ~/eRNA/12878/outs/filtered_feature_bc_matrix &
nohup bash ~/eRNA/12878/signal_eRNA.sh ~/eRNA/185021/outs/filtered_feature_bc_matrix &

nohup bash ~/eRNA/script/test_eRNA.sh ~/eRNA/12878/ ~/eRNA/12878/outs/filtered_feature_bc_matrix/data/f_dELS_reads.bed outs/filtered_feature_bc_matrix /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/pos_dELS.bed &
nohup bash ~/eRNA/script/test_eRNA.sh ~/eRNA/185021/ ~/eRNA/185021/outs/filtered_feature_bc_matrix/data/f_dELS_reads.bed outs/filtered_feature_bc_matrix /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/pos_dELS.bed &
