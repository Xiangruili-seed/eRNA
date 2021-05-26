source activate r4-base
#########filter data
mkdir -p /home/lixiangr/eRNA/hesc/rep1/hesc/outs/filtered_feature_bc_matrix/data/
##nohup Rscript ~/eRNA/12878/step1_filter.R ~/eRNA/185021_GM12878/outs/filtered_feature_bc_matrix GM185021_GM12878 &
nohup Rscript ~/eRNA/12878/step1_filter.R ~/eRNA/hesc/rep1/hesc/outs/filtered_feature_bc_matrix HESC &
ln -s ~/eRNA/hesc/rep1/hesc/outs/possorted_genome_bam.bam ~/eRNA/hesc/rep1/hesc/outs/filtered_feature_bc_matrix/data/possorted_genome_bam.bam
ln -s ~/eRNA/hesc/rep1/hesc/outs/possorted_genome_bam.bam.bai ~/eRNA/hesc/rep1/hesc/outs/filtered_feature_bc_matrix/data/possorted_genome_bam.bam.bai


source activate r-env
nohup bash ~/eRNA/12878/step_bam.sh ~/eRNA/hesc/rep1/hesc/outs/filtered_feature_bc_matrix /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/pos_dELS.bed ~/eRNA/hesc/rep1/hesc/outs/filtered_feature_bc_matrix/data/possorted_genome_bam.bam &
wait;

#################get single cell reads for eRNA
source activate r-env
nohup bash ~/eRNA/12878/signal_eRNA.sh ~/eRNA/hesc/rep1/hesc/outs/filtered_feature_bc_matrix &

wait;
source activate r-env

nohup bash ~/eRNA/script/test_eRNA.sh ~/eRNA/hesc/rep1/hesc/ ~/eRNA/hesc/rep1/hesc/outs/filtered_feature_bc_matrix/data/f_dELS_reads.bed outs/filtered_feature_bc_matrix /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/pos_dELS.bed &
wait;
mkdir -p  ~/eRNA/hesc/rep1/hesc//outs/filtered_feature_bc_matrix//pos_add/

nohup bash ~/eRNA/script/matrix_eRNA.sh  ~/eRNA/hesc/rep1/hesc//outs/filtered_feature_bc_matrix/ &

wait;
########eRNA UMAP
source activate r4-base
nohup Rscript ~/eRNA/script/umap_eRNA.R ~/eRNA/hesc/rep1/hesc//outs/filtered_feature_bc_matrix//pos_add/ HESC &
wait;
#nohup Rscript ~/eRNA/script/umap_eRNA_1.R ~/eRNA/185021/outs/filtered_feature_bc_matrix/pos/ GM18502 &
#nohup Rscript ~/eRNA/script/umap_eRNA_1.R ~/eRNA/185021/outs/filtered_feature_bc_matrix/pos_add/ GM18502 &
###################eRNA-GEN PAIRS
nohup bash ~/eRNA/script/gene_eRNA.sh ~/eRNA/hesc/rep1/hesc//outs/filtered_feature_bc_matrix//pos_add/ ~/eRNA/refdata-gex-GRCh38-2020-A/genes/genes.gtf /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/GRCh38-ccREs.dELS.bed &

############################cut cor
wait;
nohup Rscript ~/eRNA/script/cor_eRNA.R ~/eRNA/hesc/rep1/hesc//outs/filtered_feature_bc_matrix//pos_add/ &

wait;
#########################ne
nohup bash ~/eRNA/script/cor_eRNA.sh ~/eRNA/hesc/rep1/hesc//outs/filtered_feature_bc_matrix//pos_add/ &
