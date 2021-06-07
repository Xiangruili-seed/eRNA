#####example for gm12878_gm185021
#####run cellranger
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-17/SRR8551678/SRR8551678.1
wait;
nohup fastq-dump --gzip --split-files SRR8551678.1 &

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-17/SRR8551676/SRR8551676.1
wait;
nohup fastq-dump --gzip --split-files SRR8551676.1 &


wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-17/SRR8551677/SRR8551677.1
wait;
nohup fastq-dump --gzip --split-files SRR8551677.1 &
i=SRR8551678
mv ${i}.1_1*.gz ${i}_S1_L001_R1_001.fastq.gz;
mv ${i}.1_2*.gz ${i}_S1_L001_R2_001.fastq.gz;

i=SRR8551677
mv ${i}.1_1*.gz ${i}_S1_L001_R1_001.fastq.gz;
mv ${i}.1_2*.gz ${i}_S1_L001_R2_001.fastq.gz;
i=SRR8551676
mv ${i}.1_1*.gz ${i}_S1_L001_R1_001.fastq.gz;
mv ${i}.1_2*.gz ${i}_S1_L001_R2_001.fastq.gz;

nohup ~/eRNA/cellranger-4.0.0/cellranger count --id=185021 --transcriptome=/home/lixiangr/eRNA/refdata-gex-GRCh38-2020-A --fastqs=/home/lixiangr/eRNA/GM185021/data --sample=SRR8551676 --nosecondary --jobmode=local &
nohup ~/eRNA/cellranger-4.0.0/cellranger count --id=12878 --transcriptome=/home/lixiangr/eRNA/refdata-gex-GRCh38-2020-A --nosecondary --fastqs=/home/lixiangr/eRNA/GM12878/data/ --sample=SRR8551677 &
nohup ~/eRNA/cellranger-4.0.0/cellranger count --id=185021_GM12878 --transcriptome=/home/lixiangr/eRNA/refdata-gex-GRCh38-2020-A --nosecondary --fastqs=/home/lixiangr/eRNA/GM185021_GM12878//data/ --sample=SRR8551678 &


#########step 1: use seruat to filter data
################dir/barcodes.tsv dir/genes.tsv dir/matrix.mtx dir/data/bam
source activate r4-base
#########filter data
##nohup Rscript ~/eRNA/12878/step1_filter.R ~/eRNA/185021_GM12878/outs/filtered_feature_bc_matrix GM185021_GM12878 &
nohup Rscript ~/eRNA/12878/step1_filter.R ~/eRNA/12878/outs/filtered_feature_bc_matrix GM12878 &
nohup Rscript ~/eRNA/12878/step1_filter.R ~/eRNA/185021/outs/filtered_feature_bc_matrix GM185021 &
wait;
##############################run split bam by barcodes and get signal for positive enhancer 
#nohup bash ~/eRNA/12878/step_bam.sh ~/eRNA/185021_GM12878/outs/filtered_feature_bc_matrix /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/pos_dELS.bed ~/eRNA/185021_GM12878/outs/filtered_feature_bc_matrix/data/possorted_genome_bam.bam  &
nohup bash ~/eRNA/12878/step_bam.sh ~/eRNA/12878/outs/filtered_feature_bc_matrix /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/pos_dELS.bed ~/eRNA/12878/outs/filtered_feature_bc_matrix/data/possorted_genome_bam.bam dELS &
nohup bash ~/eRNA/12878/step_bam.sh ~/eRNA/185021/outs/filtered_feature_bc_matrix /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/pos_dELS.bed ~/eRNA/185021/outs/filtered_feature_bc_matrix/data/possorted_genome_bam.bam dELS &
wait;

#################get single cell reads for eRNA
source activate r-env
nohup bash ~/eRNA/12878/signal_eRNA.sh ~/eRNA/185021_GM12878/outs/filtered_feature_bc_matrix dELS &
nohup bash ~/eRNA/12878/signal_eRNA.sh ~/eRNA/12878/outs/filtered_feature_bc_matrix dELS &
nohup bash ~/eRNA/12878/signal_eRNA.sh ~/eRNA/185021/outs/filtered_feature_bc_matrix dELS &
##############get chip signal for eRNA
nohup bash ~/eRNA/script/test_eRNA.sh ~/eRNA/12878/ ~/eRNA/12878/outs/filtered_feature_bc_matrix/data/f_dELS_reads.bed outs/filtered_feature_bc_matrix /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/pos_dELS.bed &
nohup bash ~/eRNA/script/test_eRNA.sh ~/eRNA/185021/ ~/eRNA/185021/outs/filtered_feature_bc_matrix/data/f_dELS_reads.bed outs/filtered_feature_bc_matrix /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/pos_dELS.bed &
########get eRNA matrix
nohup bash ~/eRNA/script/matrix_eRNA.sh ~/eRNA/12878/outs/filtered_feature_bc_matrix &
nohup bash ~/eRNA/script/matrix_eRNA.sh ~/eRNA/185021/outs/filtered_feature_bc_matrix &
nohup bash ~/eRNA/script/matrix_eRNA.sh ~/eRNA/185021_GM12878/outs/filtered_feature_bc_matrix &
########eRNA UMAP
source activate r4-base
nohup Rscript ~/eRNA/script/umap_eRNA.R ~/eRNA/12878/outs/filtered_feature_bc_matrix/pos_add/ GM12878 &
nohup Rscript ~/eRNA/script/umap_eRNA.R ~/eRNA/185021/outs/filtered_feature_bc_matrix/pos_add/ GM18502 &
nohup Rscript ~/eRNA/script/umap_eRNA.R ~/eRNA/185021_GM12878/outs/filtered_feature_bc_matrix/pos_add GM12878_GM18502 &
#nohup Rscript ~/eRNA/script/umap_eRNA_1.R ~/eRNA/185021/outs/filtered_feature_bc_matrix/pos/ GM18502 &
#nohup Rscript ~/eRNA/script/umap_eRNA_1.R ~/eRNA/185021/outs/filtered_feature_bc_matrix/pos_add/ GM18502 &
###################eRNA-GEN PAIRS
nohup bash ~/eRNA/script/gene_eRNA.sh ~/eRNA/12878/outs/filtered_feature_bc_matrix/pos_add/ ~/eRNA/refdata-gex-GRCh38-2020-A/genes/genes.gtf /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/GRCh38-ccREs.dELS.bed &
nohup bash ~/eRNA/script/gene_eRNA.sh ~/eRNA/185021/outs/filtered_feature_bc_matrix/pos_add/ ~/eRNA/refdata-gex-GRCh38-2020-A/genes/genes.gtf /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/GRCh38-ccREs.dELS.bed &
nohup bash ~/eRNA/script/gene_eRNA.sh ~/eRNA/185021_GM12878/outs/filtered_feature_bc_matrix/pos_add/ ~/eRNA/refdata-gex-GRCh38-2020-A/genes/genes.gtf /data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/GRCh38-ccREs.dELS.bed &

############################cut cor
nohup Rscript ~/eRNA/script/cor_eRNA.R ~/eRNA/12878/outs/filtered_feature_bc_matrix/pos_add/ &
nohup Rscript ~/eRNA/script/cor_eRNA.R ~/eRNA/185021/outs/filtered_feature_bc_matrix/pos_add/ &
nohup Rscript ~/eRNA/script/cor_eRNA.R ~/eRNA/185021_GM12878/outs/filtered_feature_bc_matrix/pos_add/ &

#########################ne
nohup bash ~/eRNA/script/cor_eRNA.sh ~/eRNA/12878/outs/filtered_feature_bc_matrix/pos_add/ &
nohup bash ~/eRNA/script/cor_eRNA.sh ~/eRNA/185021/outs/filtered_feature_bc_matrix/pos_add/ &
nohup bash ~/eRNA/script/cor_eRNA.sh ~/eRNA/185021_GM12878/outs/filtered_feature_bc_matrix/pos_add/ &




