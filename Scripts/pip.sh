#####example for mesc

#########step 1: use seruat to filter data
################dir/barcodes.tsv dir/genes.tsv dir/matrix.mtx dir/data/bam
source activate r4-base
#########get fantom enhancer
wget https://fantom.gsc.riken.jp/5/datafiles/reprocessed/mm10_latest/extra/enhancer/F5.mm10.enhancers.bed.gz
gzip -d *
cd ~/eRNA/refdata-gex-mm10-2020-A
#############get ccre_dels overlapped with fantom enhancer
cat /data/projects/encode/Registry/V2/mm10/mm10-ccREs.dELS.bed|awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t-"}'>mm10-ccREs.dELS_minus.bed
bedtools intersect -a mm10-ccREs.dELS_minus.bed -b F5.mm10.enhancers.bed >pos_dELS.bed
##############################run split bam by barcodes and get signal for positive enhancer 
nohup bash ~/eRNA/script/step_bam.sh ~/eRNA/mesc/Mesendoderm ~/eRNA/refdata-gex-mm10-2020-A/pos_dELS.bed ~/eRNA/mesc/Mesendoderm/data/mesendoderm_txm.possorted_genome_bam.bam  &
nohup bash ~/eRNA/script/step_bam.sh ~/eRNA/mesc/Germ_layer ~/eRNA/refdata-gex-mm10-2020-A/pos_dELS.bed ~/eRNA/mesc/Germ_layer/data/germlayer_txm.possorted_genome_bam.bam &
nohup bash ~/eRNA/script/step_bam.sh ~/eRNA/mesc/mix ~/eRNA/refdata-gex-mm10-2020-A/pos_dELS.bed~/eRNA/mesc/mix/data/mixed_txm.possorted_genome_bam.bam &

