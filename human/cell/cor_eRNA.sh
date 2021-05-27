source activate r4-base
dir=$1

cat ${dir}/top_0.05_cor_cut.txt|awk -F "EH38" '{print "EH38"$2}'|sort >${dir}/select_cut.bed
join -1 1 -2 4 -o 2.1 2.2 2.3 2.4 ${dir}/select_cut.bed ${dir}/ccre.bed> ${dir}/select_ccre_cut.bed
join -1 1 -2 4 -o 1.1 1.2 2.1 2.2 2.3 ${dir}/select_cut.bed ${dir}/select_ccre_cut.bed|sort -u> ${dir}/final.bed
nohup Rscript ~/eRNA/script/cor_ne_eRNA.R ${dir} &
cat ${dir}/top_0.05_cor.txt|awk '{if($3>0.9)print $0}'>${dir}/top_0.05_cor_0.9.txt


cat $${dir}/top_0.05_cor_0.9.txt|awk -F "EH38" '{print "EH38"$2}'|sort >${dir}/select_cut_0.9.bed
join -1 1 -2 4 -o 2.1 2.2 2.3 2.4 ${dir}/select_cut_0.9.bed ${dir}/ccre.bed> ${dir}/select_ccre_cut_0.9.bed
join -1 1 -2 4 -o 1.1 1.2 2.1 2.2 2.3 ${dir}/select_cut_0.9.bed ${dir}/select_ccre_cut_0.9.bed|sort -u> ${dir}/final_0.9.bed
