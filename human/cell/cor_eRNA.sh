source activate r4-base
dir=$1

cat ${dir}/top_0.05_cor_cut.txt|awk -F "-" '{print $2}'|sort >${dir}/select_cut.bed
join -1 1 -2 4 -o 2.1 2.2 2.3 2.4 ${dir}/select_cut.bed ${dir}/ccre.bed> ${dir}/select_ccre_cut.bed
join -1 1 -2 4 -o 1.1 1.2 2.1 2.2 2.3 ${dir}/select_cut.bed ${dir}/select_ccre_cut.bed|sort -u> ${dir}/final.bed
nohup Rscript ~/eRNA/script/cor_ne_eRNA.R ${dir} &
