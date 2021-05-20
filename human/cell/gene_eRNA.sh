dir=$1
gtf=$2
bed=$3





cat ${dir}/top_0.05.txt|grep EH38D|sort -k 8 -u>${dir}/top_0.05_erna.txt
cat ${gtf}|awk '{if($3 =="gene") print $1"\t"$4"\t"$5"\t"$16}'|grep 'chr'>${dir}/gene.bed
cat ${bed}|sort -k 4 >${dir}/ccre.bed
cat ${dir}/top_0.05_erna.txt|awk '{print $8}'|awk -F "-" '{print $2}'|sort >${dir}/select.bed
join -1 1 -2 4 ${dir}/select.bed ${dir}/ccre.bed|awk '{print $2"\t"$3-100000"\t"$4+100000"\t"$1"\t"$6}'> ${dir}/select_ccre.bed
 bedtools intersect -a ${dir}/select_ccre.bed -b ${dir}/gene.bed -wo >${dir}/select_gene.bed
 cat ${dir}/select_gene.bed|awk '{print $1"-"$4"\t"$9}'| sed 's/\"//g'| sed 's/\;//g'|sort -k 2 >${dir}/gene_ccre.txt
sort -k 8 ${dir}/top_0.05.txt>${dir}/top_0.05_s.txt
join -1 2 -2 8 -o 1.1 1.2 ${dir}/gene_ccre.txt ${dir}/top_0.05_s.txt|sort -u >${dir}/s_gene_id.txt

