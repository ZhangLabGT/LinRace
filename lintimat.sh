ncells=$1
run=$2
mu=$3
pd=$4

echo "running experiments for ${ncells} cells. run ${run}"
profile_dir="sim/profile_${ncells}_mu_0.1_pd_1_run_${run}.txt"
tg_dir="sim/top_genes_${ncells}_mu_0.1_pd_1_run_${run}.txt"

start_time="$(date -u +%s)"

Rscript ./lintimat_gen.R "$ncells" "$run" "$mu" "$pd"

java -jar ./LinTIMaT/bin/LinTIMaT.jar \
  -i "$profile_dir" \
  -gf "$tg_dir" \
  -gc 100 \
  -ob "lintimat_results/lintimat_${ncells}_${run}_bin_tree.newick" \
  -on "lintimat_results/lintimat_${ncells}_${run}_nonbin_tree.txt" \
  -mi 20000 \
  -ci 20000

end_time="$(date -u +%s)"
elapsed="$(($end_time-$start_time))"

Rscript ./rf_lintimat.R "${ncells}" "${run}" "${elapsed}"