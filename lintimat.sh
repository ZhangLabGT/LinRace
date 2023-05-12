ncells=256

echo "running experiments for ${ncells} cells."
profile_dir="sim/profile_${ncells}_mu_0.1_pd_0_run_1.txt"
tg_dir="sim/top_genes_${ncells}_mu_0.1_pd_0_run_1.txt"

java -jar ./LinTIMaT/bin/LinTIMaT.jar \
  -i "$profile_dir" \
  -gf "$tg_dir" \
  -gc 100 \
  -ob "results/lintimat_${ncells}_bin_tree.newick" \
  -on "results/lintimat_${ncells}_nonbin_tree.txt" \
  -mi 200000 \
  -ci 200000