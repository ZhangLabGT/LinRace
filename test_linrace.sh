for mu in 0.05 0.1 0.15; do
  for ncells in 256 512 1024 2048; do
    for run in 1 2 3 4 5 6 7 8 9 10; do
      for pd in 0 1; do
        echo "starting mu=$mu ncells=$ncells pd=$pd run=$run"
        Rscript ./Linrace_compare.R "$ncells" "$mu" "$pd" "$run"
        echo "finished."
      done
    done
  done
done
