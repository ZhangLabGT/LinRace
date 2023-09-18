for ncells in 256 512 1024 2048; do
  for run in 1 2 3 4 5 6 7 8 9 10; do
    for mu in 0.1; do
      for pd in 1; do
        sh lintimat.sh "${ncells}" "${run}" "${mu}" "${pd}" &
      done
    done
  done
done