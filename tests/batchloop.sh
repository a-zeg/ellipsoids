#!/bin/bash

for aa in 0 2 # 100 200 300 400 500 600 700 800 900
  do
    scommand="sbatch jobsTurkevs.sh $aa 2"
    echo "submit command: $scommand"
    $scommand
done
