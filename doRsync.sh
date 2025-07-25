#!/bin/bash

# From where
remote="rcmp@142.90.96.93"

# Set directories
dirs=("Raw" "Cluster" "Data" "Filter" "Merger")

#Act on dirs
for d in "${dirs[@]}"; do
  echo "Rsyncing ${d}"
  rsync -e "ssh -i ~/.ssh/rcmp_at_triumf" -avz "${remote}:/home/rcmp/S2384/RootFiles/${d}/" "./RootFiles/{d}/"
done
