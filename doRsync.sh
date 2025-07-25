#!/bin/bash

# Identity file 
SSH_KEY="$1"

# From where
remote="rcmp@142.90.96.93"

# Set directories
dirs=("Raw" "Cluster" "Data" "Filter" "Merger")

#Act on dirs
for d in "${dirs[@]}"; do
  echo "Rsyncing ${d}"
  rsync -e "ssh -i ${SSH_KEY}" -avz "${remote}:/home/rcmp/S2384/RootFiles/${d}/" "./RootFiles/{d}/"
done
