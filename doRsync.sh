#!/bin/bash

# Identity file
if [ "$1" == "ibc" ]; then
  SSH_KEY="$HOME/.ssh/rcmp_at_triumf_ibc"
elif [ "$1" == "mlg" ]; then
  SSH_KEY="$HOME/.ssh/rcmp_at_triumf"
else
  echo "Error: first argument must be either ibc or mlg."
  exit 1
fi

# From where
remote="rcmp@142.90.96.93"

# Set directories
dirs=("Raw" "Cluster" "Data" "Filter" "Merger")

#Act on dirs
for d in "${dirs[@]}"; do
  echo "Rsyncing ${d}"
  rsync -e "ssh -i ${SSH_KEY}" -avz --progress "${remote}:/home/rcmp/S2384/RootFiles/${d}/" "./RootFiles/${d}/"
done
