#!/bin/bash

# Identity file
if [ "$1" == "ibc" ]; then
  SSH_KEY="$HOME/.ssh/rcmp_at_triumf_ibc"
elif [ "$1" == "mlg" ]; then
  SSH_KEY="$HOME/.ssh/rcmp_at_triumf"
elif [ "$1" == "ld" ]; then
  SSH_KEY="$HOME/.ssh/rcmp_at_triumf_ld"
else
  echo "Error: first argument must be either ibc, mlg or ld"
  exit 1
fi

# From where
remote="rcmp@142.90.96.93"

# Skip first argument
shift

# And now read directiories if passed
if [ "$#" -eq 0 ]; then
  dirs=("Raw" "Cluster" "Data" "Filter" "Merger") # default case
else
  dirs=("$@")
fi

#Act on dirs
for d in "${dirs[@]}"; do
  echo "Rsyncing ${d}"
  rsync -e "ssh -i ${SSH_KEY}" -avz --progress "${remote}:/home/rcmp/S2384/RootFiles/${d}/" "./RootFiles/${d}/"
done
