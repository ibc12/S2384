#!/bin/bash

## Set data.conf file
file="./configs/data.conf"

## Check for empty runs
if [[ $# -eq 0 ]]; then 
  echo "Specify run or list of runs (separated by spaces and in increasing order)"
  exit 1
fi

## Parse commands
runs=""
for run in "$@"; do
	echo "Setting run : $run"
  runs+="$run, "
done

## Modify file with sed
sed -i "s/Runs:.*/Runs: $runs/" $file
