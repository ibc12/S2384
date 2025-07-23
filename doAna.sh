#!/bin/bash

# echo "Travaux en cours sur le code d'analyse, ne pas executer, work in progress on the analysis code, do not attempt to launch this script"
## Do all the steps of the analysis (in chain and only if the previous finished successfully)

## 1-> Read the TPC
actroot -r tpc &&

## 2-> Read the silicons
actroot -r sil &&
#
# ## 3-> Do the filter
actroot -f &&
#
# ## 4-> Do the merging
actroot -m
