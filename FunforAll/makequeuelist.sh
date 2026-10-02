#!/bin/bash
# Builds the combined queue.list used by condor.job, from one merged input
# file (runlist.list): plain run numbers for data, "<sim> <trigger>" pairs
# for MC (one token per line = data, two tokens = MC; '#' lines and blank
# lines are comments). Replaces the old separate makequeuelist.sh /
# makeMCqueuelist.sh.
SPLIT=$1
DATASET=$2
DST=$3

touch queue.list

while IFS= read -r line || [[ -n "$line" ]]; do
  [[ -z "$line" || "$line" == "#"* ]] && continue

  read -r tok1 tok2 rest <<< "$line"

  if [[ -z "$tok2" ]]; then
    # single token -> data run number
    cat /sphenix/user/samfred/projects/queuelists/masterqueue_${DATASET}_${DST}_split${SPLIT}.list | grep "$tok1" | awk '{print "0 none", $0, "none"}' >> queue.list
  else
    # two tokens -> MC sim + trigger
    SIM="$tok1"
    TYPE="$tok2"
    for j in `cat /sphenix/user/samfred/projects/queuelists/masterqueue_${SIM}28_${TYPE}.list`; do
      echo 1 $SIM $j $TYPE >> queue.list
    done
  fi
done < "runlist.list"
