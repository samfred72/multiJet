#!/bin/bash
# Appends MC rows (isMC=1) to the shared queue.list used by condor.job.
# `rm -f queue.list` first if you want a fresh list rather than adding to
# whatever makequeuelist.sh already put there.
touch queue.list
while IFS= read -r line || [[ -n "$line" ]]; do
  if [[ $line == "#"* ]]; then
    continue
  fi
  read -r SIM TYPE remaining_text <<< "$line"
  for j in `cat /sphenix/user/samfred/projects/queuelists/masterqueue_${SIM}28_${TYPE}.list`;do
    echo 1 $SIM $j $TYPE >> queue.list
  done
done < "MCrunlist.list"
