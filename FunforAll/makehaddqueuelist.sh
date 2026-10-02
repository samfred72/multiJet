#!/bin/bash
# Builds hadd_queue.list for condor_hadd.job from the same runlist.list
# makequeuelist.sh reads (plain run numbers = data, "<sim> <trigger>" pairs
# = MC). Always writes 4 columns (isMC sim runnumber trigger) so condor_hadd
# job's Queue line has a uniform shape regardless of row type.
touch hadd_queue.list
while IFS= read -r line || [[ -n "$line" ]]; do
  [[ -z "$line" || "$line" == "#"* ]] && continue

  read -r tok1 tok2 rest <<< "$line"

  if [[ -z "$tok2" ]]; then
    # single token -> data run number
    echo "0 none $tok1 none" >> hadd_queue.list
  else
    # two tokens -> MC sim + trigger
    echo "1 $tok1 none $tok2" >> hadd_queue.list
  fi
done < "runlist.list"
