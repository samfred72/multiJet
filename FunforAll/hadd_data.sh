#!/bin/bash
# Merge the per-run hadded multiJet data trees into one data tree.
# Usage: ./hadd_data.sh [--dry-run]
# Run only after every per-run hadd job has finished: hadd merges whatever files exist.
INDIR=/sphenix/tg/tg01/jets/samfred/multiJet_hadded
PATTERN='multijet_[0-9][0-9][0-9][0-9][0-9].root'        # per-run data files only (no MC)
OUTDIR=/sphenix/tg/tg01/jets/samfred/multiJet_full_hadded
OUTFILE=$OUTDIR/multijet_Data.root
RUNLIST=/sphenix/user/samfred/projects/multiJet/FunforAll/runlist.list
MAXOPEN=500                  # hadd -n: inputs open at once (node limit: ulimit -n = 1024)

source /opt/sphenix/core/bin/sphenix_setup.sh -n new > /dev/null
set -euo pipefail   # after sourcing: the setup script reads unset variables

shopt -s nullglob
files=( $INDIR/$PATTERN )
shopt -u nullglob
nexpected=$(grep -cvE '^\s*$|^#|\S\s+\S' $RUNLIST || true)   # data runs = single-token lines

echo "input:    $INDIR/$PATTERN"
echo "files:    ${#files[@]} found, $nexpected data runs in $(basename $RUNLIST)"
echo "output:   $OUTFILE"

if (( ${#files[@]} == 0 )); then
  echo "ERROR: no input files" >&2; exit 1
fi
if (( ${#files[@]} != nexpected )); then
  echo "ERROR: file count does not match the run list - per-run hadds missing or extra; not merging." >&2
  exit 1
fi
if [[ ${1:-} == "--dry-run" ]]; then
  echo "dry run: nothing written"; exit 0
fi

mkdir -p $OUTDIR
tmp=$OUTFILE.tmp_$$
trap 'rm -f $tmp' EXIT
hadd -f -n $MAXOPEN $tmp "${files[@]}"
mv -f $tmp $OUTFILE
trap - EXIT
echo "merged ${#files[@]} files -> $OUTFILE"
