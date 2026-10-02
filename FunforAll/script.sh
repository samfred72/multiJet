#!/bin/bash
export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${USER}

source /opt/sphenix/core/bin/sphenix_setup.sh -n new
source /opt/sphenix/core/bin/setup_local.sh /sphenix/user/samfred/projects/multiJet/install/

this_script=$BASH_SOURCE
this_script=`readlink -f $this_script`
this_dir=`dirname $this_script`
echo running: $this_script $*

ISMC=$1
SIM=$2
LIST=$3
TRIGGER=$4
PROCESS=$5
TEST=0
[ -n "$6" ] && TEST=$6

if [[ $ISMC == 1 ]]; then
  QUEUEDIR="/sphenix/user/samfred/projects/filelists/${SIM}28/${SIM}_${TRIGGER}"
else
  SAMPLE="goldenruns_pp24_calojet_split15"
  QUEUEDIR="/sphenix/user/samfred/projects/filelists/${SAMPLE}"
fi

# go to condor scratch directory
if [[ $TEST == 0 ]]; then

  if [[ ! -z "$_CONDOR_SCRATCH_DIR" && -d $_CONDOR_SCRATCH_DIR ]]; then
    cd $_CONDOR_SCRATCH_DIR
  else
    echo condor scratch NOT set
    exit -1
  fi

  # Works for both list formats: data's queue.list has one DST_TOWERS path
  # per line, MC's has 5 whitespace-separated DST paths per line - `cat`
  # word-splits either into one flat token stream either way.
  files=`cat ${QUEUEDIR}/${LIST}`
  echo $files
  for file in $files; do
    echo "Copying $file"
    getinputfiles.pl $file
  done

  cp /sphenix/user/samfred/projects/multiJet/FunforAll/Fun4All_macro.C .
fi

echo "input files..."
cat ${QUEUEDIR}/${LIST}

if [[ $ISMC == 1 ]]; then
  root -l -q -b "Fun4All_macro.C(\"${QUEUEDIR}/${LIST}\",${TEST},true,\"${TRIGGER}\",\"${SIM}\")"
else
  root -l -q -b "Fun4All_macro.C(\"${QUEUEDIR}/${LIST}\",${TEST},false)"
fi

if [[ $TEST == 0 ]]; then
  BASENAME="${LIST%.*}"
  if [[ $ISMC == 1 ]]; then
    FILENAME=${BASENAME}_${SIM}_${TRIGGER}.root
    cp outtree_$FILENAME /sphenix/tg/tg01/jets/samfred/multiJet_MC/
  else
    FILENAME=$BASENAME.root
    cp outtree_$FILENAME /sphenix/tg/tg01/jets/samfred/multiJet/
  fi
  rm /sphenix/user/samfred/projects/multiJet/FunforAll/log/${PROCESS}.out
  rm /sphenix/user/samfred/projects/multiJet/FunforAll/log/${PROCESS}.err
fi

echo all done
echo "script done"
