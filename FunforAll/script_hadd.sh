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
RUNNUMBER=$3
TRIGGER=$4

if [[ $ISMC == 1 ]]; then
  hadd -f /sphenix/tg/tg01/jets/samfred/multiJet_full_hadded/multijet_${SIM}_${TRIGGER}.root /sphenix/tg/tg01/jets/samfred/multiJet_MC/*${SIM}_${TRIGGER}.root
else
  hadd -f /sphenix/tg/tg01/jets/samfred/multiJet_hadded/multijet_${RUNNUMBER}.root /sphenix/tg/tg01/jets/samfred/multiJet/*${RUNNUMBER}*
fi

echo all done
echo "script done"
