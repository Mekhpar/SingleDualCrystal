#!/bin/bash

cd /data/users/mparanjp/DD4hep_test/DD4hep/examples/SingleDualCrystal/compact/

START_TIME=`/bin/date`
echo "started at $START_TIME"


# setup software environment at UMD
#

source /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc11-opt/setup.sh
echo "ran setup"
source  /data/users/mparanjp/DD4hep_test/DD4hep/bin/thisdd4hep.sh
echo "ran thisdd4hep"
#
# run
#


ddsim --steeringFile /data/users/mparanjp/DD4hep_test/DD4hep/examples/SingleDualCrystal/compact/junkst2.py --compact /data/users/mparanjp/DD4hep_test/DD4hep/examples/SingleDualCrystal/compact/DR_e-.xml --runType batch --part.userParticleHandler='' -G --gun.position="0.,-30.,-64.33" --gun.direction "0 1 2.144" --gun.energy "10*GeV" --gun.particle="e-" -N 100 -O Root_files/PbWO4_e-/10_GeV/out_65_0309.root >& out_cndr_65.log

exitcode=$?


#
# end run
#

echo ""
END_TIME=`/bin/date`
echo "finished at $END_TIME"
exit $exitcode
