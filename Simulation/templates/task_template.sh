#!/bin/bash

export SCRATCHDIR=PROD_BASEDIR/scratch/taskid/
mkdir -p $SCRATCHDIR

source G4_INSTALLDIR/bin/geant4.sh
source ROOT_INSTALLDIR/bin/thisroot.sh

source WCSIM_INSTALLDIR/setup.sh

mkdir $SCRATCHDIR/mPMT-configfiles 
cp WCSIM_INSTALLDIR/mPMT-configfiles/mPMTconfig_Position_WCTE_RotateBarrelHalfTower.txt $SCRATCHDIR/mPMT-configfiles
cp WCSIM_INSTALLDIR/macros/daq.mac                 $SCRATCHDIR
cp WCSIM_INSTALLDIR/macros/jobOptions.mac          $SCRATCHDIR
cp WCSIM_INSTALLDIR/macros/tuning_parameters.mac   $SCRATCHDIR
cp WCSIM_INSTALLDIR/lib/libWCSimRoot_rdict.pcm     $SCRATCHDIR

cp FQTUNER_INSTALLDIR/bin/WCSim_FQTuner $SCRATCHDIR

cd $SCRATCHDIR
./WCSim_FQTuner macrofile

rm -rf $SCRATCHDIR
