#!/bin/bash

export SCRATCHDIR=PROD_BASEDIR/scratch/taskid/
mkdir -p $SCRATCHDIR

source G4_INSTALLDIR/bin/geant4.sh
source ROOT_INSTALLDIR/bin/thisroot.sh

export WCSIMDIR=WCSIM_INSTALLDIR
source $WCSIMDIR/setup.sh

mkdir $SCRATCHDIR/mPMT-configfiles 
cp $WCSIMDIR/mPMT-configfiles/mPMTconfig_Position_WCTE.txt $SCRATCHDIR/mPMT-configfiles
cp $WCSIMDIR/macros/daq.mac                 $SCRATCHDIR
cp $WCSIMDIR/macros/jobOptions.mac          $SCRATCHDIR
cp $WCSIMDIR/macros/tuning_parameters.mac   $SCRATCHDIR
cp $WCSIMDIR/lib/libWCSimRootDict_rdict.pcm $SCRATCHDIR

cp FQTUNER_INSTALLDIR/bin/WCSim_FQTuner $SCRATCHDIR

cd $SCRATCHDIR
./WCSim_FQTuner macrofile

rm -rf $SCRATCHDIR
