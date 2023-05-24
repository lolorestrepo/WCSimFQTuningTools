#!/bin/bash

export SCRATCHDIR=PROD_BASEDIR/scratch/jobid/
mkdir -p $SCRATCHDIR

source CONDA_INSTALLDIR/etc/profile.d/conda.sh
conda activate wcte
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CONDA_PREFIX/lib

source G4_INSTALLDIR/bin/geant4.sh
source ROOT_INSTALLDIR/bin/thisroot.sh

export WCSIMDIR=WCSIM_INSTALLDIR

mkdir $SCRATCHDIR/mPMT-configfiles 
cp $WCSIMDIR/mPMT-configfiles/mPMTconfig_Position_WCTE.txt $SCRATCHDIR/mPMT-configfiles
cp $WCSIMDIR/macros/daq.mac                 $SCRATCHDIR
cp $WCSIMDIR/macros/jobOptions.mac          $SCRATCHDIR
cp $WCSIMDIR/macros/tuning_parameters.mac   $SCRATCHDIR
cp $WCSIMDIR/lib/libWCSimRootDict_rdict.pcm $SCRATCHDIR

cp FQTUNER_INSTALLDIR/WCSim_FQTuner $SCRATCHDIR

cd $SCRATCHDIR
./WCSim_FQTuner macrofile

rm -rf $SCRATCHDIR
