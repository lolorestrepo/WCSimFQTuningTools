#!/bin/bash

export SCRATCHDIR=PROD_BASEDIR/scratch/taskid/
mkdir -p $SCRATCHDIR

source G4_INSTALLDIR/bin/geant4.sh
source ROOT_INSTALLDIR/bin/thisroot.sh

source WCSIM_INSTALLDIR/setup.sh

# export LD_LIBRARY_PATH=$(geant4-config --prefix)/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$THRONG_DIR/Software/HepMC3/install/lib64/

cd $SCRATCHDIR
mkdir -p data
cp WCSIM_INSTALLDIR/data/mPMT_Position_WCTE.txt    $SCRATCHDIR/data
cp WCSIM_INSTALLDIR/data/CathodeParameters.txt     $SCRATCHDIR/data
cp WCSIM_INSTALLDIR/macros/daq.mac                 $SCRATCHDIR
cp WCSIM_INSTALLDIR/macros/jobOptions.mac          $SCRATCHDIR
cp WCSIM_INSTALLDIR/macros/tuning_parameters.mac   $SCRATCHDIR
cp WCSIM_INSTALLDIR/lib/libWCSimRoot_rdict.pcm     $SCRATCHDIR

cp FQTUNER_INSTALLDIR/bin/WCSim_FQTuner $SCRATCHDIR

./WCSim_FQTuner macrofile

rm -rf $SCRATCHDIR