#!/bin/bash

export FQTUNINGTOOLS=$HOME/Software/WCSimFQTuningTools

#. $HOME/HK_Software/fiTQun/install-Darwin_arm64-gcc_15.0.0-python_3.10.13/setup.sh
. $HOME/Software/fiTQun/install-Linux_x86_64-gcc_9-python_3.10.13/setup.sh

#runfiTQun -p $FQTUNINGTOOLS/Time/fiTQun/WCTE_Parameters.dat -r out_fitqun.root $FQTUNINGTOOLS/Simulation/local/out/out_300.root

runfiTQun -p nuPRISMBeamTest_16cShort_mPMT.parameters.dat -r out_fitqun.root $FQTUNINGTOOLS/Simulation/local/out/out_300.root
