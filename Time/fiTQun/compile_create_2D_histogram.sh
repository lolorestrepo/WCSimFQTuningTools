#!/bin/bash

# Compile create_2D_histogram.cc 
export FITQUNDIR=$HOME/Software/fiTQun/install-Linux_x86_64-gcc_9-python_3.10.13/
export WCSIMDIR=$HOME/Software/WCSim/install-Linux_x86_64-gcc_9-python_3.10.13/
export ROOTDIR=$THRONG_DIR/Software/ROOT/root_v6-28-00-patches/

g++ -g -O -fpic -I. -I$ROOTDIR/include -I$FITQUNDIR/../ -I$WCSIMDIR/include/WCSim -pthread -std=c++14 -m64 -Dlinux -D__linux__ -DNOSKLIBRARIES -DHEMI_CUDA_DISABLE -o create_2D_histogram create_2D_histogram.cc -L$ROOTDIR/lib -L$FITQUNDIR/lib -L$WCSIMDIR/lib -Wl,-z -Wl,muldefs -Wl,-rpath,$ROOTDIR/lib -Wl,-rpath,$FITQUNDIR/lib -Wl,-rpath,$WCSIMDIR/lib -lm -lstdc++ -lGui -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -pthread -ldl -rdynamic -lTreePlayer -lMinuit -lnsl -lfiTQunLib -lWCSimRoot -lstdc++fs -g

# FITQUN_ROOT is needed by fiTQun, which is called in the script
export FITQUN_ROOT=$HOME/Software/fiTQun/
# Run the script
./create_2D_histogram $LUSTRE/Time/e-/out e- WCTE_Parameters.dat 1.38 1 2500 100 -5 5 125 -2. 3.
rm create_2D_histogram
