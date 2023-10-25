#!/bin/bash

# Compile create_2D_histogram.cc 
export FITQUNDIR=/pbs/home/g/gdiazlop/Software/HK_Software/fiTQun/install-Linux_x86_64-gcc_9-python_3.10.13/
export WCSIMDIR=/pbs/home/g/gdiazlop/Software/HK_Software/WCSim/install-Linux_x86_64-gcc_9-python_3.10.13/
export ROOTDIR=/pbs/home/g/gdiazlop/Software/ROOT/install/

g++ -g -O -fpic -I. -I$ROOTDIR/include -I$FITQUNDIR/../ -I$WCSIMDIR/include/WCSim -pthread -std=c++14 -m64 -Dlinux -D__linux__ -DNOSKLIBRARIES -DHEMI_CUDA_DISABLE -o create_2D_histogram create_2D_histogram.cc -L$ROOTDIR/lib -L$FITQUNDIR/lib -L$WCSIMDIR/lib -Wl,-z -Wl,muldefs -Wl,-rpath,$ROOTDIR/lib -Wl,-rpath,$FITQUNDIR/lib -Wl,-rpath,$WCSIMDIR/lib -lm -lstdc++ -lGui -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -pthread -ldl -rdynamic -lTreePlayer -lMinuit -lnsl -lfiTQunLib -lWCSimRoot -lstdc++fs -g

# FITQUN_ROOT is needed by fiTQun, which is called in the script
export FITQUN_ROOT=/pbs/home/g/gdiazlop/Software/HK_Software/fiTQun/
# Run the script
./create_2D_histogram $LUSTRE/Time/e-/out e- WCTE_Parameters.dat 1.38 1 5000 100 -10 10 125 -2. 3.
rm create_2D_histogram
