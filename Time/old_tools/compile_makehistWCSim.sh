# important export: needed to load fitQun parameters inside makehistWCSim.cc
source $HOME/Software/ROOT/install/bin/thisroot.sh
export FITQUN_ROOT=/pbs/home/g/gdiazlop/Software/HK_Software/fiTQun/install-Linux_x86_64-gcc_9-python_3.10.13/
export WCSIMDIR=/pbs/home/g/gdiazlop/Software/HK_Software/WCSim/install-Linux_x86_64-gcc_9-python_3.10.13/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$FITQUN_ROOT/lib:$WCSIMDIR/lib

make makehistWCSim

# redefine FITQUN_ROOT to find fiTQun.parameters.dat
export FITQUN_ROOT=/pbs/home/g/gdiazlop/Software/HK_Software/fiTQun/

# makehistWCSim simulation_filename, wcte_parameteroverridefile
./makehistWCSim $LUSTRE/Time/out/out_e-_1000_1.root ./WCTE_Parameters.dat
rm makehistWCSim.o makehistWCSim