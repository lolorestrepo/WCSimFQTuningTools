1)
	# important export: needed to load fitQun parameters inside makehistWCSim.cc
	source $HOME/Software/ROOT/install/bin/thisroot.sh
	export FITQUN_ROOT=/pbs/home/g/gdiazlop/Software/HK_Software/fiTQun/install-Linux_x86_64-gcc_9-python_3.10.13/
	export WCSIMDIR=/pbs/home/g/gdiazlop/Software/HK_Software/WCSim/install-Linux_x86_64-gcc_9-python_3.10.13/
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$FITQUN_ROOT/lib:$WCSIMDIR/lib

	make makehistWCSim

	# redefine FITQUN_ROOT to find fiTQun.parameters.dat
	export FITQUN_ROOT=/pbs/home/g/gdiazlop/Software/HK_Software/fiTQun/

	# get needed fitqun files
	cd $FITQUN_ROOT
	source getTuneFiles.sh WCTE_tune_files.list
	
	# produced files by myself
	cp /produced/chargePDF/filename $FITQUN_ROOT/const/cPDFpar_3inchPMTR12199_02.root
	
	# run script
	# makehistWCSim simulation_filename, wcte_parameteroverridefile
	cd -
	./makehistWCSim $HOME/Software/HK_Software/WCSimFQTuningTools/Simulation/local/out/out.root ./WCTE_Parameters.dat

2) Run it for several energies and then run combhists.cc macro
3) Fits with fitpdf.cc macro