1)

	# important export: needed to load fitQun parameters inside makehistWCSim.cc
	export FITQUN_ROOT=/pbs/home/g/gdiazlop/Software/fiTQun/fiTQun 

	make makehistWCSim

	# copy charge/angular pdf file ?

	makehistWCSim simulation_filename, wcte_parameteroverridefile

2)