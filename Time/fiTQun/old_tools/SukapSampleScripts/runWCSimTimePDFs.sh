#!/bin/bash

SETUP_SCRIPT=/home/tyoshida/nuPRISM/Analysis/Source_At_Start_nuPRISM.sh

#WCSIM_BASE_DIR=/home/tyoshida/nuPRISM/WCSim
WCSIM_BASE_DIR=/home/tyoshida/fiTQunTuning/WCSimFQTuner
#WCSIM=${WCSIM_BASE_DIR}/exe/bin/Linux-g++/WCSim
WCSIM=${WCSIM_BASE_DIR}/WCSim_FQTuner

DETECTOR_NAME=nuPRISM
#DETECTOR_NAME=nuPRISM_mPMT
PMT_NAME=PMT8inch
SUFFIX=NuPRISM_10x8_8inchPMT_40perCent


if [ $# -eq 0 ]
then
    echo "No arguments supplied, detector name is "${DETECTOR_NAME}
    echo "Need cylinder height and radius to generate particle guns!"
    exit
elif [ $# -eq 3 ]
then
    echo "Detector name given as an argument: "$1
    DETECTOR_NAME=$1
    echo "Particle gun vertex cylinder half-height given as an argument: "$2
    DETECTOR_H=$2
    echo "Particle gun vertex cylinder radius given as an argument: "$3
    DETECTOR_R=$3
    # 1510 1490 // SK
    # 1810 1684 // SK full ID volume
    # 350 150   // nuPRISM
    # 2800 3500 // 60x74 Cylinder
fi

OUT_DIR=/disk01/usr5/tyoshida/fiTQunTuning/${SUFFIX}_tbugfix/timepdf
#BASE_WCSIMMACRO=/home/tyoshida/fiTQunTuning/Utilities/timepdf/SukapSampleScripts/timePDF_base_nuPRISM_mPMT.mac
BASE_WCSIMMACRO=/home/tyoshida/fiTQunTuning/Utilities/timepdf/SukapSampleScripts/timePDF_base_nuPRISM.mac
#BASE_WCSIMMACRO=/home/cvilela/fiTQunTuning/TimePDF/Utilities/timepdf/SukapSampleScripts/timePDF_base.mac
CHART_FILE_DIR=/home/tyoshida/fiTQunTuning/Utilities/timepdf


source $SETUP_SCRIPT

echo "Script starts"

# Queue and output dirs
#QUEUE_NAME=all
QUEUE_NAME=atmpd


#STDOUT_DIR=/storage/shared/cvilela/qsub_output
#STDERR_DIR=/storage/shared/cvilela/qsub_error

# List of particles
#P_NAMES=(e- mu- proton pi+ kaon+)
#PDG_CODES=(11 13 2212 211 321)
#P_MASSES=(0.511 105.7 938.3 139.6 493.7)
#PROD_NUM=0

#P_NAMES=(e- mu- pi+)
#PDG_CODES=(11 13 211)
#P_MASSES=(0.511 105.7 139.6)
#PROD_NUM=0

#P_NAMES=(e-)
#PDG_CODES=(11)
#P_MASSES=(0.511)
#PROD_NUM=0

P_NAMES=(pi+)
PDG_CODES=(211)
P_MASSES=(139.6)
PROD_NUM=0

#P_NAMES=(mu- pi+)
#PDG_CODES=(13 211)
#P_MASSES=(105.7 139.6)
#PROD_NUM=0
#P_NAMES=(mu-)
#PDG_CODES=(13)
#P_MASSES=(105.7)
#PROD_NUM=0
#P_NAMES=(e- pi+)
#PDG_CODES=(11 211)
#P_MASSES=(0.511 139.6)
#PROD_NUM=0
#P_NAMES=(proton kaon+)
#PDG_CODES=(2212 321)
#P_MASSES=(938.3 493.7)
#PROD_NUM=0

echo "Number of particle types " ${#P_NAMES[*]} 

# Particle loop
for ((i=0; i<${#P_NAMES[*]}; i++))
do
    echo 'Starting ' ${P_NAMES[$i]}

    PDG_CODE=${PDG_CODES[$i]}
    PARTICLE_NAME=${P_NAMES[$i]}
    PARTICLE_MASS=${P_MASSES[$i]}

    CHART_FILE=${CHART_FILE_DIR}/chart_${PDG_CODE}.txt
    
    ls ${CHART_FILE}
    
    COPY_NUMBER=0
    MOM=-1
    
    # Momentum Loop
    while read line 
    do
	echo "Processing line " $line
	read -a lineArr <<< $line
	
	# Number of runs in a specific job
	N_RUNS="${lineArr[0]}"
	N_JOBS="${lineArr[1]}"
	N_EVENTS="${lineArr[2]}"
	
	echo "Runs" $N_RUNS "Jobs" $N_JOBS "Events" $N_EVENTS
	
	# Job loop
	for job in `seq 0 $(($N_JOBS-1))` 
	do
	    # Run Loop
	    #	echo " Processing job" $job 
	    for run in `seq 0 $(($N_RUNS-1))`
	    do
		if [ $job -eq 0  ] 
		then
		    echo '++++++++++++++++++= is zeroth job'
		    if [  $MOM -eq ${lineArr[$((3+$run))]} ]
		    then
			COPY_NUMBER=$(( $COPY_NUMBER + 1 ))
			echo "++++++++++++++++++++++ COPY NUMBER INCREMENTED " ${COPY_NUMBER}
		    else
			COPY_NUMBER=0
		    fi
		fi
		#	    echo "  Processing run" $run
		MOM="${lineArr[$((3+$run))]}"
		SEED1=$(( (${job}*10+$COPY_NUMBER+1)*$MOM))
		SEED2=$((($job+1)*$PDG_CODE))
		
		CURNAME=${PDG_CODE}_${MOM}_${COPY_NUMBER}_${job}_${PROD_NUM}
		
		echo "Making dir " ${OUT_DIR}/${CURNAME}
		mkdir -p ${OUT_DIR}/${CURNAME}
		
		THIS_MACRO=${OUT_DIR}/${CURNAME}/tPDF_${CURNAME}.mac
		ROOT_FNAME=${OUT_DIR}/${CURNAME}/${CURNAME}.root
		STDOUT_FNAME=${OUT_DIR}/${CURNAME}/${CURNAME}.stdout
		STDERR_FNAME=${OUT_DIR}/${CURNAME}/${CURNAME}.stderr
		SCRIPT=${OUT_DIR}/${CURNAME}/tPDF_${CURNAME}.sh
		
		# Copy and modify WCSim macro file
		cp $BASE_WCSIMMACRO ${THIS_MACRO}
		
		# Calculate kinetic energy
		
		KE=$( echo 'sqrt('${MOM}'^2 + '${PARTICLE_MASS}'^2) - '${PARTICLE_MASS} | bc -l )
		
		sed -i 's:^CHANGE_ME_GEOM*:/WCSim/WCgeom '${DETECTOR_NAME}':' $THIS_MACRO
		if [ ${DETECTOR_NAME} = "nuPRISM" ]
		    then
                    HEIGHT=$( echo ${DETECTOR_H}'*2/100' | bc -l )
                    DIAMETER=$( echo ${DETECTOR_R}'*2/100' | bc -l )
		    sed -i 's:^CHANGE_ME_PMT*:/WCSim/nuPRISM/SetPMTType '${PMT_NAME}':' $THIS_MACRO
		    sed -i 's:^CHANGE_ME_HEIGHT*:/WCSim/nuPRISM/SetDetectorHeight '${HEIGHT}' m:' $THIS_MACRO
		    sed -i 's:^CHANGE_ME_DIAMETER*:/WCSim/nuPRISM/SetDetectorDiameter '${DIAMETER}' m:' $THIS_MACRO
		fi

		sed -i 's:^CHANGE_ME_ROOT_FILE_NAME*:/WCSimIO/RootFile '${ROOT_FNAME}':' $THIS_MACRO
		sed -i 's:^CHANGE_ME_N_EVTS*:/run/beamOn '${N_EVENTS}':' $THIS_MACRO
		sed -i 's:^CHANGE_ME_PARTICLE*:/gps/particle '${PARTICLE_NAME}':' $THIS_MACRO
		sed -i 's:^CHANGE_ME_ENERGY*:/gps/energy '${KE}' MeV:' $THIS_MACRO
		sed -i 's:^CHANGE_ME_SEEDS*:/WCSim/random/seed '${SEED1}':' $THIS_MACRO
		sed -i 's:^CHANGE_ME_VTX_HALFZ*:/gps/pos/halfz '${DETECTOR_H}' cm:' $THIS_MACRO
		sed -i 's:^CHANGE_ME_VTX_RADIUS*:/gps/pos/radius '${DETECTOR_R}' cm:' $THIS_MACRO
#		sed -i 's:^CHANGE_ME_SEEDS*:/random/setSeeds '${SEED1}' '${SEED2}':' $THIS_MACRO

		# Write bash script
		echo '#!/bin/bash'              >  ${SCRIPT}
		echo 'source '${SETUP_SCRIPT}   >> ${SCRIPT}

    echo 'source /home/tyoshida/nuPRISM/geant4.10.1/geant4.10.01.p03-install/bin/geant4.sh' >> ${SCRIPT}
		echo 'source /home/tyoshida/nuPRISM/geant4.10.1/geant4.10.01.p03-build/geant4make.sh' >> ${SCRIPT}

		echo 'cd '${WCSIM_BASE_DIR}     >> ${SCRIPT}
		echo ${WCSIM}' '${THIS_MACRO}   >> ${SCRIPT}
		echo 'cd -'                     >> ${SCRIPT}
		echo 'echo "done"'              >> ${SCRIPT}
		
		cd ${OUT_DIR}/${CURNAME}
		
		echo ' -r job'${CURNAME}' -o '${STDOUT_FNAME}' -e '${STDERR_FNAME} ${SCRIPT} >> /home/tyoshida/jobManager/${QUEUE_NAME}.list

#		qsub -q ${QUEUE_NAME} -o ${STDOUT_DIR} -e ${STDERR_DIR} $SCRIPT
#		qsub -q ${QUEUE_NAME} -o /dev/null $SCRIPT
		
		cd -
	    done # Ends run loop
	done # Ends job loop
    done < $CHART_FILE # End momentum loop
done # End particle loop
