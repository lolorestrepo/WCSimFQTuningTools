#!/bin/bash


SETUP_SCRIPT=/home/tyoshida/nuPRISM/Analysis/Source_At_Start_nuPRISM.sh
WCSIM=/home/tyoshida/nuPRISM/WCSim/exe/bin/Linux-g++/WCSim

DETECTOR_NAME=NuPRISM_10x8_8inchPMT_40perCent


if [ $# -eq 0 ]
then
    echo "No arguments supplied, detector name is "${DETECTOR_NAME}
    echo "PMT type is "${PMT_NAME}
    exit
elif [ $# -eq 1 ]
then
    DETECTOR_NAME=$1
    echo "Detector name given as an argument: "${DETECTOR_NAME}
fi

OUT_DIR=/disk01/usr5/tyoshida/fiTQunTuning/${DETECTOR_NAME}_tbugfix/timepdf
TIMEPDF_DIR=/home/tyoshida/fiTQunTuning/Utilities/timepdf


source $SETUP_SCRIPT

echo "Script starts"
## List of particles

#P_NAMES=(e- mu- pi+)
#PDG_CODES=(11 13 211)
#P_MASSES=(0.511 105.7 139.6)
#PROD_NUM=0
#P_NAMES=(e- mu-)
#PDG_CODES=(11 13)
#P_MASSES=(0.511 105.7)
#PROD_NUM=0

#P_NAMES=(mu- pi+)
#PDG_CODES=(13 211)
#P_MASSES=(105.7 139.6)
#PROD_NUM=0
#P_NAMES=(e-)
#PDG_CODES=(11)
#P_MASSES=(0.511)
#PROD_NUM=0

#P_NAMES=(pi+)
#PDG_CODES=(211)
#P_MASSES=(139.6)
#PROD_NUM=0

#P_NAMES=(mu-)
#PDG_CODES=(13)
#P_MASSES=(105.7)
#PROD_NUM=0

P_NAMES=(pi+)
PDG_CODES=(211)
P_MASSES=(139.6)
PROD_NUM=0


echo "Number of particle types " ${#P_NAMES[*]} 

FAILED_JOBS=()
N_FAILED_JOBS=0

WRONG_EVENTS=()
N_WRONG_EVENTS=0

# The following commented out block will check that all root files exist and have
# the correct number of entries, etc. It takes a long time to run...
#--------------------------------------------------------------------------------
# Particle loop

#   for ((i=0; i<${#P_NAMES[*]}; i++))
#   do
#       echo 'Starting ' ${P_NAMES[$i]}
#   
#       PDG_CODE=${PDG_CODES[$i]}
#       PARTICLE_NAME=${P_NAMES[$i]}
#       PARTICLE_MASS=${P_MASSES[$i]}
#   
#       CHART_FILE=${TIMEPDF_DIR}/chart_${PDG_CODE}.txt
#       
#       ls ${CHART_FILE}
#       
#       COPY_NUMBER=0
#       MOM=-1
#       
#       # Momentum Loop
#       while read line 
#       do
#   	echo "Processing line " $line
#   	read -a lineArr <<< $line
#   	
#   	# Number of runs in a specific job
#   	N_RUNS="${lineArr[0]}"
#   	N_JOBS="${lineArr[1]}"
#   	N_EVENTS=0
#   	N_EVENTS="${lineArr[2]}"
#   	
#   	echo "Runs" $N_RUNS "Jobs" $N_JOBS "Events" $N_EVENTS
#   	
#   	# Job loop
#   	for job in `seq 0 $(($N_JOBS-1))` 
#   	do
#   	    # Run Loop
#   	    #	echo " Processing job" $job 
#   	    for run in `seq 0 $(($N_RUNS-1))`
#   	    do
#   		if [ $job -eq 0  ] 
#   		then
#   		    echo '++++++++++++++++++= is zeroth job'
#   		    if [  $MOM -eq ${lineArr[$((3+$run))]} ]
#   		    then
#   			COPY_NUMBER=$(( $COPY_NUMBER + 1 ))
#   			echo "++++++++++++++++++++++ COPY NUMBER INCREMENTED " ${COPY_NUMBER}
#   		    else
#   			COPY_NUMBER=0
#   		    fi
#   		fi
#   		
#   		MOM="${lineArr[$((3+$run))]}"
#   
#   		# Check Time PDF file exists
#   		CURNAME=${PDG_CODE}_${MOM}_${COPY_NUMBER}_${job}_${PROD_NUM}
#   		ROOT_FNAME=${OUT_DIR}/${CURNAME}/${CURNAME}.root
#   		
#   		if [ ! -f $ROOT_FNAME ]; then
#   		    echo "File " $ROOT_FNAME " not found."
#   
#   		    FAILED_JOBS[$N_FAILED_JOBS]="${ROOT_FNAME}"
#   		    ((N_FAILED_JOBS++))
#   
#   		else 
#   		    N_ENTRIES=`/home/cvilela/fiTQunTuning/TimePDF/Utilities/timepdf/SukapSampleScripts/wcSimGetEntries.py  $ROOT_FNAME`
#   		    echo ${N_EVENTS}END ${N_ENTRIES}END
#   
#   		    
#   		    if [ "$N_EVENTS" -eq "$N_ENTRIES" ]
#   #		    if (( "$N_EVENTS" == "$N_ENTRIES" ))
#   		    
#    		    then
#    			echo "Correct number of events."
#    		    else
#    			WRONG_EVENTS[$N_WRONG_EVENTS]="${ROOT_FNAME}"
#    			((N_WRONG_EVENTS++))
#    		    fi # Number of events in file is the same?
#    		fi # WCSim output file exists?
#   	    done # Ends run loop
#   	done # Ends job loop
#       done < $CHART_FILE # End momentum loop
#   done # End particle loop

#--------------------------------------------------------------------------------


# Check if all output files exists, if not print out which jobs failed
if [ ! ${#FAILED_JOBS[*]} -eq 0 ] 
then 
    echo "${#FAILED_JOBS[*]} Failed jobs: "
    for failed_job in  ${FAILED_JOBS[*]}
    do
	echo $failed_job
    done
fi

if [ ! ${#WRONG_EVENTS[*]} -eq 0 ] 
then 
    echo "${#WRONG_EVENTS[*]} jobs with the wrong number of events: "
    for failed_job in  ${WRONG_EVENTS[*]}
    do
	echo $failed_job
    done
fi

if [ ! ${#FAILED_JOBS[*]} -eq 0 ]; then
    exit
fi
if [ ! ${#WRONG_EVENTS[*]} -eq 0 ]
then
    exit
fi

echo "All WCSim jobs OK"

# Particle loop
for ((i=0; i<${#P_NAMES[*]}; i++))
do
    echo 'Starting ' ${P_NAMES[$i]}

    MOMS_LIST=()
    N_MOMS=0

    PDG_CODE=${PDG_CODES[$i]}
    PARTICLE_NAME=${P_NAMES[$i]}
    PARTICLE_MASS=${P_MASSES[$i]}

    CHART_FILE=${TIMEPDF_DIR}/chart_${PDG_CODE}.txt
    
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
	N_EVENTS=0
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
		
		MOM="${lineArr[$((3+$run))]}"

		if [[ `printf "%s\n" "${MOMS_LIST[@]}" | grep "^${MOM}$"` ]] 
		then
		    echo 'Have seen this MOM'
		else
		    MOMS_LIST[$N_MOMS]="${MOM}"
		    ((N_MOMS++))
		fi

		CURNAME=${PDG_CODE}_${MOM}_${COPY_NUMBER}_${job}_${PROD_NUM}
		ROOT_FNAME=${OUT_DIR}/${CURNAME}/${CURNAME}.root
		
		echo "Making TimePDF plot"
#		cp ${FITQUN_ROOT}/fiTQun.parameters.dat ${FITQUN_ROOT}/fiTQun.parameters.dat.backup.${DETECTOR_NAME}_tuning
#		cat ${FITQUN_ROOT}/ParameterOverrideFiles/${DETECTOR_NAME}.parameters.dat >> ${FITQUN_ROOT}/fiTQun.parameters.dat

		$TIMEPDF_DIR/makehistWCSim $ROOT_FNAME ${FITQUN_ROOT}/ParameterOverrideFiles/${DETECTOR_NAME}.parameters.dat
#		echo $TIMEPDF_DIR/checkTResVsL $ROOT_FNAME ${FITQUN_ROOT}/ParameterOverrideFiles/${DETECTOR_NAME}.parameters.dat
#		$TIMEPDF_DIR/checkTResVsL $ROOT_FNAME ${FITQUN_ROOT}/ParameterOverrideFiles/${DETECTOR_NAME}.parameters.dat


#		mv ${FITQUN_ROOT}/fiTQun.parameters.dat.backup.${DETECTOR_NAME}_tuning ${FITQUN_ROOT}/fiTQun.parameters.dat
	    done # Ends run loop
	done # Ends job loop
    done < $CHART_FILE # End momentum loop

    for i_mom in ${MOMS_LIST[*]}
    do
	echo $i_mom
	rm ${OUT_DIR}/${PDG_CODE}_${i_mom}_hist_sum.root
	hadd ${OUT_DIR}/${PDG_CODE}_${i_mom}_hist_sum.root `ls ${OUT_DIR}/${PDG_CODE}_${i_mom}_?_?_?/${PDG_CODE}_${i_mom}_?_?_?_hist.root`
#	rm ${OUT_DIR}/${PDG_CODE}_${i_mom}_checkTResVsL_sum.root
#	hadd ${OUT_DIR}/${PDG_CODE}_${i_mom}_checkTResVsL_sum.root `ls ${OUT_DIR}/${PDG_CODE}_${i_mom}_?_?_?/${PDG_CODE}_${i_mom}_?_?_?_checkTResL.root`
    done
done # End particle loop
