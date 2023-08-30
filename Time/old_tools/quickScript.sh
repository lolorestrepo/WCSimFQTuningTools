#!/bin/bash
WCSIM_CONFIG=NuPRISM_10x8_8inchPMT_40perCent

#P_NAMES=(e- mu- pi+)
#PDG_CODES=(11 13 211)
#P_NAMES=(e- mu-)
#PDG_CODES=(11 13)
#P_NAMES=(e-)
#PDG_CODES=(11)
#P_NAMES=(mu-)
#PDG_CODES=(13)
P_NAMES=(pi+)
PDG_CODES=(211)

rm *_hist_sum.root
echo rm *_hist_sum.root
for i in `ls /disk01/usr5/tyoshida/fiTQunTuning/${WCSIM_CONFIG}_tbugfix/timepdf/*_hist_sum.root`
do
    ln -s $i .
    echo    ln -s $i .
done
#exit

for ((i=0; i<${#P_NAMES[*]}; i++))
do
    echo 'Starting ' ${P_NAMES[$i]} ${PDG_CODES[$i]} 
    root -l -b -q combhists.cc'('${PDG_CODES[$i]}')' 2> /dev/null
    echo 'doing ' root -b -q fittpdf.cc'('${PDG_CODES[$i]},0,1')'
#    root -l -b -q fittpdf.cc'('${PDG_CODES[$i]},0,1')'
    root -l -b -q fittpdf_test.cc'('${PDG_CODES[$i]},0,1')'
done

#rm -r $WCSIM_CONFIG
mkdir -p $WCSIM_CONFIG
mv *_tpdf*.root ${WCSIM_CONFIG}/
