#!/bin/bash

# Run in the directory where the tarballs are!

if [ $# != 2 ]
then
  echo "Usage: ./job_gentpdf.sh [PDG] [SK Ver]"
  exit
fi

echo "SK " $2

export PDG=$1

export TPDF_TOOLDIR="/neut/data20/shimpeit/timepdf/timepdf/"

echo "PDG = " $PDG

ln -s $TPDF_TOOLDIR/makehist

#$TPDF_TOOLDIR/mergehists.pl $PDG

cd outdir

root -b -q ${TPDF_TOOLDIR}/combhists.cc\(${PDG}\)

root -b -q ${TPDF_TOOLDIR}/fittpdf.cc\(${PDG},0,1\)

mv ${PDG}_tpdfpar.root ${PDG}_tpdfpar_sk$2.root

ln -s ${PDG}_tpdfpar_sk$2.root ${PDG}_tpdfpar.root

cd ../
