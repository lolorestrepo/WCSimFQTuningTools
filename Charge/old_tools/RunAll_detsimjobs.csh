#!/bin/csh
# $1: SK version [1-4]

set wd = `pwd`

if ($#argv == 0) then
echo "Specify SK version!"
exit
endif

echo "SK " $1

set prodver = ""

if ($#argv > 1) then
set prodver = $2
endif

mkdir -p out_sk$1$prodver
cd out_sk$1$prodver

$wd/detsimjob/job_1dpdf.csh $wd/detsimjob/pt0.txt $wd/detsimjob/nodns_sk$1_80.card >& log_pt0.txt &
$wd/detsimjob/job_1dpdf.csh $wd/detsimjob/pt1.txt $wd/detsimjob/nodns_sk$1_40.card >& log_pt1.txt &
$wd/detsimjob/job_1dpdf.csh $wd/detsimjob/pt2.txt $wd/detsimjob/nodns_sk$1_20.card >& log_pt2.txt &
$wd/detsimjob/job_1dpdf.csh $wd/detsimjob/pt3.txt $wd/detsimjob/nodns_sk$1_20.card >& log_pt3.txt &
$wd/detsimjob/job_1dpdf.csh $wd/detsimjob/pt4.txt $wd/detsimjob/nodns_sk$1_20.card >& log_pt4.txt &
$wd/detsimjob/job_1dpdf.csh $wd/detsimjob/pt5.txt $wd/detsimjob/nodns_sk$1_20.card >& log_pt5.txt &
$wd/detsimjob/job_1dpdf.csh $wd/detsimjob/pt6.txt $wd/detsimjob/nodns_sk$1_20.card >& log_pt6.txt &
