#! /usr/bin/perl

use strict;
use warnings;

my $debug=0;

my $ipart=13;# PDG code

my $SKVer=3;# Set SK version!

my $prodnm=0;

my $processlist="chart_${ipart}.txt";

my $mom;

my $jobname;

my $scratchdir="/scratch/t/tanaka/shimpeit/timepdf/sk${SKVer}/${ipart}_${prodnm}";

my $ramworkdir;

my $filestotar;

system("mkdir -p $scratchdir");

# create a universal skdetsim script for thir production
open (OLD,"$ENV{SKDETSIM_ROOT}/skdetsim.sh");
open (NEW,">${scratchdir}/skdetsim.sh");
while(<OLD>) {
  s/set DIR.*/set DIR=\$SKDETSIM_ROOT/;
  if (-e "$ENV{SKDETSIM_ROOT}/skdetsim_high") {
    s/\$DIR\/skdetsim /\$DIR\/skdetsim_high /;
  }
  print NEW;
}
close (NEW);
close (OLD);
chmod 0755,"${scratchdir}/skdetsim.sh";


open (IN, $processlist);

my $nthreads=0;

my $momtmp=0;
my $imom=0;

while (<IN>) {
  
  chomp ($_);
  
  my @data = split(/\t+/,$_);
  
  for (my $i=0; $i<@data; $i++) {
    $data[$i] =~ s/(　| )+//g;
  }
  
  if ($nthreads == 0) {#new job
    if ($momtmp == $data[3]) {
      $imom=$imom + 1;
    }
    else {
      $imom=0;
    }
    $momtmp = $data[3];
    
    my $walltime=int($data[7]*1.5/3600+1);
    $jobname="job_tpdf_sk${SKVer}_${ipart}_${momtmp}_${imom}_${prodnm}";
    $ramworkdir="/dev/shm/${jobname}";
    open (SHS,">${scratchdir}/${jobname}.sh");
    
    # setup common commands for all 8 threads
    print SHS "#!/bin/bash

#PBS -l nodes=1:ppn=8,walltime=$walltime:00:00
#PBS -N ${jobname}

date

cd /home/t/tanaka/shimpeit
source setup_head_neut5142.sh

mkdir ${ramworkdir}
cd ${ramworkdir}

cp ${scratchdir}/${jobname}.tar ./

tar -xf ${jobname}.tar
rm ${jobname}.tar
";
    $filestotar="skdetsim.sh";
  }
  
  $nthreads=$nthreads+$data[1];
  
  print "$data[0], $data[1], $data[2], $data[3], $data[4], $data[5], $data[6], $data[7]\n";
  
  my $nevts=$data[2];
  
  for (my $i=0; $i<$data[1]; $i++) {
    #operations in one thread
    
    print SHS "
################################################
\(
";
    
    for (my $j=0; $j<$data[0]; $j++) {#sequential ops.
      $mom=$data[$j+3];#momentum
      print "$mom, ";
      
      my $curname="${ipart}_${mom}_${imom}_${i}_${prodnm}";
      
      ### create files for each subjob
      open (OLDCARD,"tpdf_sk${SKVer}.card");
      open (NEWCARD,">${scratchdir}/${curname}.card");
      
      my @rand = ();
      foreach my $ir (0..4) {
        $rand[$ir] = int(rand(2147483646));
        if ($ir >= 3) {
          $rand[$ir] = 0;
        }
      }
      
      while (<OLDCARD>) {
        s/VECT-NEVT.*/VECT-NEVT $nevts/;
        s/VECT-RAND.*/VECT-RAND $rand[0] $rand[1] $rand[2] $rand[3] $rand[4]/;
        s/VECT-PART.*/VECT-PART $ipart/;
        s/VECT-MOM.*/VECT-MOM $mom/;
        print NEWCARD;
      }
      
      close (NEWCARD);
      close (OLDCARD);
      
      $filestotar.=" ${curname}.card";
      
      print SHS "./skdetsim.sh ${curname}.card ${curname}.dat
";
      
    }
    print SHS 
"\)&
";
  }
  
  
  if ($nthreads == 8) {
    #change node
    
    # finish up
    print SHS "
################################################
wait
cd $ramworkdir
tar -czf ${scratchdir}/${jobname}_out.tar *.dat
cd
rm -rf ${ramworkdir}
date
";
    
    #        print "\$ijob = $ijob, closing sh file\n";
    close SHS;
    
    chmod 0755,"${scratchdir}/${jobname}.sh";
    
    system("cd ${scratchdir}; tar -cf ${jobname}.tar ${filestotar}");
    
    print "\n";
    $nthreads=0;
    
    # submit the job
    if ($debug) {
	    system("qsub ${scratchdir}/${jobname}.sh -q debug -j oe -o ${scratchdir}/${jobname}.log");
    } else {
	    system("qsub ${scratchdir}/${jobname}.sh -j oe -o ${scratchdir}/${jobname}.log");
    }
  }
}

close (IN);
