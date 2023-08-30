#! /usr/bin/perl

use strict;
use warnings;
use Cwd;

my $wd = Cwd::getcwd();

my $hstcmd = "${wd}/makehist";# set the location of "makehist"

my $ipart=$ARGV[0];

my $fPCflg=1;

my $nwtr=1.38;

print "nwtr=", $nwtr, "\n";

my $nfls;
my $fnbase;

system("date");

print $wd, "\n";

system("mkdir -p ${wd}/tmp");
my @list = <job_tpdf_sk?_${ipart}_*_out.tar>;
$nfls=@list;
if ($nfls > 0) {
  chdir("${wd}/tmp");
  for (my $ifile=0; $ifile<$nfls; $ifile++) {
    system("echo ../$list[$ifile]");
    system("tar -xvf ../$list[$ifile]");
  }
}

chdir("${wd}");
system("mkdir -p ${wd}/outdir");
chdir("${wd}/tmp");
for (my $mom=1; $mom<11000; $mom++) {
  my @datlst = <${ipart}_${mom}_*.dat>;
  $nfls=@datlst;
  if ($nfls > 0) {#available momentum!
    print "${mom}\n";
    for (my $ifile=0; $ifile<$nfls; $ifile++) {
      system("echo $datlst[$ifile]");
      system("${hstcmd} $datlst[$ifile] ${nwtr} ${fPCflg}");
    }
    system("hadd ${ipart}_${mom}_hist_sum.root ${ipart}_${mom}_*_hist.root");
    system("mv ${ipart}_${mom}_hist_sum.root ${wd}/outdir/");
  }
}

system("date");
