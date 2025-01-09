#!/bin/bash
#source perl-5.30.0
#source gcc-5.2.0
#source python-3.8.3
#source interproscan-5.45
source nlrtracker-1.0.0


dir=`pwd`
fasta=$1
tmpdir=$2
threads=$(($3 - 2))
out=$4


interproscan-5.72-103.0.sh -i $fasta \
  -f gff3 \
  -T $tmpdir \
  -cpu $threads \
  -o $out \
  -appl Pfam,Gene3D,SUPERFAMILY,PRINTS,SMART,CDD,ProSiteProfiles
