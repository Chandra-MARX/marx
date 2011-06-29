#! /bin/sh
#set -v
#MARX_DATA_DIR=../data; export MARX_DATA_DIR

use_valgrind=1
use_gdb=0

#bindir=/tmp/marxroot/bin
bindir=./objs
#bindir=/nfs/cxc/a1/i686/opt/packages/marx-4.3.0/bin
#bindir=$HOME/sys/linux/opt/marx4.5/bin
#bindir=/nfs/cxc/h1/davis/sys/linux/test/bin
#bindir=/nfs/cxc/h1/davis/sys/linux/test/bin
#bindir=/nfs/cxc/h1/davis/sys/linux-x86_64/opt/marx5.0/bin

output_dir="/tmp/marxout1"
mkdir $output_dir
mkdir $output_dir/uparm
valgrind_dir=$output_dir/valgrind
mkdir $valgrind_dir
valgrind="valgrind --tool=memcheck --leak-check=yes --error-limit=no --num-callers=25 --log-file=$valgrind_dir/%q{VGOUTPUT}"
if [ $use_valgrind = 0 ] 
then
  unset valgrind
fi

if [ $use_gdb = 1 ]
then
  valgrind="gdb --args"
fi

PFILES=$output_dir/uparm

UPARM=$PFILES
export PFILES UPARM
cp ../par/marxasp.par ../par/marx.par ../par/marxpileup.par $output_dir/uparm 

marx="$bindir/marx NumRays=1000 dNumRays=200 ExposureTime=0 \
DitherModel=INTERNAL SourceRA=45 SourceDEC=30 RA_Nom=45 Dec_Nom=30 Roll_Nom=20 \
OutputVectors=ETXYZ123DxyMPOabcdSrB#
"
marx2fits="$bindir/marx2fits"
marxasp="$bindir/marxasp"
marxpileup="$bindir/marxpileup"

#purify="./solarisobjs/purify.log"

detectors="ACIS-S ACIS-I HRC-S HRC-I"
gratings="NONE HETG LETG"
count=1
#VALGRIND="valgrind -v --leak-check=yes"
#valgrind="gdb --args"
#use_valgrind=0
for det in $detectors
do
   for grat in $gratings
   do
      for x in yes no
      do
         export VGOUTPUT="marx-$count.log"
         $valgrind $marx OutputDir=$output_dir/marx$count DetectorType=$det GratingType=$grat UseGratingEffFiles=$x
         #$marx OutputDir=$output_dir/marx$count DetectorType=$det GratingType=$grat UseGratingEffFiles=$x
	 if [ "$?" != "0" ]
	 then
	    if [ $use_valgrind = 0 ]; then exit 1; fi
	 fi

#	 /bin/cp $purify $output_dir/marx$count/marx.purify

         export VGOUTPUT="marx2fits-$count.log"
	 $valgrind $marx2fits $output_dir/marx$count $output_dir/marx$count/events.fits
	 if [ "$?" != "0" ]
	 then
	    if [ $use_valgrind = 0 ]; then exit 1; fi
	 fi
#	 /bin/cp $purify $output_dir/marx$count/marx2fits.purify

         export VGOUTPUT="marxasp-$count.log"
	 $valgrind $marxasp MarxDir=$output_dir/marx$count OutputFile=$output_dir/marx$count/pcad.fits
	 if [ "$?" != "0" ]
	 then
	    if [ $use_valgrind = 0 ]; then exit 1; fi
	 fi
#	 /bin/cp $purify $output_dir/marx$count/marxasp.purify
        count=`expr $count + 1`
	 if [ $use_gdb != "0" ]; then exit 0; fi
	if [ $count = "xxxxx" ]
	then
	   valgrind="gdb --args"
	   use_gdb=1
	fi
      done
   done
done

export VGOUTPUT="marxpileup.log"
$valgrind $marxpileup MarxOutputDir="$output_dir"/marx1 Verbose=1
export VGOUTPUT="marx2fits-pileup.log"
$valgrind $marx2fits --pileup "$output_dir"/marx1/pileup "$output_dir"/marx1/pileup.fits

exit 0
