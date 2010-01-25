#! /bin/sh
#set -v
MARX_DATA_DIR=../data; export MARX_DATA_DIR

use_valgrind=0

#bindir=/tmp/marxroot/bin
bindir=./${ARCH}objs
#bindir=/nfs/cxc/a1/i686/opt/packages/marx-4.3.0/bin
#bindir=$HOME/sys/linux/test/bin

output_dir="/tmp/marxout"
mkdir $output_dir
mkdir $output_dir/uparm
valgrind_dir=$output_dir/valgrind
mkdir $valgrind_dir
valgrind="valgrind --tool=memcheck --leak-check=yes --error-limit=no --num-callers=25 --log-file=$valgrind_dir"

PFILES=$output_dir/uparm

UPARM=$PFILES
export PFILES UPARM
cp ../par/marxasp.par ../par/marx.par ../par/marxpileup.par $output_dir/uparm 

marx="$bindir/marx NumRays=1000 dNumRays=200 ExposureTime=0 DitherModel=INTERNAL"
marx2fits="$bindir/marx2fits"
marxasp="$bindir/marxasp"
marxpileup="$bindir/marxpileup"

if [ "$use_valgrind" = "1" ]
then
  valgrindmarx="$valgrind/marx"
  valgrindmarx2fits="$valgrind/marx2fits"
  valgrindmarxasp="$valgrind/marxasp"
  valgrindpileup="$valgrind/marxpileup"
else
  valgrindmarx=""
  valgrindmarx2fits=""
  valgrindmarxasp=""
  valgrindpileup=""
fi

#purify="./solarisobjs/purify.log"

detectors="ACIS-S ACIS-I HRC-S HRC-I"
gratings="NONE HETG LETG"
count=1
#VALGRIND="valgrind -v --leak-check=yes"
for det in $detectors
do
   for grat in $gratings
   do
      for x in yes no
      do
         #$valgrindmarx $marx OutputDir=$output_dir/marx$count DetectorType=$det GratingType=$grat UseGratingEffFiles=$x
         gdb --args $marx OutputDir=$output_dir/marx$count DetectorType=$det GratingType=$grat UseGratingEffFiles=$x
	 exit 1
	 if [ "$?" != "0" ]
	 then
	    if [ $use_valgrind = 0 ]; then exit 1; fi
	 fi
#	 /bin/cp $purify $output_dir/marx$count/marx.purify
	 $valgrindmarx2fits $marx2fits $output_dir/marx$count $output_dir/marx$count/events.fits
	 if [ "$?" != "0" ]
	 then
	    if [ $use_valgrind = 0 ]; then exit 1; fi
	 fi
#	 /bin/cp $purify $output_dir/marx$count/marx2fits.purify
	 $valgrindmarxasp $marxasp MarxDir=$output_dir/marx$count OutputFile=$output_dir/marx$count/pcad.fits
	 if [ "$?" != "0" ]
	 then
	    if [ $use_valgrind = 0 ]; then exit 1; fi
	 fi
#	 /bin/cp $purify $output_dir/marx$count/marxasp.purify
        count=`expr $count + 1`
      done
   done
done

$valgrindpileup $marxpileup MarxOutputDir="$output_dir"/marx1 Verbose=1
$valgrindmarx2fits $marx2fits --pileup "$output_dir"/marx1/pileup "$output_dir"/marx1/pileup.fits

exit 0
