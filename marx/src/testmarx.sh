#! /bin/sh
#set -v
MARX_DATA_DIR=../data
export MARX_DATA_DIR

use_valgrind=1

output_dir="/tmp/marx"
mkdir $output_dir
mkdir /tmp/marx/uparm
valgrind_dir=$output_dir/valgrind
mkdir $valgrind_dir
valgrind="valgrind --leak-check=yes --logfile=$valgrind_dir"

PFILES=/tmp/marx/uparm

UPARM=$PFILES
export PFILES UPARM
cp ../par/marxasp.par ../par/marx.par ../par/pileup.par /tmp/marx/uparm 

marx="./$ARCH""objs/marx NumRays=1000 dNumRays=200 ExposureTime=0 DitherModel=INTERNAL"
marx2fits="./$ARCH""objs/marx2fits"
marxasp="./$ARCH""objs/marxasp"
pileup="./$ARCH""objs/pileup"

if [ "$use_valgrind" = "1" ]
then
  valgrindmarx="$valgrind/marx"
  valgrindmarx2fits="$valgrind/marx2fits"
  valgrindmarxasp="$valgrind/marxasp"
  valgrindpileup="$valgrind/pileup"
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
         $valgrindmarx $marx OutputDir=$output_dir/marx$count DetectorType=$det GratingType=$grat UseGratingEffFiles=$x
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

$valgrindpileup $pileup MarxOutputDir="$output_dir"/marx1 Verbose=1
$valgrindmarx2fits $marx2fits --pileup "$output_dir"/marx1/pileup "$output_dir"/marx1/pileup.fits

exit 0
