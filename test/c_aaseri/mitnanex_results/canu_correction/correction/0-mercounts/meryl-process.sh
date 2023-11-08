#!/bin/sh


#  Path to Canu.

bin="/Users/jjpc/miniforge3/envs/mitnanex/bin"

#  Report paths.

echo ""
echo "Found perl:"
echo "  " `which perl`
echo "  " `perl --version | grep version`
echo ""
echo "Found java:"
echo "  " `which /Users/jjpc/miniforge3/envs/mitnanex/lib/jvm/bin/java`
echo "  " `/Users/jjpc/miniforge3/envs/mitnanex/lib/jvm/bin/java -showversion 2>&1 | head -n 1`
echo ""
echo "Found canu:"
echo "  " $bin/canu
echo "  " `$bin/canu -version`
echo ""


#  Environment for any object storage.

export CANU_OBJECT_STORE_CLIENT=
export CANU_OBJECT_STORE_CLIENT_UA=
export CANU_OBJECT_STORE_CLIENT_DA=
export CANU_OBJECT_STORE_NAMESPACE=
export CANU_OBJECT_STORE_PROJECT=




#  Discover the job ID to run, from either a grid environment variable and a
#  command line offset, or directly from the command line.
#
if [ x$CANU_LOCAL_JOB_ID = x -o x$CANU_LOCAL_JOB_ID = xundefined -o x$CANU_LOCAL_JOB_ID = x0 ]; then
  baseid=$1
  offset=0
else
  baseid=$CANU_LOCAL_JOB_ID
  offset=$1
fi
if [ x$offset = x ]; then
  offset=0
fi
if [ x$baseid = x ]; then
  echo Error: I need CANU_LOCAL_JOB_ID set, or a job index on the command line.
  exit
fi
jobid=`expr -- $baseid + $offset`
if [ x$baseid = x0 ]; then
  echo Error: jobid 0 is invalid\; I need CANU_LOCAL_JOB_ID set, or a job index on the command line.
  exit
fi
if [ x$CANU_LOCAL_JOB_ID = x ]; then
  echo Running job $jobid based on command line options.
else
  echo Running job $jobid based on CANU_LOCAL_JOB_ID=$CANU_LOCAL_JOB_ID and offset=$offset.
fi

if [ $jobid -gt 1 ]; then
  echo Error: Only 1 job, you asked for $jobid.
  exit 1
fi

#  If the meryl ignore files exst, then we're done.

if [ -e ./c_aaseri_ERR10466724_collected_reads.ms16.histogram -a -e ./c_aaseri_ERR10466724_collected_reads.ms16.dump -a -e ./c_aaseri_ERR10466724_collected_reads.ms16.ignore.gz ] ; then
  exit 0
fi

#  If those exist in the object store, we're also done.


if [ -e c_aaseri_ERR10466724_collected_reads.ms16.histogram ]; then
  exists1=true
else
  exists1=false
fi

if [ -e c_aaseri_ERR10466724_collected_reads.ms16.dump ]; then
  exists2=true
else
  exists2=false
fi

if [ -e c_aaseri_ERR10466724_collected_reads.ms16.ignore.gz ]; then
  exists3=true
else
  exists3=false
fi
if [ $exists1 = true -a $exists2 = true -a $exists3 = true ] ; then
  echo "Output files 'c_aaseri_ERR10466724_collected_reads.ms16.histogram', 'c_aaseri_ERR10466724_collected_reads.ms16.dump' and 'c_aaseri_ERR10466724_collected_reads.ms16.ignore.gz' exist in 'correction/0-mercounts'."
  exit 0
fi


#  Nope, not done.  Fetch all the intermediate meryl databases.


#
#  Merge counting jobs, strip out unique kmers.
#

if [ ! -e ./c_aaseri_ERR10466724_collected_reads.ms16/merylIndex ] ; then
  /Users/jjpc/miniforge3/envs/mitnanex/bin/meryl threads=4 memory=2 \
    greater-than 1 \
      output c_aaseri_ERR10466724_collected_reads.ms16.WORKING \
      union-sum  \
        ./c_aaseri_ERR10466724_collected_reads.01.meryl \
  && \
  mv -f ./c_aaseri_ERR10466724_collected_reads.ms16.WORKING ./c_aaseri_ERR10466724_collected_reads.ms16

  #  Fail if there is no meryl database.
  if [ ! -e ./c_aaseri_ERR10466724_collected_reads.ms16/merylIndex ] ; then
    echo meryl merge failed.
    exit 1
  fi

  #  Remove meryl intermediate files.
  rm -rf ./c_aaseri_ERR10466724_collected_reads.01.meryl ./c_aaseri_ERR10466724_collected_reads.01.meryl.err
fi

#
#  Dump a histogram, 'cause they're useful.
#

if [ ! -e ./c_aaseri_ERR10466724_collected_reads.ms16.histogram ] ; then
  /Users/jjpc/miniforge3/envs/mitnanex/bin/meryl threads=1 memory=1 \
    statistics ./c_aaseri_ERR10466724_collected_reads.ms16 \
  > ./c_aaseri_ERR10466724_collected_reads.ms16.histogram.WORKING \
  && \
  mv ./c_aaseri_ERR10466724_collected_reads.ms16.histogram.WORKING ./c_aaseri_ERR10466724_collected_reads.ms16.histogram
fi

#
#  Dump frequent mers.
#
#  The indenting of the at-least options is misleading.  'print'
#  takes input from the first 'at-least', which that takes input from
#  the second 'at-least'.  The effect is the same as taking the
#  'intersection' of all the 'at-least' filters -- logically, it is
#  doing 'at-least X AND at-least Y AND at-least Z'.
#

if [ ! -e ./c_aaseri_ERR10466724_collected_reads.ms16.dump ] ; then
  /Users/jjpc/miniforge3/envs/mitnanex/bin/meryl threads=4 memory=2 \
    print ./c_aaseri_ERR10466724_collected_reads.ms16.##.dump \
      at-least threshold=1000 \
      at-least word-frequency=0.0000001 \
        ./c_aaseri_ERR10466724_collected_reads.ms16

  cat ./c_aaseri_ERR10466724_collected_reads.ms16.??.dump > ./c_aaseri_ERR10466724_collected_reads.ms16.dump
  rm -f ./c_aaseri_ERR10466724_collected_reads.ms16.??.dump
fi

#
#  Convert the dumped kmers into a mhap ignore list.
#
#    numKmers - number of kmers we're filtering
#    totKmers - total number of kmers in the dataset

if [ ! -e ./c_aaseri_ERR10466724_collected_reads.ms16.ignore.gz ] ; then
  numKmers=`wc -l < ./c_aaseri_ERR10466724_collected_reads.ms16.dump`
  totKmers=`/Users/jjpc/miniforge3/envs/mitnanex/bin/meryl statistics ./c_aaseri_ERR10466724_collected_reads.ms16 | grep present | awk '{ print $2 }'`


  ./meryl-make-ignore.pl $numKmers $totKmers < ./c_aaseri_ERR10466724_collected_reads.ms16.dump | gzip -1c > ./c_aaseri_ERR10466724_collected_reads.ms16.ignore.gz
fi


exit 0
