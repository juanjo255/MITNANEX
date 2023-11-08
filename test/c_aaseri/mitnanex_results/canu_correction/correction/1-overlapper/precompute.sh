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

if [ $jobid -eq 1 ] ; then
  rge="-r 1-8231"
  job="000001"
fi


if [ x$job = x ] ; then
  echo Job partitioning error.  jobid $jobid is invalid.
  exit 1
fi

if [ ! -d ./blocks ]; then
  mkdir -p ./blocks
fi


if [ -e blocks/$job.dat ]; then
  exists=true
else
  exists=false
fi
if [ $exists = true ] ; then
  echo Job previously completed successfully.
  exit
fi

#  Grab the kmers.ignore.gz file.
mkdir -p ../0-mercounts/
cd       ../0-mercounts/
cd -

#  Extract the sequences we want to compute on.
$bin/sqStoreDumpFASTQ \
  -S ../../c_aaseri_ERR10466724_collected_reads.seqStore \
  $rge \
  -nolibname \
  -noreadname \
  -fasta \
  -o ./blocks/$job.input \
|| \
mv -f ./blocks/$job.input.fasta ./blocks/$job.input.fasta.FAILED

if [ ! -e ./blocks/$job.input.fasta ] ; then
  echo Failed to extract fasta.
  exit 1
fi

echo ""
echo Starting mhap precompute.
echo ""

#  Remove any previous result.
rm -f ./blocks/$job.input.dat

#  So mhap writes its output in the correct spot.
cd ./blocks

/Users/jjpc/miniforge3/envs/mitnanex/lib/jvm/bin/java  -XX:ParallelGCThreads=8 -server -Xms5530m -Xmx5530m \
  -jar  $bin/../share/java/classes/mhap-2.1.3.jar  \
  --repeat-weight 0.9 --repeat-idf-scale 10 -k 16 \
  --store-full-id \
  --num-hashes 256 \
  --num-min-matches 3 \
  --ordered-sketch-size 1000 \
  --ordered-kmer-size 14 \
  --threshold 0.8 \
  --filter-threshold 0.0000001 \
  --min-olap-length 500 \
  --num-threads 8 \
  -f  ../../0-mercounts/c_aaseri_ERR10466724_collected_reads.ms16.ignore.gz  \
  -p  ./$job.input.fasta  \
  -q  .  \
&& \
mv -f ./$job.input.dat ./$job.dat

if [ ! -e ./$job.dat ] ; then
  echo Mhap failed.
  exit 1
fi

#  Clean up, remove the fasta input
rm -f ./$job.input.fasta


exit 0
