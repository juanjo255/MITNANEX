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
  blk="000001"
  slf=""
  qry="000001"
fi


if [ x$qry = x ]; then
  echo Error: Job index out of range.
  exit 1
fi

if [ -e ./results/$qry.ovb ]; then
  echo Job previously completed successfully.
  exit
fi


if [ -e ./queries.tar -a ! -d ./queries ] ; then
  tar -xf ./queries.tar
fi


if [ ! -d ./results ]; then  mkdir -p ./results;  fi
if [ ! -d ./blocks  ]; then  mkdir -p ./blocks;   fi

for ii in `ls ./queries/$qry` ; do
  echo Fetch blocks/$ii
done

echo ""
echo Running block $blk in query $qry
echo ""

if [ ! -e ./results/$qry.mhap.ovb ] ; then
  #  Make a fifo so we can check return status on both
  #  mhap and mhapConvert, and still pipe results so we
  #  stop filling up disks.
  rm -f  $qry-pipe
  mkfifo $qry-pipe

  #  Start up the consumer.
  $bin/mhapConvert \
    -S ../../c_aaseri_ERR10466724_collected_reads.seqStore \
    -o ./results/$qry.mhap.ovb.WORKING \
    -minlength 500 \
    $qry-pipe \
  && \
  touch ./results/$qry.mcvt.success &

  #  Start up the producer.
  /Users/jjpc/miniforge3/envs/mitnanex/lib/jvm/bin/java  -XX:ParallelGCThreads=8 -server -Xms5530m -Xmx5530m \
    -jar  $bin/../share/java/classes/mhap-2.1.3.jar  \
    --repeat-weight 0.9 --repeat-idf-scale 10 -k 16 \
    --store-full-id \
    --num-hashes 256 \
    --num-min-matches 3 \
    --threshold 0.8 \
    --filter-threshold 0.0000001 \
    --ordered-sketch-size 1000 \
    --ordered-kmer-size 14 \
    --min-olap-length 500 \
    --num-threads 8 \
    -s  ./blocks/$blk.dat $slf  \
    -q  queries/$qry  \
  > $qry-pipe \
  && \
  touch ./results/$qry.mhap.success

  #  Wait for mhapConvert to finish.
  wait

  #  Now that they're done, check status.
  if [ -e ./results/$qry.mhap.success -a -e ./results/$qry.mcvt.success ] ; then
    mv ./results/$qry.mhap.ovb.WORKING ./results/$qry.mhap.ovb
    rm -f ./results/$qry.mhap.success
    rm -f ./results/$qry.mcvt.success
  fi

  #  And destroy the pipe.
  rm -f  $qry-pipe
fi

if [ -e ./results/$qry.mhap.ovb ] ; then
  mv -f ./results/$qry.mhap.ovb ./results/$qry.ovb
fi


exit 0
