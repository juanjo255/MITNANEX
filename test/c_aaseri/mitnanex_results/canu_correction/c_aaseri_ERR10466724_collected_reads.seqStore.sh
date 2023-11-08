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





/Users/jjpc/miniforge3/envs/mitnanex/bin/sqStoreCreate \
  -o ./c_aaseri_ERR10466724_collected_reads.seqStore.BUILDING \
  -minlength 1000 \
  -genomesize 40000 \
  -coverage   200 \
  -bias       0 \
  -raw -nanopore c_aaseri_ERR10466724_collected_reads.filtlong /Users/jjpc/Documents/TESIS/MITNANEX_PROJECT/test/c_aaseri/mitnanex_results/c_aaseri_ERR10466724_collected_reads.filtlong.fastq \
&& \
mv ./c_aaseri_ERR10466724_collected_reads.seqStore.BUILDING ./c_aaseri_ERR10466724_collected_reads.seqStore \
&& \
exit 0

exit 1
