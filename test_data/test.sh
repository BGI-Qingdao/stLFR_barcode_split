#!/bin/bash

###
# remove output files
###
rm -rf diff_cache split_reads.1.fq.gz split_reads.2.fq.gz
rm -rf split_reads.1.fq split_reads.2.fq
rm -rf barcode_freq.txt split_stat_read.log

###
# run command
###
../split_barcode.sh stLFR_raw.1.fq.gz stLFR_raw.2.fq.gz
###
# diff output files
###
gzip -dc split_reads.1.fq.gz >split_reads.1.fq
gzip -dc split_reads.2.fq.gz >split_reads.2.fq
diff expected_split_reads.1.fq split_reads.1.fq >diff_cache
diff expected_split_reads.2.fq split_reads.2.fq >>diff_cache
diff barcode_freq.txt expected_barcode_freq.txt >>diff_cache
diff split_read_stat.log expected_split_stat_read.log >>diff_cache
###
# check difffiles
###
diff_info=`wc -l diff_cache`

if [[ "$diff_info" == "0 diff_cache" ]] ; then 
    echo "test pass !"
else 
    echo "\"$diff_info\""
    echo "test failed ! see diff_cache for more details." 
fi
