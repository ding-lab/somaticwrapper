# Test case
# running just vep.  This is unnecessary and we exit early

echo Output of run_vep\(\) step is ignored in the pipeline.  Skipping step.
exit

# However, if you really want to run_vep(), just comment out the stuff above and run this:
CONFIG=/data/data/SWtest/sw.config

perl SomaticWrapper.pl /data/data 6 $CONFIG

