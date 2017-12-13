# run_varscan step

#echo Output of run_vep\(\) step is ignored in the pipeline.  Skipping step.
#exit

# However, if you really want to run_vep(), just comment out the stuff above and run this:

# NOTE: to define configuration file, set on command line,
#   export CONFIG="/gscuser/mwyczalk/projects/SomaticWrapper/config/01BR001.config"
# before running this script.
# This will overwrite the default configuration file below
if [ -z $CONFIG ]; then
    CONFIG=/data/data/SWtest/sw.config
fi

perl SomaticWrapper.pl 10 $CONFIG

