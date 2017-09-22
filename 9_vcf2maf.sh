# vcf2maf step

# NOTE: to define configuration file, set on command line,
#   export CONFIG="/gscuser/mwyczalk/projects/SomaticWrapper/config/01BR001.config"
# before running this script.
# This will overwrite the default configuration file below
if [ -z $CONFIG ]; then
    CONFIG=/data/data/SWtest/sw.config
fi
perl SomaticWrapper.pl /data/data 9 $CONFIG
