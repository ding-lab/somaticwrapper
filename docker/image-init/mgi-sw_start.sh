# This script to be executed within new BPS container running in MGI context
# replicates behavior of running BPS in native docker context

/bin/bash --rcfile /home/sw/mgi-sw.bashrc
cd /usr/local/somaticwrapper  # TODO: this should be modified to correspond to SW_HOME.C
