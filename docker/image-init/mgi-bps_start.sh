# This script to be executed within new BPS container running in MGI context
# replicates behavior of running BPS in native docker context

/bin/bash --rcfile /home/bps/mgi-bps.bashrc
cd /usr/local/somaticwrapper
