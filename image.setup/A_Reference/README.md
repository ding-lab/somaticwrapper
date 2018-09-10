Scripts here run from inside docker container to download and index reference files.

* `1_get_GRCh37.sh` processes GRCh37-lite, which has chromosomes labeled '1'
* `2_get_GRCh37.p13.sh` processes GRCh37.p13, which has chromosomes labeled 'chr1'.  See notes in that script for details
* `3_prepare_strelka_demo.sh` processes the `demo20.fa` reference used for the strelka demo work.  This reference needs to be generated 
outside of the host image with the script `somaticwrapper/docker/2_setup_data.sh`.

Processes only the references necessary for your work.  Note that steps 1 and 2 are time and storage intensive.

TODO: add scripts to download and process GRCh38
TODO: get Mouse: ftp://ftp.ensembl.org/pub/release-93/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz

