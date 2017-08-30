
# Usage: bsub_strelka($sample_name )
sub bsub_strelka {
	my $sample_name = shift;
    my $sample_full_path = shift;
	my $job_files_dir = shift;
	my $bsub = shift;
    my $STRELKA_DIR = shift;
	my $h37_REF = shift;

    $current_job_file = "j1_streka_".$sample_name.".sh"; 
    my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
    my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";
    if (! -e $IN_bam_T) {#make sure there is a input fasta file 
        print $red,  "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
        print "Warning: Died because there is no input bam file for bwa:\n";
        print "File $IN_bam_T does not exist!\n";
        die "Please check command line argument!", $normal, "\n\n";

    }
    if (! -s $IN_bam_T) {#make sure input fasta file is not empty
        print $red, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
        die "Warning: Died because $IN_bam_T is empty!", $normal, "\n\n";
    }
    if (! -e $IN_bam_N) {#make sure there is a input fasta file 
        print $red,  "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
        print "Warning: Died because there is no input bam file for bwa:\n";
        print "File $IN_bam_N does not exist!\n";
        die "Please check command line argument!", $normal, "\n\n";

    }
    if (! -s $IN_bam_N) {#make sure input fasta file is not empty
        print $red, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
        die "Warning: Died because $IN_bam_N is empty!", $normal, "\n\n";
    }
	my $strelka_out = "$job_files_dir/$current_job_file";
	print("Writing to $strelka_out\n");
    open(STREKA, ">$strelka_out") or die $!;
    print STREKA "#!/bin/bash\n";
    print STREKA "#BSUB -n 1\n";
    print STREKA "#BSUB -R \"rusage[mem=30000]\"","\n";
    print STREKA "#BSUB -M 30000000\n";
    print STREKA "#BSUB -J $current_job_file\n";
    print STREKA "scr_t0=\`date \+\%s\`\n";
    print STREKA "TBAM=".$sample_full_path."/".$sample_name.".T.bam\n";
    print STREKA "NBAM=".$sample_full_path."/".$sample_name.".N.bam\n";
    print STREKA "myRUNDIR=".$sample_full_path."/strelka\n";
    print STREKA "STATUSDIR=".$sample_full_path."/status\n";
    print STREKA "mkdir -p \$STATUSDIR\n";
    print STREKA "RESULTSDIR=".$sample_full_path."/results\n";
    print STREKA "SG_DIR=".$sample_full_path."/strelka\n"; 
    print STREKA "RUNDIR=".$sample_full_path."\n";
    print STREKA "STRELKA_OUT=".$sample_full_path."/strelka/strelka_out"."\n";
    print STREKA "CONFDIR=/usr/local/somaticwrapper/config\n";
    print STREKA "export SAMTOOLS_DIR=/usr/local/bin\n";
    print STREKA "export JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64/jre\n";
    print STREKA "export JAVA_OPTS=\"-Xms256m -Xmx512m\"\n";
    print STREKA "export PATH=\${JAVA_HOME}/bin:\${PATH}\n";
    print STREKA "if [ ! -d \${myRUNDIR} ]\n";
    print STREKA "then\n";
    print STREKA "mkdir \${myRUNDIR}\n";
    print STREKA "fi\n";
    print STREKA "if [ -d \${STRELKA_OUT} ]\n";
    print STREKA "then\n";
    print STREKA "rm -rf \${STRELKA_OUT}\n";
    print STREKA "fi\n";
    print STREKA "if \[ -z \"\$LD_LIBRARY_PATH\" \] \; then\n"; 
    print STREKA "export LD_LIBRARY_PATH=\${JAVA_HOME}/lib\n";
    print STREKA "else\n";
    print STREKA "export LD_LIBRARY_PATH=\${JAVA_HOME}/lib:\${LD_LIBRARY_PATH}\n";
    print STREKA "fi\n";
    print STREKA "put_cmd=\"ln -s\"\n";
    print STREKA "del_cmd=\"rm -f\"\n";
    print STREKA "del_local=\"rm -f\"\n";
    print STREKA "statfile=incomplete.strelka\n";
    print STREKA "localstatus=".$sample_full_path."/status/\$statfile\n";
    print STREKA "touch \$localstatus\n";
    print STREKA "   ".$STRELKA_DIR."/bin/configureStrelkaWorkflow.pl --normal \$NBAM --tumor \$TBAM --ref ". $h37_REF." --config \$CONFDIR/strelka.ini --output-dir \$STRELKA_OUT\n";
    print STREKA "cd \$STRELKA_OUT\n";
    print STREKA "make -j 16\n";
    close STREKA;
    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
	print($bsub_com."\n");
die();
    system ( $bsub_com );
}

1;
