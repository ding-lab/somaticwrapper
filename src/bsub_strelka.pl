
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
	my $outfn = "$job_files_dir/$current_job_file";
	print("Writing to $outfn\n");
    open(OUT, ">$outfn") or die $!;

    print OUT "#!/bin/bash\n";
    print OUT "scr_t0=\`date \+\%s\`\n";
    print OUT "TBAM=".$sample_full_path."/".$sample_name.".T.bam\n";
    print OUT "NBAM=".$sample_full_path."/".$sample_name.".N.bam\n";
    print OUT "myRUNDIR=".$sample_full_path."/strelka\n";
    print OUT "STATUSDIR=".$sample_full_path."/status\n";
    print OUT "mkdir -p \$STATUSDIR\n";
    print OUT "RESULTSDIR=".$sample_full_path."/results\n";
    print OUT "SG_DIR=".$sample_full_path."/strelka\n"; 
    print OUT "RUNDIR=".$sample_full_path."\n";
    print OUT "STRELKA_OUT=".$sample_full_path."/strelka/strelka_out"."\n";
    print OUT "CONFDIR=/usr/local/somaticwrapper/config\n";
    print OUT "export SAMTOOLS_DIR=/usr/local/bin\n";
    print OUT "export JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64/jre\n";
    print OUT "export JAVA_OPTS=\"-Xms256m -Xmx512m\"\n";
    print OUT "export PATH=\${JAVA_HOME}/bin:\${PATH}\n";
    print OUT "if [ ! -d \${myRUNDIR} ]\n";
    print OUT "then\n";
    print OUT "mkdir \${myRUNDIR}\n";
    print OUT "fi\n";
    print OUT "if [ -d \${STRELKA_OUT} ]\n";
    print OUT "then\n";
    print OUT "rm -rf \${STRELKA_OUT}\n";
    print OUT "fi\n";
    print OUT "if \[ -z \"\$LD_LIBRARY_PATH\" \] \; then\n"; 
    print OUT "export LD_LIBRARY_PATH=\${JAVA_HOME}/lib\n";
    print OUT "else\n";
    print OUT "export LD_LIBRARY_PATH=\${JAVA_HOME}/lib:\${LD_LIBRARY_PATH}\n";
    print OUT "fi\n";
    print OUT "put_cmd=\"ln -s\"\n";
    print OUT "del_cmd=\"rm -f\"\n";
    print OUT "del_local=\"rm -f\"\n";
    print OUT "statfile=incomplete.strelka\n";
    print OUT "localstatus=".$sample_full_path."/status/\$statfile\n";
    print OUT "touch \$localstatus\n";
    print OUT "   ".$STRELKA_DIR."/bin/configureStrelkaWorkflow.pl --normal \$NBAM --tumor \$TBAM --ref ". $h37_REF." --config \$CONFDIR/strelka.ini --output-dir \$STRELKA_OUT\n";
    print OUT "cd \$STRELKA_OUT\n";
    print OUT "make -j 16\n";
    close OUT;
    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";

	print($bsub_com."\n");
    system ( $bsub_com );
}

1;
