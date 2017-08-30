
sub bsub_varscan{
	my $sample_name = shift;
    my $sample_full_path = shift;
	my $job_files_dir = shift;
	my $bsub = shift;
	my $h37_REF = shift;

    $current_job_file = "j2_varscan_".$sample_name.".sh";
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
    print OUT "myRUNDIR=".$sample_full_path."/varscan\n";
    print OUT "STATUSDIR=".$sample_full_path."/status\n";
    print OUT "RESULTSDIR=".$sample_full_path."/varscan_results\n";
    print OUT "RUNDIR=".$sample_full_path."\n";

    print OUT "CONFDIR=/usr/local/somaticwrapper/config\n";  # is this necessary?
    print OUT "GENOMEVIP_SCRIPTS=/gscmnt/gc2525/dinglab/rmashl/Software/bin/genomevip\n";
    print OUT "export VARSCAN_DIR=/gscmnt/gc2525/dinglab/rmashl/Software/bin/varscan/2.3.8\n";
    print OUT "export SAMTOOLS_DIR=/gscmnt/gc2525/dinglab/rmashl/Software/bin/samtools/1.2/bin\n";
    print OUT "export JAVA_HOME=/gscmnt/gc2525/dinglab/rmashl/Software/bin/jre/1.8.0_60-x64\n";
    print OUT "export JAVA_OPTS=\"-Xms256m -Xmx512m\"\n";
    print OUT "export PATH=\${JAVA_HOME}/bin:\${PATH}\n";
    print OUT "if [ ! -d \${myRUNDIR} ]\n";
    print OUT "then\n";
    print OUT "mkdir \${myRUNDIR}\n";
    print OUT "fi\n";
    print OUT "if \[\[ -z \"\$LD_LIBRARY_PATH\" \]\] \; then\n";
    print OUT "export LD_LIBRARY_PATH=\${JAVA_HOME}/lib\n";
    print OUT "else\n";
    print OUT "export LD_LIBRARY_PATH=\${JAVA_HOME}/lib:\${LD_LIBRARY_PATH}\n";
    print OUT "fi\n";
    print OUT "put_cmd=\"ln -s\"\n";
    print OUT "del_cmd=\"rm -f\"\n";
    print OUT "del_local=\"rm -f\"\n";
    print OUT "statfile=incomplete.vs_som_snvindels\n";
    print OUT "localstatus=\${RUNDIR}\/status\/\${statfile}\n";
    print OUT "if [ ! -d \${myRUNDIR}\/status ]\n";
    print OUT "then\n";
    print OUT "mkdir \${myRUNDIR}\/status\n";
    print OUT "fi\n";
    print OUT "touch \${localstatus}\n";
    print OUT "cd \${RUNDIR}/varscan\n";
    print OUT "TMPBASE=.\/varscan.out.som\n";
    print OUT "LOG=\${TMPBASE}.log\n";
    print OUT "snvoutbase=\${TMPBASE}_snv\n";
    print OUT "indeloutbase=\${TMPBASE}_indel\n";
    print OUT "BAMLIST=\${RUNDIR}/varscan/bamfilelist.inp\n";
    print OUT "if [ ! -e \${BAMLIST} ]\n";
    print OUT "then\n";
    print OUT "rm \${BAMLIST}\n";
    print OUT "fi\n";
    print OUT "echo \"$IN_bam_N\" > \${BAMLIST}\n"; 
    print OUT "echo \"$IN_bam_T\" >> \${BAMLIST}\n";  
    print OUT "ncols=\$(echo \"3*( \$(wc -l < \$BAMLIST) +1)\"|bc)\n";
    print OUT "\${SAMTOOLS_DIR}/samtools mpileup -q 1 -Q 13 -B -f $h37_REF -b \${BAMLIST} | awk -v ncols=\$ncols \'NF==ncols\' | java \${JAVA_OPTS} -jar \${VARSCAN_DIR}/VarScan.jar somatic - \${TMPBASE} --mpileup 1 --p-value 0.99 --somatic-p-value 0.05 --min-coverage-normal 30 --min-coverage-tumor 22 --min-var-freq 0.08 --min-freq-for-hom 0.75 --normal-purity 1.00 --tumor-purity 1.00 --strand-filter 1 --min-avg-qual 15 --output-vcf 1 --output-snp \${snvoutbase} --output-indel \${indeloutbase} &> \${LOG}\n";
    close OUT;	
    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
	print($bsub_com."\n");

	print("Aborting.\n");
	die();
    system ( $bsub_com );

}

1;
