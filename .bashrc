# .bashrc

# added by myself
alias ls='ls -G --color=force'
alias ll='ls -l'
alias la='ls -la'
alias dk='docker-interactive registry.gsc.wustl.edu/genome/genome_perl_environment'
#alias perlvep='LSF_DOCKER_PRESERVE_ENVIRONMENT=true bsub -Is -R"select[mem>15000] rusage[mem=15000]" -M 32000000 -q docker-interactive -a "docker(quay.io/biocontainers/perl-bioperl-core:1.007002--pl526_2)" /bin/bash'
#alias dockerg='sh ~/scripts/docker'
alias dockerhpc='LSF_DOCKER_PRESERVE_ENVIRONMENT=true bsub -Is -R"select[mem>15000] rusage[mem=15000]" -M 32000000 -q research-hpc -a "docker(registry.gsc.wustl.edu/genome/genome_perl_environment)" /bin/bash'
alias home='cd /gscmnt/gc2518/dinglab/scao/home'
alias mm='cd /gscmnt/gc2737/ding/scao/mmy'
alias dockeri='LSF_DOCKER_PRESERVE_ENVIRONMENT=true bsub -Is -R"select[mem>15000] rusage[mem=15000]" -M 32000000 -q docker-interactive -a "docker(registry.gsc.wustl.edu/genome/genome_perl_environment)" /bin/bash'
alias gscp='bsub -Is -q docker-interactive -a "docker(lbwang/dailybox)" /bin/bash --noprofile -l'
alias cao='LSF_DOCKER_PRESERVE_ENVIRONMENT=true bsub -Is -R"select[mem>15000] rusage[mem=15000]" -M 32000000 -q docker-interactive -a "docker(scao/dailybox)" /bin/bash'
alias caohpc='bsub -Is -q research-hpc -a "docker(scao/dailybox)" /bin/bash --noprofile -l'
alias dockerchris='bsub -Is -q docker-interactive -a "docker(chrisamiller/docker-genomic-analysis) /bin/bash --noprofile -l' 
export PS1="[\u@\h \W]\$ "
case $TERM in
    screen*)
        # This is the escape sequence ESC k \w ESC
        # Use current dir as the title
        SCREENTITLE='\[\ek\W\e\\\]'
        PS1="${SCREENTITLE}${PS1}"
        ;;
    *)
        ;;
esac

#export PATH=/gscuser/scao/bin:$PATH
#export PATH=/gscuser/scao/tools/cmake-2.8.10.2-Linux-i386/bin:$PATH
#export LD_LIBRARY_PATH=/gsc/lib:$LD_LIBRARY_PATH
#export mfold_datdir=/gscuser/scao/share/mfold:$mfold_datdir
#export DATAPATH=/gscuser/scao/tools/RNAstructure/data_tables/:$DATAPATH
#export PATH=/gscuser/scao/gc2524/dinglab/hotspot3d/hotspot3d/bin:$PATH
export PATH=/gscmnt/gc2524/dinglab/scao/conda_root/envs/vep/bin/:$PATH
export JAVA_HOME=/gscmnt/gc2525/dinglab/rmashl/Software/bin/jre/1.8.0_60-x64
export JAVA_OPTS="-Xms256m -Xmx512m"
export PATH=${JAVA_HOME}/bin:${PATH}
if [[ -z "$LD_LIBRARY_PATH" ]] ; then
export LD_LIBRARY_PATH=${JAVA_HOME}/lib
else
export LD_LIBRARY_PATH=${JAVA_HOME}/lib:${LD_LIBRARY_PATH}
fi

#PATH=/usr/bin:$PATH
#export PERL5LIB=/usr/lib/perl5/:$PERL5LIB

#export PATH=/opt/conda/bin

#PATH=/gsc/bin:$PATH
#somaticwrapper
#export PERL_PATH=/gsc
#export PATH=$PERL_PATH/bin:$PATH
#export PERL_BIN=$PERL_PATH/bin/perl
#export PERL5LIB=$PERL_PATH/lib/perl5/5.8.7/:$PERL5LIB


## gmt or somaticwrapper ###

#export PERL_PATH=/gscmnt/gc2525/dinglab/rmashl/Software/perl/perl-5.22.0
#export PATH=$PERL_PATH/bin:$PATH
#export PERL_BIN=$PERL_PATH/bin/perl
#export PERL5LIB=$PERL_PATH/lib/perl5:$PERL_PATH/lib/perl5/x86_64-linux:$PERL5LIB

### vep102
export PERL_PATH=/gscmnt/gc2524/dinglab/scao/conda_root/envs/vep/
export PATH=$PERL_PATH/bin:$PATH
export PERL_BIN=$PERL_PATH/bin/perl
export PERL5LIB=$PERL_PATH/lib/site_perl/5.26.2:$PERL_PATH/site_perl/5.26.2/x86_64-linux-thread-multi:$PERL5LIB

#/gscmnt/gc2524/dinglab/scao/conda_root/envs/vep/lib/

##repeatmasker ##
#export PERL_PATH=/gscmnt/gc2518/dinglab/scao/home/tools/anaconda2
#export PATH=$PERL_PATH/bin:$PATH
#export PERL_BIN=$PERL_PATH/bin/perl
#export PERL5LIB=$PERL_PATH/lib/5.26.2

#music2
#PATH=$PATH:/gscuser/mbailey/perl5/bin
#export PERL5LIB=/gscuser/mbailey/perl5/lib/perl5:$PERL5LIB

#STAR
PATH=/gscmnt/gc2741/ding/qgao/tools/STAR/source:$PATH

#STAR-Fusion
PATH=/gscmnt/gc2741/ding/qgao/tools/STAR-Fusion:$PATH

#export PATH=$PATH:~/perl5/bin/
#export PERL5LIB=PERL5LIB:~/perl5/lib/perl5/:$PERL5LIB

#export PERL5LIB=/gscuser/rjayasin/perl5/lib/perl5:$PERL5LIB
#export PERL5LIB=/gscuser/rjayasin/lib/:$PERL5LIB

#source /opt/lsf9/conf/profile.lsf 
# added by Anaconda2 5.3.1 installer
# >>> conda init >>>
# !! Contents within this block are managed by 'conda init' !!
### optitype, perl environment will change ##
#__conda_setup="$(CONDA_REPORT_ERRORS=false '/gscmnt/gc2518/dinglab/scao/home/tools/anaconda2/bin/conda' shell.bash hook 2> /dev/null)"
#if [ $? -eq 0 ]; then
#    \eval "$__conda_setup"
#else
#    if [ -f "/gscmnt/gc2518/dinglab/scao/home/tools/anaconda2/etc/profile.d/conda.sh" ]; then
# . "/gscmnt/gc2518/dinglab/scao/home/tools/anaconda2/etc/profile.d/conda.sh"  # commented out by conda initialize
#        CONDA_CHANGEPS1=false conda activate base
#    else
#       \export PATH="/gscmnt/gc2518/dinglab/scao/home/tools/anaconda2/bin:$PATH"
#    fi
#fi
#unset __conda_setup

# <<< conda init <<<

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/gscmnt/gc2518/dinglab/scao/home/tools/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/gscmnt/gc2518/dinglab/scao/home/tools/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/gscmnt/gc2518/dinglab/scao/home/tools/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/gscmnt/gc2518/dinglab/scao/home/tools/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

