#!/bin/bash 
echo
echo '****************************************'
echo '*** NeuroMiner Elessar DEV           ***'
echo '*** SGE joblist manager:             ***'
echo '*** Preprocess features              ***'
echo '*** (c) 2017 N. Koutsouleris         ***'
echo '*** Updated Sep2020: Dom             ***'
echo '****************************************'
echo '     This is the STABLE      version.   '
echo '         Expect bugs and report         '
echo '****************************************'
echo   

# compiled with matlab R2018b so MCR main is v95. Needs to change if different MCR is used.
export LD_LIBRARY_PATH=/opt/matlab/v95/runtime/glnxa64:/opt/matlab/v95/bin/glnxa64:/opt/matlab/v95/sys/os/glnxa64:/opt/matlab/v95/sys/opengl/lib/glnxa64
export JOB_DIR=$PWD
export NEUROMINER=/opt/NM/NeuroMinerMCCMain_1.05_v95/for_testing 
export ACTION=preproc
read -e -p 'Path to NM structure: ' datpath
if [ ! -f $datpath ]; then
 	echo $datpath' not found.'
 	exit 
fi

read -e -p 'Path to job directory ['$JOB_DIR']: ' tJOB_DIR
if [ "$tJOB_DIR" != '' ]; then
	if [ -d $tJOB_DIR ]; then 
		export JOB_DIR=$tJOB_DIR
	else
		echo $tJOB_DIR' not found.'
		exit
	fi
fi 

read -e -p 'Change path to compiled NM directory ['$NEUROMINER']: ' tNEUROMINER
if [ "$tNEUROMINER" != '' ]; then     
  if [ -d $tNEUROMINER ]; then  
    export NEUROMINER=$tNEUROMINER
  else
    echo $tNEUROMINER' not found.'
    exit
  fi    
fi
read -p 'Provide your email address: ' EMAIL
echo '-----------------------'
echo 'PATH definitions:'
echo 'LOG directory: '$JOB_DIR
echo 'NeuroMiner directory: '$NEUROMINER
echo '-----------------------'
read -p 'Index to analysis container (NM.analysis{<index>}): ' analind
if [ "$analind" = '' ] ; then
	echo 'An analysis index is mandatory! Exiting program.'
	exit   
fi
read -p 'CV2 grid start row: ' CV2x1
read -p 'CV2 grid end row: ' CV2x2
read -p 'CV2 grid start column: ' CV2y1
read -p 'CV2 grid end column: ' CV2y2
read -p 'No. of SGE jobs: ' numCPU
read -p 'Server to use [any=1, psy0cf20=2, mitnvp1=3]: ' sind
if [ "$sind" = '1' ]; then
        SERVER_ID='all.q'
        echo "WARNING: if it is a high RAM job then please use psy0cf20"
elif [ "$sind" = '2' ]; then
        SERVER_ID='psy0cf20'
        echo "Please estimate RAM accurately"
elif [ "$sind" = '3' ]; then
        SERVER_ID='mitnvp1-2'
        echo "WARNING: if it is a high RAM job then please use psy0cf20"
else
        echo "Enter a number between 1-3"
fi
read -p 'xxxx MB RAM / SGE job: ' MB
MB=$MB'M'
read -p 'Use OpenMP [yes = 1 | no = 0]: ' pLibsvm
if [ "$pLibsvm" = '1' ] ; then
	read -p 'Specify number of CPUs assigned to MATLAB job [4, 8, 16, 32]: ' pnum
	PMODE='#\$-pe omp '$pnum
else
	PMODE=''
	pnum=1
fi
read -p 'Overwrite existing PreprocData files [yes = 1 | no = 2] ' ovrwrt
read -p 'Submit jobs immediately [y/n]: ' todo
for curCPU in $(seq $((numCPU)))
do
SD='_CPU'$curCPU
ParamFile=$JOB_DIR/Param_NM_$ACTION$SD
SGEFile=$JOB_DIR/NM_$ACTION$SD
echo 'Generate parameter file: NM_'$ACTION$SD' => '$ParamFile
# Generate parameter file
cat > $ParamFile <<EOF
$NEUROMINER
$datpath
$analind
$curCPU
$numCPU
$CV2x1
$CV2x2
$CV2y1
$CV2y2
$ovrwrt
EOF
cat > $SGEFile <<EOF
#!/bin/bash
#\$-o $JOB_DIR/\$JOB_IDnm$ACTION$SD -j y
#\$-N nm$ACTION$SD
#\$-S /bin/bash
#\$-M $EMAIL
#\$-m ae
#\$-l mem_total=$MB
#\$-q $SERVER_ID
$PMODE
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH
export OMP_NUM_THREADS=$pnum
cd $NEUROMINER 
./NeuroMinerMCCMain $ACTION $ParamFile
EOF
chmod u+x $SGEFile
datum=`date +"%Y%m%d"`
if [ "$todo" = 'y' -o "$todo" = 'Y' ] ; then
qsub $SGEFile >> NeuroMiner_Preproc_$datum.log
fi
done
