#!/bin/bash 
  
echo
echo '****************************************'
echo '*** NeuroMiner Elessar               ***'
echo '*** SGE joblist manager:             ***'
echo '*** Train and crossvalidate models   ***'
echo '*** (c) 2017 N. Koutsouleris         ***'
echo '*** Updated Sep20 Dom                ***'
echo '****************************************'
echo '     This is the STABLE  NM  version.   '
echo '         Expect bugs and report         '
echo ' This .sh script is experimental.       ' 
echo ' it is designed to run multiple NM analyses simultaneously '
echo '****************************************'
echo   

# compiled with matlab R2018b so MCR main is v95. Needs to change if different MCR is used.
export LD_LIBRARY_PATH=/opt/matlab/v95/runtime/glnxa64:/opt/matlab/v95/bin/glnxa64:/opt/matlab/v95/sys/os/glnxa64:/opt/matlab/v95/sys/opengl/lib/glnxa64

export JOB_DIR=$PWD
export NEUROMINER=/opt/NM/NeuroMinerMCCMain_1.05_v95/for_testing 
export ACTION=train
read -e -p 'Path to NM structure: ' datpath
if [ ! -f $datpath ]; then
 	echo $datpath' not found.'
 	exit 
fi

export masterpath=''  
read -p 'Define processing mode [ 1 = from scratch, 2 = using PreprocData Master file ] ' procmode
if [ "$procmode" = '2' ] ; then	
	read -e -p 'Path to PreprocData Master file: ' masterpath
	if [ ! -f $masterpath ]; then
		echo $masterpath' not found.' exit
	fi
else
	masterpath=NaN
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

read -e -p 'Change path to NeuroMiner directory ['$NEUROMINER']: ' tNEUROMINER
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

read -p 'Provide indices to analysis containers (NM.analysis{<index>})' analind

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
read -p 'Overwrite existing CVdatamats [yes = 1 | no = 2]: ' ovrwrt
read -p 'Submit jobs immediately (y): ' todo

for curanalind in $analind
do
for curCPU in $(seq $((numCPU)))
do
CA='_ANALYSIS'$curanalind
SD='_CPU'$curCPU
ParamFile=$JOB_DIR/Param_NM$ACTION$SD$CA
SGEFile=$JOB_DIR/NM_$ACTION$SD$CA
echo 'Generate parameter file: svm'$SD$CA' => '$ParamFile
# Generate parameter file
cat > $ParamFile <<EOF
$NEUROMINER
$datpath
$masterpath
$curanalind
$ovrwrt
$curCPU
$numCPU
$CV2x1
$CV2x2
$CV2y1
$CV2y2
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
qsub $SGEFile >> NeuroMiner_Train_$datum.log
fi
done
done
