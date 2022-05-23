#!/bin/bash
# set -e
trap 'kill $(jobs -p); exit 1' SIGINT
ERA=$1
CHANNEL=$2
TAG=$3
TAG_2=$4
CUT_FRAC=$5
NODENUM=$6
LAYERNUM=$7
RM_FRAC=$8
# export altpara=$8

echo "$0 $@"

export NODE_NUM=${NODENUM}
export LAYER_NUM=${LAYERNUM}

outdir=output/ml/${TAG_2}/${ERA}_${CHANNEL}_${TAG}
if [[ ! -d ${outdir} ]]; then
  mkdir ${outdir}
fi 

if [[ "${ABS_NORM_WEIGHTS}" == "True" ]]; then
  ABS_NORM_WEIGHTS_STRING="--abw 1"
fi

# PERC_NE=0.1

if [[ ! -z ${PERC_NE} ]]; then
  PERC_NE_STRING="--addn ${PERC_NE}"
fi

if [[ ! -z ${CUT_FRAC} ]]; then
  EXTEND_BATCH="--extbatch ${CUT_FRAC}"
fi

if [[ ! -z ${RM_FRAC} ]]; then
  REMOVE_BATCH="--sample_frac ${RM_FRAC}"
fi

# source utils/setup_cvmfs_sft_individual.sh requirements_training.txt 101cuda
# cat requirements_training.txt
if [[ ${_CONDOR_JOB_IWD} == "/var/lib/condor/execute/dir"* ]]; then
    echo "Running on HTCondor"
    echo "Install start"
    date +"%T"
    python3 -m pip install --target=$_CONDOR_JOB_IWD/packages -U -r requirements_training.txt #> ${outdir}/install_log.txt
    echo "Install end"
    date +"%T"
    export PYTHONPATH=$_CONDOR_JOB_IWD/packages:$PYTHONPATH;
else
    echo "Running locally"
    python3 -m pip install -U -r requirements_training.txt #> ${outdir}/install_log.txt
fi
source utils/bashFunctionCollection.sh

python3 -m pip freeze

# if uname -a | grep ekpdeepthought; then
#     source ~mscham/p/bms3/root-dev/build2/bin/thisroot.sh
#     if [[ -z USE_CPU ]]; then
#       source utils/setup_cuda.sh
#       export PYTHONUSERBASE=$HOME/.local/pylibs-$(hostname)
#       export PYTHONPATH=$HOME/.local/pylibs-$(hostname)/lib/python2.7/site-packages:$PYTHONPATH
#       export PATH=$HOME/.local/pylibs-$(hostname)/bin:$PATH
#     fi
#     source utils/setup_cvmfs_sft.sh
# else
#     source utils/setup_cvmfs_sft.sh
#     export PYTHONUSERBASE=$HOME/.local/pylibs-$(hostname)
#     export PYTHONPATH=$HOME/.local/pylibs-$(hostname)/lib/python2.7/site-packages:$PYTHONPATH
#     export PATH=$HOME/.local/pylibs-$(hostname)/bin:$PATH
# fi

# set number of CPUs to 12 if no other value has been set 
# export OMP_NUM_THREADS=1
if [[ -z ${OMP_NUM_THREADS} ]]; then
  export OMP_NUM_THREADS=12
fi
export TF_GPU_THREAD_MODE="gpu_private"
# export CUDA_VISIBLE_DEVICES=""

CPU_START=20
CPU_END=$(( ${CPU_START}+${OMP_NUM_THREADS}-1 ))

# python3 monitor.py --live-update --columns "name,create_time,cores,cpu_usage,status,nice,memory_usage,read_bytes,write_bytes,n_threads,username" > ${outdir}/monitor_logs_training.txt &
# BACK_PID=$!
# echo "Start monitoring with pid ${BACK_PID}."


# rm ${outdir}/monitor_logs_training*.txt
mkdir -p $outdir
if [[ $ERA == *"all"* ]]
then
  python3 htt-ml/training/keras_training.py ${outdir}/dataset_config.yaml \
    0 --balance-batches 1 --conditional 1 ${ABS_NORM_WEIGHTS_STRING} ${PERC_NE_STRING} ${EXTEND_BATCH} ${REMOVE_BATCH} #& #--randomization 1
  # python3 htt-ml/training/keras_training.py ${outdir}/dataset_config.yaml \
  #   0 --conditional 1 ${ABS_NORM_WEIGHTS_STRING} ${PERC_NE_STRING} ${EXTEND_BATCH} ${REMOVE_BATCH} #& #--randomization 1
  # taskset -c ${CPU_START}-${CPU_END} python3 htt-ml/training/keras_training.py ${outdir}/dataset_config.yaml \
  #   0 --balance-batches 1 --conditional 1 ${ABS_NORM_WEIGHTS_STRING} ${PERC_NE_STRING} ${EXTEND_BATCH} ${REMOVE_BATCH} & #--randomization 1
  # TRAIN_PID=$!
  # echo "This is the training: ${TRAIN_PID}"
  # (while [[ True ]]; do ps -aux | grep " ${TRAIN_PID} " | grep training >> ${outdir}/monitor_logs_training0.txt; sleep 10; done;) &
  # MONITOR_PID=$!
  # echo "Start monitoring with pid ${MONITOR_PID}."
  # wait ${TRAIN_PID}
  # kill -9 ${MONITOR_PID}

  python3 htt-ml/training/keras_training.py ${outdir}/dataset_config.yaml \
    1 --balance-batches 1 --conditional 1 ${ABS_NORM_WEIGHTS_STRING} ${PERC_NE_STRING} ${EXTEND_BATCH} ${REMOVE_BATCH} #& #--randomization 1
  # python3 htt-ml/training/keras_training.py ${outdir}/dataset_config.yaml \
  #   1 --conditional 1 ${ABS_NORM_WEIGHTS_STRING} ${PERC_NE_STRING} ${EXTEND_BATCH} ${REMOVE_BATCH} #& #--randomization 1
  # taskset -c ${CPU_START}-${CPU_END} python3 htt-ml/training/keras_training.py ${outdir}/dataset_config.yaml \
  #   1 --balance-batches 1 --conditional 1 ${ABS_NORM_WEIGHTS_STRING} ${PERC_NE_STRING} ${EXTEND_BATCH} ${REMOVE_BATCH} & #--randomization 1
  # TRAIN_PID=$!
  # echo "This is the training: ${TRAIN_PID}"
  # (while [[ True ]]; do ps -aux | grep " ${TRAIN_PID} " | grep training >> ${outdir}/monitor_logs_training1.txt; sleep 10; done;) &
  # MONITOR_PID=$!
  # echo "Start monitoring with pid ${MONITOR_PID}."
  # wait ${TRAIN_PID}
  # kill -9 ${MONITOR_PID}
else
  python3 htt-ml/training/keras_training.py ${outdir}/dataset_config.yaml 0 \
  --balance-batches 1 ${ABS_NORM_WEIGHTS_STRING} ${PERC_NE_STRING} ${EXTEND_BATCH} ${REMOVE_BATCH} &
  TRAIN_PID=$!
  echo "This is the training: ${TRAIN_PID}"
  (while [[ True ]]; do ps -aux | grep " ${TRAIN_PID} " | grep training >> ${outdir}/monitor_logs_training0.txt; sleep 10; done;) &
  MONITOR_PID=$!
  echo "Start monitoring with pid ${MONITOR_PID}."
  wait ${TRAIN_PID}
  kill -9 ${MONITOR_PID}
  python3 htt-ml/training/keras_training.py ${outdir}/dataset_config.yaml 1 \
  --balance-batches 1 ${ABS_NORM_WEIGHTS_STRING} ${PERC_NE_STRING} ${EXTEND_BATCH} ${REMOVE_BATCH} &
  TRAIN_PID=$!
  echo "This is the training: ${TRAIN_PID}"
  (while [[ True ]]; do ps -aux | grep " ${TRAIN_PID} " | grep training >> ${outdir}/monitor_logs_training1.txt; sleep 10; done;) &
  MONITOR_PID=$!
  echo "Start monitoring with pid ${MONITOR_PID}."
  wait ${TRAIN_PID}
  kill -9 ${MONITOR_PID}
fi
# kill -9 ${BACK_PID}
# echo "End monitoring."