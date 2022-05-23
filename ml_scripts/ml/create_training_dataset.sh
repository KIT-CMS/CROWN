#!/bin/bash
set -e

# source utils/setup_cvmfs_sft.sh
source utils/setup_python.sh
source utils/bashFunctionCollection.sh

echo "$@"
ERA=$1
CHANNEL=$2
MASS=$3
BATCH=$4
TAG_NAME=$5

tauEstimation=emb
jetEstimation=ff


TRAIN_STAGE_ARG="--nmssm"

# Create directories if needed
if [[ ! -d output/log/logandrun ]]; then
  echo "create output/log/logandrun"
  mkdir -p output/log/logandrun
fi

source utils/setup_samples.sh $ERA 

outdir=output/ml/${TAG_NAME}/${ERA}_${CHANNEL}_${MASS}_${BATCH}

[[ -d $outdir ]] ||  mkdir -p $outdir
echo $FF_Friends
if [[ ${CHANNEL} != 'em' ]]
then
  FRIENDS="${SVFit_Friends} ${HHKinFit_Friends} ${FF_Friends}"
else
  FRIENDS="${SVFit_Friends} ${HHKinFit_Friends}"
fi


# Write dataset config
logandrun python ml/write_dataset_config.py \
  --era ${ERA} \
  --channel ${CHANNEL} \
  --base-path $ARTUS_OUTPUTS \
  --friend-paths $FRIENDS \
  --database $KAPPA_DATABASE \
  --output-path $outdir \
  --mass $MASS \
  --batch $BATCH \
  --output-filename training_dataset.root \
  --tree-path ${CHANNEL}_nominal/ntuple \
  --event-branch event \
  --training-weight-branch training_weight \
  --training-z-estimation-method $tauEstimation \
  --training-jetfakes-estimation-method $jetEstimation \
  --output-config $outdir/dataset_config.yaml \
  $TRAIN_STAGE_ARG

# Create dataset files from config
./htt-ml/dataset/create_training_dataset.py $outdir/dataset_config.yaml

#  Reweight STXS stage 1 signals so that each stage 1 signal is weighted equally but
#  conserve the overall weight of the stage 0 signal
#  python ml/reweight_stxs_stage1.py \
#     $outdir \
#     $outdir/fold0_training_dataset.root \
#     $outdir/fold1_training_dataset.root

# reunite the dataset
# logandrun hadd -f $outdir/combined_training_dataset.root \
#     $outdir/fold0_training_dataset.root \
#     $outdir/fold1_training_dataset.root

python ./ml/sum_training_weights.py \
  --era ${ERA} \
  --channel ${CHANNEL} \
  --dataset $outdir/fold0_training_dataset.root \
            $outdir/fold1_training_dataset.root \
  --dataset-config-file "$outdir/dataset_config.yaml" \
  --training-template "ml/templates/${ERA}_${CHANNEL}_training.yaml" \
  --channel $CHANNEL \
  --write-weights True

# Remove intermediate root files
rm $outdir/merge_fold*.root -f