
CHANNEL=$1
MASS=$2
BATCH=$3
BASE_PATH=$4

TAG="${MASS}_${BATCH}"

outdir=${BASE_PATH}/all_eras_${CHANNEL}_${TAG}

# source utils/setup_cvmfs_sft.sh
python ml/create_combined_config.py  --tag ${TAG} --input_base_path ${BASE_PATH} --channel $CHANNEL --output_dir ${outdir}


# cat $outdir/dataset_config.yaml