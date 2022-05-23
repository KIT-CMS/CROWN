#!/bin/bash
set -e

ERA=$1
TAG=$2

if [[ ! -z $3 ]]; then
    OUTPUTDIR=$3
elif [[ -d output/friend_trees ]];then
    OUTPUTDIR=$( cd output; pwd -P)
fi

# Kappa database
KAPPA_DATABASE=datasets/datasets.json

#### ERA specific part. If a sample is not available comment it out here.
# Samples Run2016
ARTUS_OUTPUTS_2016="/ceph/jbechtel/nmssm/ntuples/2016/+CH+/"
SVFit_Friends_2016="/ceph/jbechtel/nmssm/friends/2016/+CH+/SVFit/"
FF_Friends_2016="/ceph/jbechtel/nmssm/friends/2016/+CH+/FakeFactors_nmssm/"
HHKinFit_Friends_2016="/ceph/jbechtel/nmssm/friends/2016/+CH+/HHKinFit/"
#NNScore_Friends_2016="/ceph/jbechtel/nmssm/friends/2016/+CH+/NNScore_nmssm_v5/NNScore_workdir/NNScore_collected/"
NNScore_Friends_2016="/ceph/jbechtel/nmssm/friends/2016/+CH+/NNScore_train_all/NNScore_workdir/+MASS+_+BATCH+/NNScore_workdir/NNScore_collected/"

# Samples Run2017
ARTUS_OUTPUTS_2017="/ceph/jbechtel/nmssm/ntuples/2017/+CH+/"
SVFit_Friends_2017="/ceph/jbechtel/nmssm/friends/2017/+CH+/SVFit/"
FF_Friends_2017="/ceph/jbechtel/nmssm/friends/2017/+CH+/FakeFactors_nmssm/"
#FF_Friends_2017="/ceph/jbechtel/nmssm/friends/2017/+CH+/FakeFactors_medium_v5/FakeFactors_workdir/FakeFactors_collected/"
HHKinFit_Friends_2017="/ceph/jbechtel/nmssm/friends/2017/+CH+/HHKinFit/"
#NNScore_Friends_2017="/ceph/jbechtel/nmssm/friends/2017/+CH+/NNScore_nmssm_v5/NNScore_workdir/NNScore_collected/"
NNScore_Friends_2017="/ceph/jbechtel/nmssm/friends/2017/+CH+/NNScore_train_all/NNScore_workdir/+MASS+_+BATCH+/NNScore_workdir/NNScore_collected/"

# Samples Run2018
ARTUS_OUTPUTS_2018="/ceph/jbechtel/nmssm/ntuples/2018/+CH+/"
SVFit_Friends_2018="/ceph/jbechtel/nmssm/friends/2018/+CH+/SVFit/"
#FF_Friends_2018="/ceph/jbechtel/nmssm/friends/2018/+CH+/FakeFactors_medium/"
FF_Friends_2018="/ceph/jbechtel/nmssm/friends/2018/+CH+/FakeFactors_nmssm/"
#FF_Friends_2018="/ceph/jbechtel/nmssm/friends/2018/+CH+/FakeFactors_medium_v7/FakeFactors_workdir/FakeFactors_collected/"
HHKinFit_Friends_2018="/ceph/jbechtel/nmssm/friends/2018/+CH+/HHKinFit/"
#NNScore_Friends_2018="/ceph/jbechtel/nmssm/friends/2018/+CH+/NNScore_nmssm_v5/NNScore_workdir/NNScore_collected/"
NNScore_Friends_2018="/ceph/jbechtel/nmssm/friends/2018/+CH+/NNScore_train_all/NNScore_workdir/+MASS+_+BATCH+/NNScore_workdir/NNScore_collected/"


# ERA handling
if [[ $ERA == *"2016"* ]]
then
    ARTUS_OUTPUTS=$ARTUS_OUTPUTS_2016
    NNScore_Friends=$NNScore_Friends_2016
    SVFit_Friends=$SVFit_Friends_2016
    MELA_Friends=$MELA_Friends_2016
    FF_Friends=$FF_Friends_2016
    HHKinFit_Friends=$HHKinFit_Friends_2016
elif [[ $ERA == *"2017"* ]]
then
    ARTUS_OUTPUTS=$ARTUS_OUTPUTS_2017
    NNScore_Friends=$NNScore_Friends_2017
    SVFit_Friends=$SVFit_Friends_2017
    MELA_Friends=$MELA_Friends_2017
    FF_Friends=$FF_Friends_2017
    HHKinFit_Friends=$HHKinFit_Friends_2017
elif [[ $ERA == *"2018"* ]]
then
    ARTUS_OUTPUTS=$ARTUS_OUTPUTS_2018
    NNScore_Friends=$NNScore_Friends_2018
    SVFit_Friends=$SVFit_Friends_2018
    MELA_Friends=$MELA_Friends_2018
    FF_Friends=$FF_Friends_2018
    HHKinFit_Friends=$HHKinFit_Friends_2018
fi

ARTUS_FRIENDS_FAKE_FACTOR=$FF_Friends
