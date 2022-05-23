ERA=$1

for MASS_H3 in  240 280 320 360 400 450 500 550 600 700 800 900 1000 1200 1400 1600 1800 2000 2500 3000
do


python ../scripts/SkimManager.py --nick "NMSSMM${MASS_H3}h1M125tautauh2M.*_RunIIFall17MiniAODv2_PU2017_13TeV_MINIAOD_madgraph-pythia8_v1" -w /portal/ekpbms2/home/jbechtel/kappa_skim_workdir/nmssm_skims_${ERA}_${MASS_H3} -b freiburg --status-gc 


python ../scripts/SkimManager.py --nick "NMSSMM${MASS_H3}h1M125tautauh2M.*_RunIIFall17MiniAODv2_PU2017_13TeV_MINIAOD_madgraph-pythia8_v1" -w /portal/ekpbms2/home/jbechtel/kappa_skim_workdir/nmssm_skims_${ERA}_${MASS_H3} -b freiburg --create-filelist

done
