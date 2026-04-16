#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=reixesreis@gmail.com

#SBATCH --qos=cu_hpc
#SBATCH --partition=cpugpu
#SBATCH --job-name=Wholeprocess
#SBATCH --output=LOG/sbatch_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G

# Source the CMS environment
source /work/app/cms/cmsset_default.sh
function cmsset() {
  cd ~/binary/CMSSW_14_1_0_pre5/src && source /work/app/share_env/hepsw.sh && cmsenv && cd -
}
cmsset

#./run_signal_MuTauE.sh
#./run_backgrounds_MuTauE.sh
#./run_signal_ETauMu.sh
#./run_backgrounds_ETauMu.sh

#python3 makecard.py --in-dir outdir_mutaue --lumi-pb 1000000 --out-root ./datacards_mutaue/merged.root --out-card ./datacards_mutaue ; python3 makecard.py --in-dir outdir_etaumu --lumi-pb 1000000 --out-root ./datacards_etaumu/merged.root --out-card ./datacards_etaumu ; cd datacards_mutaue/ ; sbatch slurm_submit.slurm ; cd ../datacards_etaumu ; sbatch slurm_submit.slurm


# Get bash arguments
while getopts "m:" opt; do
  case ${opt} in
    m )
      SIGNAL_TYPE=$OPTARG
      ;;
    \? )
      echo "Usage: cmd [-m]"
      ;;
  esac
done

#./runall.sh
# python3 -u runall_twoRegion.py -p 32 -o ZH_110To220_integerBin --skip-sbatch --skip-run
if [ "$SIGNAL_TYPE" == "VBF" ]; then
  python3 -u runall_twoRegion_NEW.py -p 32 -o VBF_110To220_integerBin --signal-type VBF --skip-run #--skip-sbatch
else
  python3 -u runall_twoRegion_NEW.py -p 32 -o ZH_110To220_integerBin --signal-type ZH --skip-run #--skip-sbatch
fi
# VBF sample
#python3 -u runall_twoRegion_VBF.py -p 21 -o VBF_165To240_17Mar26 --skip-run

