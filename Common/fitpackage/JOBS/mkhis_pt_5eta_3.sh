source $VO_CMS_SW_DIR/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc530
export INITDIR=/storage_mnt/storage/user/hanwen/work/CMSSW_8_0_26_patch1/src/shears_old/fit/test
cd $INITDIR
eval `scramv1 runtime -sh`
hostname ;
date;
./makehistos pt=5 eta=3;
mv histos_muons_data_pt5_eta3.root Outputs;
