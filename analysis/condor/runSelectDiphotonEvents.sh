#!/bin/bash

#configPath="/afs/cern.ch/work/r/rchudasa/private/LightByLight2018/analysis/configs/efficiencies_eleNoIsolation.md"
#configPath="/afs/cern.ch/work/r/rchudasa/private/LightByLight2018/analysis/configs/efficiencies_eleNoIsolation_newThresholds.md"
#configPath="/afs/cern.ch/work/j/jniedzie/private/LightByLight2018/analysis/configs/efficiencies_eleNoIsolation_newThresholdsEta2p2.md"
configPath="/afs/cern.ch/user/p/pjana/private/LightByLight/LByL/LightByLight2018/analysis/configs/efficiencies_eleNoIsolation_newThresholdsEta2p2.md"

inputPath=""
outputPath=""
sampleName=""
#Ruchi
basePath="/eos/cms/store/group/phys_heavyions/rchudasa/lbyl_2018"
#Pranati
#basePath="/eos/cms/store/cmst3/group/lightbylight/Final_afterTrigger/Data_29thJune/merged"
#basePathOutput="/eos/cms/store/group/phys_heavyions/rchudasa/lbyl_2018/analysis_ruchi/diphoton_data-MC_plots"
basePathOutput="/eos/user/p/pjana/LByL/Test_Ruchintuple/31stJan_TriggerDiphoton"
suffix="_diphoton_latestNtuple"

if [ $2 -eq 0 ]
then
  sampleName="Data" # 10400 files
 # inputPath="${basePath}/merged_ntuples_${1}.root"
  inputPath=`sed "${1}q;d" /eos/cms/store/cmst3/group/lightbylight/Pranati/HIForest_allConversions/HIForward/ntuples_data_ext/230618_233239/Data_HIForest_2018.txt` 
 #Vettex info
# inputPath="/eos/cms/store/cmst3/group/lightbylight/Pranati/vtxInfo_Pranati/data_doubleEG2_full_lumi_11thFeb/merged/merged_ntuples_${1}.root"
#  inputPath="/eos/cms/store/cmst3/group/lightbylight/Pranati/Final_afterTrigger/Data/Data_29thJune/merged/merged_ntuples_${1}.root"
  #EmptyBX
  #inputPath="${basePath}/HIEmptyBX/ntuples/HIEmptyBX/ntuples_emptyBx_castor_pixel_branches/210803_040945/0000/data_emptyBx_HiForwardAOD_${1}.root"
  #inputPath=`sed "${1}q;d" /afs/cern.ch/work/r/rchudasa/private/LightByLight2018/analysis/input_list.txt`
 # inputPath=`sed "${1}q;d" /afs/cern.ch/user/p/pjana/private/LightByLight/LByL/LightByLight2018/analysis/input_list.txt`
  # SD Muon
 # inputPath="${basePath}/HIForward/ntuples_data_ext/221121_170303/0000/data_HiForwardAOD_${1}.root"
 # Standalone Muon, Pranati's ntuple  
  #inputPath=`sed "${1}q;d" /afs/cern.ch/user/p/pjana/private/LightByLight/LByL/LightByLight2018/analysis/input_list_STMuon.txt`
#  inputPath=`sed "${1}q;d" /eos/cms/store/cmst3/group/lightbylight/Pranati/HIForest_allConversions/HIForward/ntuples_data_ext/230618_233239/Data_HIForest_2018.txt`
  #outputPath="${basePathOutput}/data${suffix}"
 # outputPath="/eos/user/p/pjana/LByL/Final_Data/Data/Data${suffix}"
  outputPath="/eos/cms/store/cmst3/group/lightbylight/Pranati/CutFlow_study/Data_TriggerWithDiphoton${suffix}"
elif [ $2 -eq 1 ]
then
  sampleName="QED_SC" # last chunk numer: 255
  #inputPath="${basePath}/skimmed_ntuples/mc_qed_sc_doubleEG2_full_lumi/ntuples_loose_selections_${1}.root"
  #inputPath="${basePath}/skimmed_ntuples/mc_qed_sc_FSR_doubleEG2_full_lumi/merged/merged_ntuples_${1}.root"
  #inputPath="${basePath}/skimmed_ntuples/mc_qedMG5_FSR_doubleEG2_full_lumi/merged/merged_ntuples_${1}.root"
  #inputPath="${basePath}/skimmed_ntuples/mc_qedMG5_2FSR_doubleEG2_full_lumi/merged/merged_ntuples_${1}.root"
  #TauTau Check
  #inputPath="/eos/user/p/pjana/LByL/TauTau_files_afterApplySelection/mc_TauTau_SC__16thFeb_Arash/ntuples_loose_selections_${1}.root"
  #inputPath="/eos/user/p/pjana/LByL/TauTau_files_afterApplySelection/mc_TauTau_GammaUPC_UPCEDFFkTSmearing_26thFeb_GK/ntuples_loose_selections_${1}.root"
  #inputPath="/eos/user/p/pjana/LByL/TauTau_files_afterApplySelection/GammaUPC/ntuples_loose_selections_${1}.root"
#  inputPath="/eos/user/p/pjana/LByL/TauTau_files_afterApplySelection/GammaUPC/afterTrigger/GammaGammatoTauTau_5p02TeV_gammaUPCEDFFkTSmearing-pLHE-v1_2/ntuples_loose_selections_${1}.root"
#inputPath="/eos/user/p/pjana/LByL/TauTau_files_afterApplySelection/mc_TauTau_SC__15thFeb_OfficialGk/ntuples_loose_selections_${1}.root" 
   inputPath="/eos/cms/store/cmst3/group/lightbylight/Pranati/Final_afterTrigger/QED_SC_Photos/ntuples_loose_selections_${1}.root" 
  #inputPath="/eos/cms/store/cmst3/group/lightbylight/Pranati/HIForest_QED_Photos_New_SC/24thJuly/QED_Photosplpl_5p02TeV_Superchic3/reco_mc_SC/230724_190910/0000/mc_HiForestAOD_${1}.root"
  outputPath="/eos/cms/store/cmst3/group/lightbylight/Pranati/CutFlow_study/TriggeredPassed_DiphotonEvents/QED_PHOTOS_${suffix}"
  #outputPath="/eos/cms/store/cmst3/group/lightbylight/Pranati/CutFlow_study/FirstTriggerToDiPhoSelectection/QED_PHOTOS_${suffix}"
  #outputPath="${basePathOutput}/mc_qedMG5_2FSR${suffix}"
  #TauTau Check
  #outputPath="/eos/user/p/pjana/LByL/TauTau_files_afterApplySelection/afterDiphoton/mc_TauTau_SC${suffix}"
  #outputPath="/eos/user/p/pjana/LByL/QED_POTOS_Official/Diphoton_Selection${suffix}"
  #outputPath="/eos/user/p/pjana/LByL/TauTau_files_afterApplySelection/GammaUPC/afterTrigger_diphoton/mc_TauTau_GammaUPC${suffix}"
elif [ $2 -eq 2 ]
then
  sampleName="QED_SL" # last chunk numer: 253
  #inputPath="${basePath}/mc_qed/ntuples_sl_full_lumi_v5/QEDGammaGamma_5p02TeV_STARlight/reco_mc_qed_sl_full_lumi_v5/200929_094304/0000/mc_HiForestAOD_${1}.root"
  inputPath="${basePath}/skimmed_ntuples/mc_qed_sl_doubleEG2_full_lumi/ntuples_loose_selections_${1}.root"
  outputPath="${basePathOutput}/mc_qed_sl${suffix}"
elif [ $2 -eq 3 ]
then
  sampleName="LbL" # 3 files only
  inputPath="${basePath}/skimmed_ntuples/mc_lbl_sc_doubleEG2_full_lumi/ntuples_loose_selections_${1}.root"
#  inputPath="$/eos/cms/store/group/phys_diffraction/lbyl_2018/skimmed_ntuples/mc_lbl_sc_doubleEG2_full_lumi/ntuples_loose_selections_${1}.root"
  #inputPath="/eos/user/p/pjana/LByL/LByL_ZS_Test_2023/after_TriggerOnly/19thAugust/merged/merged_ntuples_${1}.root" 
 # inputPath="/eos/user/p/pjana/LByL/LByL_ZS_Test_2023/after_TriggerOnly/ZS_lblConfig/merged/merged_ntuples_${1}.root"
 # inputPath="/eos/user/p/pjana/LByL/LByL_ZS_Test_2023/after_TriggerOnly/Default_LByL2018MC_GK_ForCheck/merged/merged_ntuples_${1}.root"
#  inputPath="/eos/user/p/pjana/LByL/LByL_ZS_Test_2023/after_TriggerOnly/Default_LbyL2018MC_25thAugust_GK/ntuples_loose_selections_${1}.root"
# inputPath="/eos/user/p/pjana/LByL/LByL_ZS_Test_2023/Official_LByL_MC_Default_20thAugust/LbyLSignal_5p02TeV_SuperChic/reco_mc_lbl_try2/230821_061454/0000/mc_HiForestAOD_${1}.root"
# inputPath="/eos/cms/store/cmst3/group/lightbylight/Pranati/Final_afterTrigger/LByLMC_SC//ntuples_loose_selections_${1}.root"
 outputPath="/eos/cms/store/cmst3/group/lightbylight/Pranati/CutFlow_study/TriggeredPassed_DiphotonEvents/mc${suffix}" 
# outputPath="/eos/cms/store/cmst3/group/lightbylight/Pranati/CutFlow_study/mc${suffix}"
#  outputPath="/eos/user/p/pjana/LByL/LByL_ZS_Test_2023/after_DiphotonSelection/mc${suffix}"
  #outputPath="/eos/user/p/pjana/LByL/LByL_ZS_Test_2023/after_DiphotonSelection/mc${suffix}"
 #outputPath="${basePathOutput}/mc_lbl_sc${suffix}"
  #outputPath="/eos/user/p/pjana/LByL/LByL_ZS_Test_2023/after_DiphotonSelection/mc_SC_ZS_2023${suffix}"
  #outputPath="/eos/user/p/pjana/LByL/LByL_ZS_Test_2023/after_DiphotonSelection/mc_SC_ZS_2023_lblConfig${suffix}"
#  outputPath="/eos/user/p/pjana/LByL/LBLMC_afterDiphotonSelection/mc_lbl_sc${suffix}"
elif [ $2 -eq 4 ]
then
  sampleName="CEP" # 3 files only
 inputPath="${basePath}/skimmed_ntuples/mc_cep_sc_doubleEG2_full_lumi/ntuples_loose_selections_${1}.root"
   
  outputPath="${basePathOutput}/mc_cep_sc${suffix}"
elif [ $2 -eq 5 ]
then
  sampleName="CEP" # 3 files only
  inputPath="${basePath}/mcForests/mc_alps/AxionLikeParticles_M-5_5p02TeV_SuperChic/mc_alps5GeV/211208_073352/0000/mc_HiForestAOD_${1}.root"
  outputPath="${basePathOutput}/mc_alps_5GeV${suffix}"
elif [ $2 -eq 6 ]
then
  sampleName="CEP" # 3 files only
  inputPath="${basePath}/mcForests/mc_alps/AxionLikeParticles_M-6_5p02TeV_SuperChic/mc_alps6GeV/211208_081512/0000/mc_HiForestAOD_${1}.root"
  outputPath="${basePathOutput}/mc_alps_6GeV${suffix}"
elif [ $2 -eq 9 ]
then
  sampleName="CEP" # 3 files only
  inputPath="${basePath}/mcForests/mc_alps/AxionLikeParticles_M-9_5p02TeV_SuperChic/mc_alps9GeV/211208_081631/0000/mc_HiForestAOD_${1}.root"
  outputPath="${basePathOutput}/mc_alps_9GeV${suffix}"
elif [ $2 -eq 11 ]
then
  sampleName="CEP" # 3 files only
  inputPath="${basePath}/mcForests/mc_alps/AxionLikeParticles_M-11_5p02TeV_SuperChic/mc_alps11GeV/211208_081715/0000/mc_HiForestAOD_${1}.root"
  outputPath="${basePathOutput}/mc_alps_11GeV${suffix}"
elif [ $2 -eq 14 ]
then
  sampleName="CEP" # 3 files only
  inputPath="${basePath}/mcForests/mc_alps/AxionLikeParticles_M-14_5p02TeV_SuperChic/mc_alps14GeV/211208_081746/0000/mc_HiForestAOD_${1}.root"
  outputPath="${basePathOutput}/mc_alps_14GeV${suffix}"
elif [ $2 -eq 16 ]
then
  sampleName="CEP" # 3 files only
  inputPath="${basePath}/mcForests/mc_alps/AxionLikeParticles_M-16_5p02TeV_SuperChic/mc_alps16GeV/211208_081833/0000/mc_HiForestAOD_${1}.root"
  outputPath="${basePathOutput}/mc_alps_16GeV${suffix}"
elif [ $2 -eq 22 ]
then
  sampleName="CEP" # 3 files only
  inputPath="${basePath}/mcForests/mc_alps/AxionLikeParticles_M-22_5p02TeV_SuperChic/mc_alps22GeV/211208_081937/0000/mc_HiForestAOD_${1}.root"
  outputPath="${basePathOutput}/mc_alps_22GeV${suffix}"
elif [ $2 -eq 30 ]
then
  sampleName="CEP" # 3 files only
  inputPath="${basePath}/mcForests/mc_alps/AxionLikeParticles_M-30_5p02TeV_SuperChic/mc_alps30GeV/211208_082015/0000/mc_HiForestAOD_${1}.root"
  outputPath="${basePathOutput}/mc_alps_30GeV${suffix}"
elif [ $2 -eq 50 ]
then
  sampleName="CEP" # 3 files only
  inputPath="${basePath}/mcForests/mc_alps/AxionLikeParticles_M-50_5p02TeV_SuperChic/mc_alps50GeV/211208_082053/0000/mc_HiForestAOD_${1}.root"
  outputPath="${basePathOutput}/mc_alps_50GeV${suffix}"
elif [ $2 -eq 90 ]
then
  sampleName="CEP" # 3 files only
  inputPath="${basePath}/mcForests/mc_alps/AxionLikeParticles_M-90_5p02TeV_SuperChic/mc_alps90GeV/211208_082354/0000/mc_HiForestAOD_${1}.root"
  outputPath="${basePathOutput}/mc_alps_90GeV${suffix}"
fi


mkdir -p $outputPath
output="${outputPath}/diphoton_${1}.root"

if [ -s ${output} ]
then
  echo "File already exists, skipping"
else
  echo "File doesn't exist or is empty - running"
  /afs/cern.ch/user/p/pjana/private/LightByLight/LByL/LightByLight2018/analysis/selectDiphotonEvents $configPath $inputPath $output $sampleName
fi
