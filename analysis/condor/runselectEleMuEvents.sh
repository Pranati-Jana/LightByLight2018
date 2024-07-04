#!/bin/bash

userBasePath="/afs/cern.ch/user/p/pjana/private/DiTau2018/analysis/"
#configPath="${userBasePath}/LightByLight2018/analysis/configs/efficiencies_eleNoIsolation.md"
#configPath="${userBasePath}/LightByLight2018/analysis/configs/efficiencies_eleNoIsolation_newThresholds.md"
#configPath="${userBasePath}/LightByLight2018/analysis/configs/efficiencies_HFVeto_withEleIsolatioNewThresholds.md"
#configPath="${userBasePath}/LightByLight2018/analysis/configs/efficiencies_HFVeto_withEleIsolatioNewThresholdsEta2p2.md"
configPath="/afs/cern.ch/user/p/pjana/private/DiTau2018/analysis/configs/efficiencies_eleNoIsolation_newThresholdsEta2p2.md"

inputPath=""
outputPath=""
sampleName=""


basePath="/eos/user/p/pjana/GammaGammaToTauTau"
#basePath="/dpm/indiacms.res.in/home/cms/store/user/pjana/t3store3/Data"

suffix="_ForHFstudy_EleMu"

if [ $2 -eq 0 ]
then
  sampleName="Data" # 10400 files
#  inputPath=`sed "${1}q;d" ${userBasePath}/LightByLight2018/analysis/input_list.txt`
  #inputPath=`sed "${1}q;d" ${userBasePath}/LightByLight2018/analysis/input_list_Vtx_Pranati.txt`
  inputPath=`sed "${1}q;d" /eos/cms/store/cmst3/group/lightbylight/Pranati/HIForest_allConversions/HIForward/ntuples_data_ext/230618_233239/Data_HIForest_2018.txt`
  #inputPath=`sed "${1}q;d" /eos/cms/store/group/phys_heavyions/pjana/Pranati_ForDiTau_Analysis/Data_OneFourthLumi/HIForward/ntuples_data_ext/240427_053549/Data_OnefourthLumi.txt`
  outputPath="/eos/user/p/pjana/RootFileTauTau/AllSelections_FromSelectEleMuScript/CutFlowHistAdded_18thMarch/Data_FullLumi_NO_NE${suffix}"
elif [ $2 -eq 1 ]
then 
  sampleName="QED_SC"
  #inputPath="${basePath}/skimmed_ntuples/mc_qedMG5_FSR_doubleEG2_full_lumi/merged/merged_ntuples_${1}.root"
  #inputPath="${basePath}/ntuples/MC_ntuples/GammaUPC/12thMay_CH_kTSmear/GammaGammatoTauTau_5p02TeV_gammaUPCChFFkTSmearing-pLHE-v1/reco_mc_SC/230512_141327/0000/mc_HiForestAOD_${1}.root"
  #inputPath="${basePath}/ntuples/MC_ntuples/GammaUPC/21stJune_SignalMC/GammaGammatoTauTau_5p02TeV_gammaUPCChFFkTSmearing-pLHE-v1/reco_mc_SC/230621_144737/0000/mc_HiForestAOD_${1}.root"
 # inputPath="${basePath}/ntuples/MC_ntuples/GammaUPC/9thJuly_TauGen_daughterinfo/GammaGammatoTauTau_5p02TeV_gammaUPCChFFkTSmearing-pLHE-v1/reco_mc_SC/230710_223838/0000/mc_HiForestAOD_${1}.root"
#  inputPath= "/eos/cms/store/cmst3/group/lightbylight/Pranati/Pranati_ForDiTau_Analysis/QED_tautau_atau0_TuneCP5_5p02TeV_UPCgen_pythia8/reco_mc_UPCgen/240124_173810/0000/mc_HiForestAOD_${1}.root" 
 #outputPath="${basePath}/ElectronMuon/MC_Signal_Post_Selection/UPCGamma/mc_GammaUPC${suffix}"
  #outputPath="/eos/user/p/pjana/GammaGammaToTauTau/ElectronMuon/MC_Signal_Post_Selection/UPCGamma/mc_GammaUPC${suffix}"
#  inputPath="/eos/cms/store/cmst3/group/lightbylight/Pranati/Pranati_ForDiTau_Analysis/QED_tautau_atauPlus6_TuneCP5_5p02TeV_UPCgen_pythia8/reco_mc_UPCgen/240126_042154/0000/mc_HiForestAOD_${1}.root" 
inputPath="/eos/user/p/pjana/GammaGammaToTauTau/ntuples/MC_ntuples/SC_Arash/MC_SCArash_GenInfo_2ndAugust/ggTauTau_TuneCP5_5p02TeV_SuperChic_pythia8/reco_mc_SC/230802_053701/0000/mc_HiForestAOD_${1}.root" 
#
outputPath="/eos/user/p/pjana/RootFileTauTau/AllSelections_FromSelectEleMuScript/CutFlowHistAdded_18thMarch/mc_SC${suffix}"
elif [ $2 -eq 2 ]
then
  sampleName="QED_SL"
  inputPath="${basePath}/skimmed_ntuples/mc_qed_sl_doubleEG2_full_lumi/ntuples_loose_selections_${1}.root"
  outputPath="${basePath}/analysis_ruchi/qed_data-MC_plots/mc_qed_sl${suffix}"
elif [ $2 -eq 3 ]
then
  sampleName="GAMMA_UPC"

#inputPath="/eos/user/p/pjana/GammaGammaToTauTau/ntuples/MC_ntuples/GammaUPC/MC_GAMMA_UPC_GenInfo_23rdJuly_nDaughterInfo/GammaGammatoTauTau_5p02TeV_gammaUPCChFFkTSmearing-pLHE-v1/reco_mc_SC/230724_092329/0000/mc_HiForestAOD_${1}.root"
inputPath="/eos/cms/store/group/phys_heavyions/pjana/Pranati_ForDiTau_Analysis/a_TauPoints/QED_tautau_atau0_TuneCP5_5p02TeV_UPCgen_pythia8/reco_mc_UPCgen/240124_173810/0000/mc_HiForestAOD_${1}.root"
#inputPath="/eos/cms/store/group/phys_heavyions/pjana/Pranati_ForDiTau_Analysis/a_TauPoints/QED_tautau_atauPlus1_TuneCP5_5p02TeV_UPCgen_pythia8/reco_mc_UPCgen/240226_094822/0000/mc_HiForestAOD_${1}.root"
#inputPath="/eos/cms/store/group/phys_heavyions/pjana/Pranati_ForDiTau_Analysis/a_TauPoints/QED_tautau_atauPlus2_TuneCP5_5p02TeV_UPCgen_pythia8/reco_mc_UPCgen/240226_151939/0000/mc_HiForestAOD_${1}.root"
#inputPath="/eos/cms/store/group/phys_heavyions/pjana/Pranati_ForDiTau_Analysis/a_TauPoints/QED_tautau_atauPlus3_TuneCP5_5p02TeV_UPCgen_pythia8/reco_mc_UPCgen/240227_163759/0000/mc_HiForestAOD_${1}.root"
#inputPath="/eos/cms/store/group/phys_heavyions/pjana/Pranati_ForDiTau_Analysis/a_TauPoints/QED_tautau_atauPlus4_TuneCP5_5p02TeV_UPCgen_pythia8/reco_mc_UPCgen/240227_070118/0000/mc_HiForestAOD_${1}.root"
#inputPath="/eos/cms/store/group/phys_heavyions/pjana/Pranati_ForDiTau_Analysis/a_TauPoints/QED_tautau_atauPlus5_TuneCP5_5p02TeV_UPCgen_pythia8/reco_mc_UPCgen/240228_044851/0000/mc_HiForestAOD_${1}.root"
#inputPath="/eos/cms/store/group/phys_heavyions/pjana/Pranati_ForDiTau_Analysis/a_TauPoints/QED_tautau_atauPlus6_TuneCP5_5p02TeV_UPCgen_pythia8/reco_mc_UPCgen/240316_185055/0000/mc_HiForestAOD_${1}.root"
#inputPath="/eos/cms/store/group/phys_heavyions/pjana/Pranati_ForDiTau_Analysis/a_TauPoints/QED_tautau_atauPlus7_TuneCP5_5p02TeV_UPCgen_pythia8/reco_mc_UPCgen/240126_091217/0000/mc_HiForestAOD_${1}.root"
#inputPath="/eos/cms/store/group/phys_heavyions/pjana/Pranati_ForDiTau_Analysis/a_TauPoints/QED_tautau_atauPlus8_TuneCP5_5p02TeV_UPCgen_pythia8/reco_mc_UPCgen/240126_134315/0000/mc_HiForestAOD_${1}.root"
#inputPath="/eos/cms/store/group/phys_heavyions/pjana/Pranati_ForDiTau_Analysis/a_TauPoints/QED_tautau_atauPlus9_TuneCP5_5p02TeV_UPCgen_pythia8/reco_mc_UPCgen/240126_185217/0000/mc_HiForestAOD_${1}.root"
#inputPath="/eos/cms/store/group/phys_heavyions/pjana/Pranati_ForDiTau_Analysis/a_TauPoints/QED_tautau_atauPlus10_TuneCP5_5p02TeV_UPCgen_pythia8/reco_mc_UPCgen/240127_045202/0000/mc_HiForestAOD_${1}.root"
#inputPath="/eos/cms/store/group/phys_heavyions/pjana/Pranati_ForDiTau_Analysis/a_TauPoints/QED_tautau_atauMinus1_TuneCP5_5p02TeV_UPCgen_pythia8/reco_mc_UPCgen/240127_124533/0000/mc_HiForestAOD_${1}.root"
#inputPath="/eos/cms/store/group/phys_heavyions/pjana/Pranati_ForDiTau_Analysis/a_TauPoints/QED_tautau_atauMinus3_TuneCP5_5p02TeV_UPCgen_pythia8/reco_mc_UPCgen/240128_022957/0000/mc_HiForestAOD_${1}.root"
#inputPath="/eos/cms/store/group/phys_heavyions/pjana/Pranati_ForDiTau_Analysis/a_TauPoints/QED_tautau_atauMinus4_TuneCP5_5p02TeV_UPCgen_pythia8/reco_mc_UPCgen/240128_122502/0000/mc_HiForestAOD_${1}.root"
#inputPath="/eos/cms/store/group/phys_heavyions/pjana/Pranati_ForDiTau_Analysis/a_TauPoints/QED_tautau_atauMinus10_TuneCP5_5p02TeV_UPCgen_pythia8/reco_mc_UPCgen/240131_083928/0000/mc_HiForestAOD_${1}.root"
outputPath="/eos/user/p/pjana/RootFileTauTau/AllSelections_FromSelectEleMuScript/CutFlowHistAdded_18thMarch/gammaUPC_${suffix}"
elif [ $2 -eq 4 ]
then
  sampleName="MUMU_FSR"
  inputPath="/eos/user/p/pjana/RootFileTauTau/ntuples/24thAugust_mumuFSR_Arash_gUPC/gammaUPCmumuFSR/reco_mc_SC/230824_094901/0000/mc_HiForestAOD_${1}.root"
  #outputPath="/eos/user/p/pjana/RootFileTauTau/AllSelections_FromSelectEleMuScript/mc_mumuFSR${suffix}"
  outputPath="/eos/user/p/pjana/RootFileTauTau/AllSelections_FromSelectEleMuScript/CutFlowHistAdded_18thMarch/mc_mumuFSR${suffix}"

fi

mkdir -p $outputPath
output="${outputPath}/EleMu_${1}.root"

echo "Input: ${inputPath}"

if [ -s ${output} ]
then
  echo "File already exists, skipping"
else
  echo "File doesn't exist or is empty - running"
  ${userBasePath}/selectEleMuEvents $configPath $inputPath $output $sampleName
fi
