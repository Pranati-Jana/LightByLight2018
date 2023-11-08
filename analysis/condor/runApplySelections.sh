#!/bin/bash

userBasePath="/afs/cern.ch/user/p/pjana/private/LightByLight/LByL/"
configPath="${userBasePath}/LightByLight2018/analysis/configs/efficiencies_eleNoIsolation_newThresholdsEta2p2.md"
basePath="/eos/cms/store/group/phys_diffraction/lbyl_2018"
#basePath="/eos/user/p/pjana/LByL/ntuples/QED_SC_Photos_New_24thJuly"
sampleName=""

# for the data:
inputPath=""
outputPathReco=""
outputPathTrigger=""
outputPathHFveto=""
outputPathExclusivity=""
outputPathLbLsignal=""
outputPathQEDsignal=""
outputPathLowAco=""
outputPathHighAco=""
outputPath=""

echo "Params: ${1}, ${2}, ${3}"

if [ $2 -eq 0 ] # data
then
  sampleName="Data"
  inputPath=`sed "${1}q;d" ${userBasePath}/LightByLight2018/analysis/input_list.txt`
  outputPathReco="/eos/cms/store/group/phys_diffraction/lbyl_2018/skimmed_ntuples/data_forRecoEff_tmp"
  outputPathTrigger="/eos/cms/store/group/phys_diffraction/lbyl_2018/skimmed_ntuples/data_forTriggerEff_tmp"
  outputPathHFveto="/eos/cms/store/group/phys_diffraction/lbyl_2018/skimmed_ntuples/data_forHFvetoEff_tmp"
  outputPathExclusivity="/eos/cms/store/group/phys_diffraction/lbyl_2018/skimmed_ntuples/data_forExclusivityEff_tmp"
  outputPathLbLsignal="/eos/cms/store/group/phys_diffraction/lbyl_2018/skimmed_ntuples/data_forLbLsignal_CHE_500MeV"
  outputPathQEDsignal="/eos/cms/store/group/phys_diffraction/lbyl_2018/skimmed_ntuples/data_forQEDsignal_tmp"
  mkdir -p $outputPathReco
  mkdir -p $outputPathTrigger
  mkdir -p $outputPathHFveto
  mkdir -p $outputPathExclusivity
  mkdir -p $outputPathLbLsignal
  mkdir -p $outputPathQEDsignal
elif [ $2 -eq 1 ] # MC
then
  # QED
#   sampleName="QED_SC"
#  inputPath="/eos/cms/store/group/phys_diffraction/lbyl_2018/mc_qed/ntuples_superchic_1034/ntuples_sc_1034/ntuples_sc_1034/191113_105005/0000/HiForestAOD_LbyL_full_sample_lbyl_reco_${1}.root"
  
  # LbL
#  sampleName="LbL"
#  inputPath="/eos/cms/store/group/phys_diffraction/lbyl_2018/mc_lbl/ntuples_1034/ntuples_lbl_1034/ntuples_lbl_1034/200207_114802/0000/HiForestAOD_LbyL_${1}.root"
  
  # CEP
  sampleName="CEP" inputPath="/eos/cms/store/group/phys_diffraction/lbyl_2018/mc_cep/ntuples_1034/ntuples_cep_1034/ntuples_cep_1034/200211_054704/0000/HiForestAOD_cep_${1}.root"
  
  outputPathReco="/eos/cms/store/group/phys_diffraction/lbyl_2018/mc_cep_sc_forRecoEff"
  outputPathTrigger="/eos/cms/store/group/phys_diffraction/lbyl_2018/mc_cep_sc_forTriggerEff"
  outputPathHFveto="/eos/cms/store/group/phys_diffraction/lbyl_2018/mc_cep_sc_forHFvetoEff"
  outputPathExclusivity="/eos/cms/store/group/phys_diffraction/lbyl_2018/mc_cep_sc_forExclusivityEff"
  outputPathLbLsignal="/eos/cms/store/group/phys_diffraction/lbyl_2018/mc_cep_sc_forLbLsignal"
  outputPathQEDsignal="/eos/cms/store/group/phys_diffraction/lbyl_2018/mc_cep_sc_forQEDsignal"
  mkdir -p $outputPathReco
  mkdir -p $outputPathTrigger
  mkdir -p $outputPathHFveto
  mkdir -p $outputPathExclusivity
  mkdir -p $outputPathLbLsignal
  mkdir -p $outputPathQEDsignal
elif [ $2 -eq 2 ] # Data passing LbL selections
then
  sampleName="Data"
  inputPath=`sed "${1}q;d" ${userBasePath}/LightByLight2018/analysis/input_list.txt`
  outputPathLowAco="/eos/cms/store/group/phys_diffraction/lbyl_2018/skimmed_ntuples/data_passingLbL_lowAco"
  outputPathHighAco="/eos/cms/store/group/phys_diffraction/lbyl_2018/skimmed_ntuples/data_passingLbL_highAco"
  mkdir -p $outputPathLowAco
  mkdir -p $outputPathHighAco
elif [ $2 -eq 3 ] # Loose selections
then
  echo "Applying loose selections, ${3}"

  if [ ${3} -eq 0 ] # data, 10400 chunks
  then
    echo "Data"
    sampleName="Data"
    inputPath=`sed "${1}q;d" ${userBasePath}/LightByLight2018/analysis/input_list_withAllConversions.txt`
    #inputPath=""
    #outputPath="${basePath}/skimmed_ntuples/data_doubleEG2_full_lumi"
    #outputPath="/eos/cms/store/cmst3/group/lightbylight/Test_Ruchi/data_doubleEG2_full_lumi_13thFeb"
    outputPath="/eos/cms/store/cmst3/group/lightbylight/Final_afterTrigger/Data_29thJune"
  elif [ ${3} -eq 1 ] # QED SC, 255 chunks, max chunk number: 255
  then
    echo "QED SC"
    sampleName="QED_SC"
    #inputPath="${basePath}/mc_qed/ntuples_superchic_1034/ntuples_sc_1034/ntuples_sc_1034/191113_105005/0000/HiForestAOD_LbyL_full_sample_lbyl_reco_${1}.root"
    #inputPath="${basePath}/mc_qed/ntuples_sc_full_lumi/QEDGammaGamma_5p02TeV_SuperChic/reco_mc_qed_sc_full_lumi/200807_100412/0000/mc_HiForestAOD_${1}.root"
    #inputPath="${basePath}/mc_qed/QEDGammaGamma_5p02TeV_SuperChic/reco_mc_qed_sc_full_lumi/210417_081453/0000/mc_HiForestAOD_${1}.root"
    #inputPath="${basePath}/mcForests/mc_qed/QEDGammaGamma_5p02TeV_SuperChic/reco_mc_qed_sc_CastorInfo/210906_102658/0000/mc_HiForestAOD_${1}.root"
    #inputPath=`sed "${1}q;d" ${userBasePath}/LightByLight2018/analysis/input_qed_fsr.txt`
    #inputPath=`sed "${1}q;d" ${userBasePath}/LightByLight2018/analysis/qedMG5.txt`
    #inputPath=`sed "${1}q;d" ${userBasePath}/LightByLight2018/analysis/qedMG5_2FSR.txt`
    #TauTau Superchic sample
    #inputPath="/eos/user/p/pjana/GammaGammaToTauTau/ntuples/MC_ntuples/MC_ntuples_Arash/ggTauTau_TuneCP5_5p02TeV_SuperChic_pythia8/reco_mc_SC/220707_112036/0000/mc_HiForestAOD_${1}.root"  #Date:21/01/2023
    #GammaUPC sample
    # inputPath="/eos/user/p/pjana/GammaGammaToTauTau/ntuples/MC_ntuples/GammaUPC/24thFeb_UPCEDFFkTSmearing/GammaGammatoTauTau_5p02TeV_gammaUPCEDFFkTSmearing-pLHE-v1/reco_mc_SC/230226_184903/0000/mc_HiForestAOD_${1}.root"
    #inputPath="/eos/user/p/pjana/LByL/TauTau_files_afterApplySelection/GammaUPC/GammaGammatoTauTau_5p02TeV_gammaUPCEDFFkTSmearing-pLHE-v1/reco_mc_SC/230507_145915/0000/mc_HiForestAOD_${1}.root"
#    inputPath="/eos/cms/store/cmst3/group/lightbylight/HIForest_QED_Photos_New_SC/24thJuly/QED_Photosplpl_5p02TeV_Superchic3/reco_mc_SC/230724_190910/0000/mc_HiForestAOD_${1}.root"
    inputPath=`sed "${1}q;d" /eos/user/p/pjana/LByL/LByL_ZS_Test_2023/19thAugust/LbyLSignal_5p02TeV_SuperChic-v2/reco_mc_SC/230819_085059/LByL_New_ZS_ECAL_Th.txt`
    #New ZS ECAL Th 2023
    outputPath="/eos/user/p/pjana/LByL/LByL_ZS_Test_2023/after_TriggerOnly/19thAugust"
    #outputPath="${basePath}/skimmed_ntuples/mc_qed_sc_FSR_doubleEG2_full_lumi"
    #outputPath="${basePath}/skimmed_ntuples/mc_qedMG5_2FSR_doubleEG2_full_lumi"
    #outputPath="/eos/user/p/pjana/LByL/GammaUPC"
    #TauTau SC
    #outputPath="/eos/user/p/pjana/LByL/TauTau_files_afterApplySelection/GammaUPC/afterTrigger/GammaGammatoTauTau_5p02TeV_gammaUPCEDFFkTSmearing-pLHE-v1_2/"
    #outputPath="/eos/cms/store/cmst3/group/lightbylight/Final_afterTrigger/QED_SC_Photos"
  elif [ ${3} -eq 2 ] # CEP SC, 3 chunks, max chunk number: 3
  then
    echo "QED SL"
    sampleName="QED_SL"
    #inputPath="${basePath}/mc_qed/ntuples_sl_full_lumi/QEDGammaGamma_5p02TeV_STARlight/reco_mc_qed_sl_full_lumi/200702_082621/0000/mc_HiForestAOD_${1}.root"
    #inputPath="${basePath}/mc_qed/QEDGammaGamma_5p02TeV_STARlight/reco_mc_qed_sl_full_lumi/210417_080949/0000/mc_HiForestAOD_${1}.root"
    inputPath="${basePath}/mcForests/mc_qed/QEDGammaGamma_5p02TeV_STARlight/reco_mc_qed_sl_CastorInfo/210906_102343/0000/mc_HiForestAOD_${1}.root"
    outputPath="${basePath}/skimmed_ntuples/mc_qed_sl_doubleEG2_full_lumi"
  elif [ ${3} -eq 3 ] # LbL SC, 3 chunks, max chunk number: 3
  then
    echo "LbL"
    sampleName="LbL"
    #inputPath="${basePath}/mc_lbl/ntuples_1034/ntuples_lbl_1034/ntuples_lbl_1034/200207_114802/0000/HiForestAOD_LbyL_${1}.root"
    #inputPath="${basePath}/mc_lbl/LbyLSignal_5p02TeV_SuperChic/reco_mc_lbl_try2/210420_063700/0000/mc_HiForestAOD_${1}.root"
    #inputPath="${basePath}/mcForests/mc_lbl/LbyLSignal_5p02TeV_SuperChic/reco_mc_lbl_CastorInfo/210906_102904/0000/mc_HiForestAOD_${1}.root"
    #inputPath=`sed "${1}q;d" /eos/user/p/pjana/LByL/LByL_ZS_Test_2023/19thAugust/LbyLSignal_5p02TeV_SuperChic-v2/reco_mc_SC/230819_085059/LByL_New_ZS_ECAL_Th.txt`
    #inputPath= "/eos/user/p/pjana/LByL/LByL_ZS_Test_2023/Official_LByL_MC_Default_20thAugust/LbyLSignal_5p02TeV_SuperChic/reco_mc_lbl_try2/230820_192436/0000/mc_HiForestAOD_${1}.root"
   #inputPath=`sed "${1}q;d" /eos/user/p/pjana/LByL/LByL_ZS_Test_2023/ZS_21stAugust_lblConfig/LbyLSignal_5p02TeV_SuperChic-v2/reco_mc_lbl_try2/230821_045152/ZS_lblConfig.txt` 
   #inputPath="/eos/user/p/pjana/LByL/LByL_ZS_Test_2023/Official_LByL_MC_Default_20thAugust/LbyLSignal_5p02TeV_SuperChic/reco_mc_lbl_try2/230821_061454/0000/mc_HiForestAOD_${1}.root"
    #inputPath=`sed "${1}q;d" /eos/user/p/pjana/LByL/LByL_ZS_Test_2023/Default_LByL2018MC_GK_ForCheck/LbyLSignal_5p02TeV_SuperChic-v2/reco_mc_lbl_try2/230821_134937/Default_LByL2018MC_GK_ForCheck.txt`
    #inputPath="/eos/user/p/pjana/LByL/LByL_ZS_Test_2023/Default_LByL2018MC_GK_ForCheck_2/LbyLSignal_5p02TeV_SuperChic-v2/reco_mc_lbl_try2/230821_180224/0000/mc_HiForestAOD_${1}.root"
    #inputPath="/eos/user/p/pjana/LByL/LByL_ZS_Test_2023/Default_LbyL2018MC_25thAugust_GK/LbyLSignal_5p02TeV_SuperChic-v4/reco_mc_lbl_try2/230825_052534/0000/mc_HiForestAOD_${1}.root"
    #outputPath="/eos/user/p/pjana/LByL/LByL_ZS_Test_2023/after_TriggerOnly/Default_LbyL2018MC_25thAugust_GK"
    #inputPath="/eos/user/p/pjana/LByL/LByL_ZS_Test_2023/ZS_2023_GK_28thAugust_Realistic/LbyLSignal_5p02TeV_SuperChic-v4/reco_mc_lbl_try2/230828_123806/0000/mc_HiForestAOD_${1}.root"
    #outputPath="/eos/user/p/pjana/LByL/LByL_ZS_Test_2023/after_TriggerOnly/ZS_2023_GK_28thAugust_Realistic"
#   inputPath="/eos/user/p/pjana/LByL/LByL_ZS_Test_2023/ZS_2023_GK_28thAugust_ZSHigh/LbyLSignal_5p02TeV_SuperChic-v4/reco_mc_lbl_try2/230828_155616/0000/mc_HiForestAOD_${1}.root" 
   #inputPath="/eos/user/p/pjana/LByL/LByL_ZS_Test_2023/HIForest/DoSR_LTH8p0_HTH10p0_NoZS_MIEB8p0_MIEE8p0/LbyLSignal_5p02TeV_SuperChic-v4/reco_mc_lbl_try2/230919_170916/0000/mc_HiForestAOD_${1}.root"
   #inputPath="/eos/user/p/pjana/LByL/LByL_ZS_Test_2023/HIForest/DoSR_LTH10p0_HTH14p0_NoZS_MIEB8p0_MIEE8p0/LbyLSignal_5p02TeV_SuperChic-v4/reco_mc_lbl_try2/230920_150633/0000/mc_HiForestAOD_${1}.root"
   #inputPath="/eos/user/p/pjana/LByL/LByL_ZS_Test_2023/HIForest/DoZS_MIEB8p0_MIEE8p0_HIEB2p5_HIEE3p0_CenterTow0/LbyLSignal_5p02TeV_SuperChic-v4/reco_mc_lbl_try2/230920_171135/0000/mc_HiForestAOD_${1}.root"
   inputPath="/eos/user/p/pjana/LByL/LByL_ZS_Test_2023/HIForest/DoSR_LTH4p0_HTH8p0_NoZS_MIEB8p0_MIEE8p0/LbyLSignal_5p02TeV_SuperChic-v4/reco_mc_lbl_try2/230921_054545/0000/mc_HiForestAOD_${1}.root"
  # inputPath=""
 #  outputPath="/eos/cms/store/cmst3/group/lightbylight/Pranati/Final_afterTrigger/LByLMC_SC/"   
outputPath="/eos/user/p/pjana/LByL/LByL_ZS_Test_2023/after_TriggerOnly/GK_From_Chris/DoSR_LTH4p0_HTH8p0_NoZS_MIEB8p0_MIEE8p0"
  # outputPath="/eos/user/p/pjana/LByL/LByL_ZS_Test_2023/after_TriggerOnly/GK_From_Chris/DoZS_MIEB8p0_MIEE8p0_HIEB2p5_HIEE3p0_CenterTow0"
   #outputPath="/eos/user/p/pjana/LByL/LByL_ZS_Test_2023/after_TriggerOnly/GK_From_Chris/DoSR_LTH10p0_HTH14p0_NoZS_MIEB8p0_MIEE8p0"
#    outputPath="/eos/user/p/pjana/LByL/LByL_ZS_Test_2023/after_TriggerOnly/ZS_2023_GK_28thAugust_ZSHigh"
    #outputPath="/eos/user/p/pjana/LByL/LByL_ZS_Test_2023/after_TriggerOnly/Default_LByL2018MC_GK_ForCheck_1"
    #outputPath="/eos/user/p/pjana/LByL/LByL_ZS_Test_2023/after_TriggerOnly/21stAugust_Official_LbyL_Signal"
#   outputPath="/eos/user/p/pjana/LByL/LByL_ZS_Test_2023/after_TriggerOnly/ZS_lblConfig"   
#   outputPath=""
# outputPath="${basePath}/skimmed_ntuples/mc_lbl_sc_doubleEG2_full_lumi"
  elif [ ${3} -eq 4 ] # QED SL, 253 chunks, max chunk number: 253
  then
    echo "CEP"
    sampleName="CEP"
    #inputPath="${basePath}/mc_cep/ntuples_1034/ntuples_cep_1034/ntuples_cep_1034/200211_054704/0000/HiForestAOD_cep_${1}.root"
    #inputPath="${basePath}/mc_cep/ntuples_full_lumi/QCDDiphoton_5p02TeV_SuperChic/reco_mc_cep_full_lumi/200811_121848/0000/mc_HiForestAOD_${1}.root"
    #inputPath="${basePath}/mc_cep/QCDDiphoton_5p02TeV_SuperChic/reco_mc_cep_tryv2/210420_064351/0000/mc_HiForestAOD_${1}.root"
#    inputPath="${basePath}/mcForests/mc_cep/QCDDiphoton_5p02TeV_SuperChic/reco_mc_cep_CastorInfo/210906_103054/0000/mc_HiForestAOD_${1}.root"
 #   outputPath="${basePath}/skimmed_ntuples/mc_cep_sc_doubleEG2_full_lumi"
     #inputPath="/eos/cms/store/group/phys_diffraction/lbyl_2018/mcForests/mc_cepIncoh/gen_sim_cep_SC_incoh/reco_mc_cepIncoh_CastorInfo/211125_053225/0000/mc_HiForestAOD_${1}.root"
     inputPath="/eos/cms/store/group/phys_diffraction/lbyl_2018/mcForests/mc_cepIncoh/gen_sim_cep_SC_incoh/reco_mc_cepIncoh_CastorInfo/211125_053225/0000/CEPIncoh_forest.root "
     outputPath="/eos/cms/store/cmst3/group/lightbylight/Pranati/Final_afterTrigger/CEP_Incoherent/mc_cepIncoh_sc_DoubleEG2"
 
 elif [ ${3} -eq 5] # CEP incoh
  then
    echo "CEPIncoh"
    sampleName="CEPIncoh"
    #inputPath="${basePath}/mc_cep/ntuples_1034/ntuples_cep_1034/ntuples_cep_1034/200211_054704/0000/HiForestAOD_cep_${1}.root"
    #inputPath="${basePath}/mc_cep/ntuples_full_lumi/QCDDiphoton_5p02TeV_SuperChic/reco_mc_cep_full_lumi/200811_121848/0000/mc_HiForestAOD_${1}.root"
    #inputPath="${basePath}/mc_cep/QCDDiphoton_5p02TeV_SuperChic/reco_mc_cep_tryv2/210420_064351/0000/mc_HiForestAOD_${1}.root"
    inputPath="${basePath}/mc_cepIncoh/gen_sim_cep_SC_incoh/reco_mc_cepIncoh_CastorInfo/211125_053225/0000/mc_HiForestAOD_${1}.root"   
    outputPath="/eos/cms/store/cmst3/group/lightbylight/Pranati/Final_afterTrigger/CEP_Incoherent/mc_cepIncoh_sc_DoubleEG2" 
 #  outputPath="/eos/cms/store/group/phys_diffraction/lbyl_2018/skimmed_ntuples/mc_cepIncoh_sc_doubleEG2_full_lumi"
  else
    echo "Unknown option: ${3}"
  fi
  mkdir -p $outputPath
fi

outputReco="${outputPathReco}/ntuples_forRecoEff_${1}.root"
outputTrigger="${outputPathTrigger}/ntuples_forTriggerEff_${1}.root"
outputHFveto="${outputPathHFveto}/ntuples_forHFvetoEff_${1}.root"
outputExclusivity="${outputPathExclusivity}/ntuples_forExclusivityEff_${1}.root"
outputLbLsignal="${outputPathLbLsignal}/ntuples_forLbLsignal_${1}.root"
outputQEDsignal="${outputPathQEDsignal}/ntuples_forQEDsignal_${1}.root"
outputLowAco="${outputPathLowAco}/ntuples_passingLbL_lowAco_${1}.root"
outputHighAco="${outputPathHighAco}/ntuples_passingLbL_highAco_${1}.root"
output="${outputPath}/ntuples_loose_selections_${1}.root"

echo "Config: ${configPath}"
echo "Input: ${inputPath}"
echo "Output reco: ${outputReco}"
echo "Output trigger: ${outputTrigger}"
echo "Output HF veto: ${outputHFveto}"
echo "Output exclusivity: ${outputExclusivity}"
echo "Output LbL signal: ${outputLbLsignal}"
echo "Output QED signal: ${outputQEDsignal}"
echo "Output low aco: ${outputLowAco}"
echo "Output high aco: ${outputHighAco}"
echo "Output: ${output}"

if [ $2 -eq 0 ] # data
then
  ${userBasePath}/LightByLight2018/analysis/applySelections $configPath $inputPath $outputReco $outputTrigger $outputHFveto $outputExclusivity $outputLbLsignal $outputQEDsignal $sampleName
elif [ $2 -eq 1 ] # MC
then
  ${userBasePath}/LightByLight2018/analysis/applySelections $configPath $inputPath $outputReco $outputTrigger $outputHFveto $outputExclusivity $outputLbLsignal $outputQEDsignal $sampleName
elif [ $2 -eq 2 ] # Data passing LbL selections
then
  ${userBasePath}/LightByLight2018/analysis/applySelections $configPath $inputPath $outputLowAco $outputHighAco $sampleName
elif [ $2 -eq 3 ] # Data passing loose selections
then
#  ${userBasePath}/LightByLight2018/analysis/applySelections $configPath $inputPath $output $sampleName
  if [ -s ${output} ]
  then
    echo "File already exists, skipping"
  else
    echo "File doesn't exist or is empty - running"
    ${userBasePath}/LightByLight2018/analysis/applySelections $configPath $inputPath $output $sampleName
  fi

  
  
fi


