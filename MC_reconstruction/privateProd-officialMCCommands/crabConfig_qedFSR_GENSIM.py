from CRABClient.UserUtilities import config
config = config()

#config.section_('General')
config.General.requestName = 'gen_sim_qed_FSR_starlightv2'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

#config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '1ststep_GEN_SIM.py'
#config.JobType.maxMemoryMB = 4000
config.Data.outputPrimaryDataset = 'qed_FSR_starlight'
config.Data.userInputFiles = open('qed_fsr_SL.txt').readlines()
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#NJOBS = 2000  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
#config.Data.totalUnits = config.Data.unitsPerJob * NJOBS

#config.section_('Data')
config.Data.outLFNDirBase = '/store/group/phys_diffraction/lbyl_2018/mc_qed/gen_sim_FSR_starlight'
config.Data.allowNonValidInputDataset = True
config.Data.publication = True
config.Data.outputDatasetTag = 'qed_FSR_starlight'
config.Site.storageSite = 'T2_CH_CERN'