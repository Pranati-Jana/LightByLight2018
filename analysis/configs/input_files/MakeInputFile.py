# MC data at /eos/cms/store/group/phys_top/gkrintir/LbL/GammaGamma2MuMu_GENSIM_woPtCut_PbPb5TeV_v1/crab_GammaGamma2MuMu_GENSIM_woPtCut/200317_235208/0000/
# Real Data at/eos/cms/store/group/phys_diffraction/lbyl_2018/HIForward_Reco/ntuples/ntuples_data/HIForward/ntuples_data_lbl/200617_140125/0000/data_HiForwardAOD_1-1.root

import sys
import csv
import os.path

python_config_file = sys.argv[1]
input_file = sys.argv[2]
batch_number = sys.argv[3]
job_number = sys.argv[4]

ending =  batch_number + "_" + job_number
runfile = "/afs/cern.ch/user/p/pjana/private/CMSSW_10_3_5/src/LightByLight2018/analysis/configs/input_files/condor_inputfiles/tmp_runfile_" + ending + ".txt"


with open(python_config_file) as f:
    python_config = f.readlines()
	
	
config = python_config[0].rstrip()
SigEle = python_config[1].rstrip()
MuEle = python_config[2].rstrip()
MuMu = python_config[3].rstrip()

headers = ["Batch_Number", "Job_Number","config_file", "input_file", "Output_file_1", "Output_file_2", "Output_file_3"]
dic = {}

dic["Batch_Number"] = batch_number
dic["Job_Number"] = job_number
dic["input_file"] = input_file
dic["config_file"] = config
dic["Output_file_1"] = SigEle +  ".root"
dic["Output_file_2"] = MuEle + ending + ".root"
dic["Output_file_3"] = MuMu + ending + ".root"

f = open(runfile,"w")
f.write( dic["config_file"]  + "\n")
f.write( dic["input_file"]  + "\n")
f.write( dic["Output_file_1"]  + "\n")
f.write( dic["Output_file_2"]  + "\n")
f.write( dic["Output_file_3"]  + "\n")
f.close()
