#!/bin/bash

alias cmake='/cvmfs/sft.cern.ch/lcg/contrib/CMake/3.13.4/Linux-x86_64/bin/cmake'


. /cvmfs/sft.cern.ch/lcg/releases/LCG_104/ROOT/6.28.04/x86_64-el9-gcc11-opt/bin/thisroot.sh
#. /cvmfs/sft.cern.ch/lcg/releases/LCG_99/ROOT/v6.22.06/x86_64-centos7-gcc10-opt/bin/thisroot.sh




export lcgenv=/cvmfs/sft.cern.ch/lcg/releases/LCG_105/lcgenv/1.3.22/x86_64-el9-gcc12-opt/lcgenv

. /cvmfs/sft.cern.ch/lcg/releases/LCG_105/vdt/0.4.4/x86_64-el9-gcc11-opt/vdt-env.sh

. /cvmfs/sft.cern.ch/lcg/contrib/gcc/11/x86_64-el9/setup.sh

. /cvmfs/sft.cern.ch/lcg/releases/LCG_105/tbb/2021.10.0/x86_64-el9-gcc11-opt/tbb-env.sh

. /cvmfs/sft.cern.ch/lcg/releases/LCG_105/Davix/0.8.4/x86_64-el9-gcc11-opt/Davix-env.sh

export LD_LIBRARY_PATH=/eos/cms/store/group/phys_diffraction/lbyl_2018/libPNG/libpng-1.6.37/install/lib/:$LD_LIBRARY_PATH

. /cvmfs/sft.cern.ch/lcg/releases/LCG_105/GSL/2.7/x86_64-el9-gcc11-opt/GSL-env.sh

. /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc11-opt/setup.sh

#. /cvmfs/sft.cern.ch/lcg/releases/LCG_104/ROOT/6.28.04/x86_64-el9-gcc13-opt/bin/thisroot.sh

