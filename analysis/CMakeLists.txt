cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project (LightByLight)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

find_package(ROOT REQUIRED COMPONENTS Minuit Eve Geom)
include(${ROOT_USE_FILE})

SET(CMAKE_SKIP_BUILD_RPATH FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE) 
SET(CMAKE_INSTALL_RPATH "/cvmfs/sft.cern.ch/lcg/releases/vdt/0.4.3-992df/x86_64-centos7-gcc8-opt/lib:/cvmfs/sft.cern.ch/lcg/releases/tbb/2019_U7-ba3eb/x86_64-centos7-gcc8-opt/lib:/cvmfs/sft.cern.ch/lcg/releases/Davix/0.7.3-d94fa/x86_64-centos7-gcc8-opt/lib64/:/cvmfs/sft.cern.ch/lcg/releases/png/1.6.37-9c2fe/x86_64-centos7-gcc8-opt/lib/:/cvmfs/sft.cern.ch/lcg/releases/GSL/2.5-32fc5/x86_64-centos7-gcc8-opt/lib/")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


add_executable(selectEleMuEvents  selectEleMuEvents.cpp)
#add_executable(selectGenEleMuEvents  selectGenEleMuEvents.cpp)
#add_executable(EleOneProng  EleOneProng.cpp)
#add_executable(OnlyTriggerPass  OnlyTriggerPass.cpp)
#add_executable(prepareCutflow  prepareCutflow.cpp)
#add_executable(selectElectrons  selectElectrons.cpp)
#add_executable(TunfoldTree  TunfoldTree.cpp)
add_executable(selectMuMuEvents  selectMuMuEvents.cpp)


set_property(TARGET selectEleMuEvents   PROPERTY CXX_STANDARD 17)
#set_property(TARGET selectGenEleMuEvents   PROPERTY CXX_STANDARD 17)
#set_property(TARGET selectElectrons   PROPERTY CXX_STANDARD 17)
#set_property(TARGET TunfoldTree   PROPERTY CXX_STANDARD 17)
set_property(TARGET selectMuMuEvents   PROPERTY CXX_STANDARD 17)

include_directories("${PROJECT_SOURCE_DIR}/include")

add_library(Analysis  src/Event.cpp
                      src/EventProcessor.cpp
                      src/PhysObject.cpp
                      src/PhysObjectProcessor.cpp
                      src/ConfigManager.cpp
                      src/DumpEvent.cpp
)

set_property(TARGET Analysis PROPERTY CXX_STANDARD 17)


#target_link_libraries(prepareCutflow Analysis ${ROOT_LIBRARIES})
target_link_libraries(selectEleMuEvents Analysis ${ROOT_LIBRARIES})
#target_link_libraries(selectGenEleMuEvents Analysis ${ROOT_LIBRARIES})
#target_link_libraries(selectElectrons Analysis ${ROOT_LIBRARIES})
#target_link_libraries(TunfoldTree Analysis ${ROOT_LIBRARIES})
target_link_libraries(selectMuMuEvents Analysis ${ROOT_LIBRARIES})
