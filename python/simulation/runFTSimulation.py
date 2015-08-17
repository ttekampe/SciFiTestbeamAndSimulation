import GaudiPython as GP
from GaudiConf import IOHelper
from Configurables import LHCbApp, ApplicationMgr, DataOnDemandSvc
from Configurables import SimConf, DigiConf, DecodeRawEvent
from Configurables import CondDB

import array

from LinkerInstances.eventassoc import *

import ROOT as R


def resetSipmVals(sipimValPtr):
  for layer in sipimValPtr:
    for adcID in layer:
      for adcChan in layer[adcID]:
        adcChan[0] = 0


LHCbApp().Simulation = True
CondDB().Upgrade = True


LHCbApp().DDDBtag = "dddb-20150424"
LHCbApp().CondDBtag = "sim-20140204-vc-md100"

#work around for bug in DB
CondDB().LoadCALIBDB = 'HLT1'


# Configure all the unpacking, algorithms, tags and input files
appConf = ApplicationMgr()
appConf.ExtSvc+= [
                  'ToolSvc'
                  ,'DataOnDemandSvc'
                  #,'NTupleSvc'
                  ]
appConf.TopAlg += [
                   "MCFTDepositCreator"
                   ,"MCFTDigitCreator"
                   ,"FTClusterCreator"
                   #,"FTNtupleMaker"
                   ]


from Configurables import SiPMResponse
SiPMResponse().useNewResponse = 2

from Configurables import MCFTDepositCreator
MCFTDepositCreator().SpillVector = ["/"]
MCFTDepositCreator().SpillTimes = [0.0]
#MCFTDepositCreator().AttenuationToolName = "Julians nice tool"

from Configurables import MCFTDigitCreator
#MCFTDigitCreator().Irradiation = 0.0 #disables thermal noise
MCFTDigitCreator().Force2bitADC = 0

#could set
#MCFTDigitCreator().SiPMGain

tof = 25.4175840541


MCFTDigitCreator().IntegrationOffset = [26 - tof, 28 - tof, 30 - tof]
MCFTDigitCreator().SimulateNoise = False
#MCFTDigitCreator().Temperature =
MCFTDigitCreator().Irradiation = 0.0
MCFTDigitCreator().SiPMGain = sipm_gain = 1000.
MCFTDigitCreator().CrossTalkProbability = 0.13

f_5l = 0.57595865315812876038
f_6l = 0.63798496965208108843

MCFTDigitCreator().PhotoElectronsPerMeV = 120. * f_5l

#from Configurables import FTClusterCreator
#MCFTClusterCreator().ADCThreshold = #neighbor threshold
#MCFTClusterCreator().ClusterMinADCPeak = #seed threshold
#MCFTClusterCreator().ClusterMinCharge = #sum threshold
#MCFTClusterCreator().ClusterMinWidth
#MCFTClusterCreator().ClusterMaxWidth

### jacos config for craeting nTuple
#from Configurables import GaudiSequencer, FTNtupleMaker, NTupleSvc
#NTupleSvc().Output = ["FILE1 DATAFILE='mytupleFile.root' TYP='ROOT' OPT='NEW'"]

#GaudiSequencer("MoniFTSeq").Members += [FTNtupleMaker()]

#### end

s = SimConf()
SimConf().Detectors = ['VP', 'UT', 'FT', 'Rich1Pmt', 'Rich2Pmt', 'Ecal', 'Hcal', 'Muon']
SimConf().EnableUnpack = True
SimConf().EnablePack = False

d = DigiConf()
DigiConf().Detectors = ['VP', 'UT', 'FT', 'Rich1Pmt', 'Rich2Pmt', 'Ecal', 'Hcal', 'Muon']
DigiConf().EnableUnpack = True
DigiConf().EnablePack = False

dre = DecodeRawEvent()
dre.DataOnDemand = True

lhcbApp = LHCbApp()
lhcbApp.Simulation = True


import sys
inputFiles = [sys.argv[-1]]
IOHelper('ROOT').inputFiles(inputFiles)


# Configuration done, run time!
appMgr = GP.AppMgr()
evt = appMgr.evtsvc()
det = appMgr.detsvc()

#quarter 0 - 3
#sipm id 0 - 15
#sipm chanel 0 - 127


resultPath = "/fhgfs/users/ttekampe/SciFi/testbeamData/simulated/boole"

fileName = sys.argv[1].split("/")[-1][:".sim"] + ".root"

print("Outputfile: " + fileName)

exit()

#outputFile = R.TFile(resultPath + "simulationResponse_fibMatVolCor_newTags_PosA.root", "RECREATE")
outputFile = R.TFile(resultPath + fileName, "RECREATE")
#IOHelper('ROOT').outputFiles(resultPath + "simulationResponse.root")
nLayer = 12
sipmIDs = [1, 2, 3, 4]
sipmValPtr = []
outputTrees = []
outputFile.cd()
for layerNumber in xrange(nLayer):
  outputTrees.append(R.TTree("layer_" + str(layerNumber), "layer_" + str(layerNumber) ) )
  sipmValPtr_thisLayer = {}
  for sipmID in sipmIDs:
    arr = []
    for sipmChan in xrange(128):
      arr.append(array.array("f", [0]))
    sipmValPtr_thisLayer[sipmID] = arr
    for adcChan in xrange(128):
      outputTrees[-1].Branch("Uplink_" + str(sipmID) +"_adc_" + str(adcChan+1), sipmValPtr_thisLayer[sipmID][adcChan] ,"Uplink_" + str(sipmID) +"_adc_" + str(adcChan+1) + "/F")
  sipmValPtr.append(sipmValPtr_thisLayer)


while True:
  appMgr.run(1)
  if not evt['MC/Particles']:
    print "no more particles"
    break

  digits = evt['/Event/MC/FT/Digits'].containedObjects()
  for digit in digits:
    if digit.channelID().sipmId() in sipmIDs and digit.channelID().module() == 1 and digit.channelID().quarter() == 3:
      sipmValPtr[digit.channelID().layer()][digit.channelID().sipmId()][digit.channelID().sipmCell()][0] = digit.adcCount() / sipm_gain

  for t in outputTrees:
    t.Fill()
  resetSipmVals(sipmValPtr)



outputFile.cd()
for t in outputTrees:
  t.Write()
outputFile.Close()
