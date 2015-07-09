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

stereo = 5

LHCbApp().Simulation = True
CondDB().Upgrade = True
t = {
    "DDDB": "dddb-20131025",
    "CondDB": "sim-20130830-vc-md100",
    "Others": ["VP_UVP+RICH_2019+UT_UUT",
               "FT_StereoAngle%s"%stereo,
               "Muon_NoM1", "Calo_NoSPDPRS"],
    }
LHCbApp().DDDBtag = t['DDDB']
LHCbApp().CondDBtag = t['CondDB']
CondDB().AllLocalTagsByDataType = t['Others']

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
#appConf.addTool("SiPMResponse")


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


#MCFTDigitCreator().IntegrationOffset = [0.6,2.6,4.6]
MCFTDigitCreator().IntegrationOffset = [26 - tof, 28 - tof, 30 - tof]
MCFTDigitCreator().SimulateNoise = False
#MCFTDigitCreator().Temperature =
MCFTDigitCreator().Irradiation = 0.0
#MCFTDigitCreator().
MCFTDigitCreator().SiPMGain = sipm_gain = 1000.

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


resultPath = "/net/storage03/data/users/ttekampe/SciFi/testbeamSimu/PosA/"


outputFile = R.TFile(resultPath + "simulationResponse_fibMatVolCor_PosA.root", "RECREATE")
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



n=0
#with open("testOutput.txt", "w") as f:
while True:
  appMgr.run(1)
  if not evt['MC/Particles']:
    print "no more particles"
    break

  digits = evt['/Event/MC/FT/Digits'].containedObjects()
  for digit in digits:
    if digit.channelID().sipmId() in sipmIDs and digit.channelID().module() == 1 and digit.channelID().quarter() == 3:
      sipmValPtr[digit.channelID().layer()][digit.channelID().sipmId()][digit.channelID().sipmCell()][0] = digit.adcCount() / sipm_gain
  #  f.write("layer: " + str(digit.channelID().layer()) + "\tmodule: " +str(digit.channelID().module()) +"\tquarter: " + str(digit.channelID().quarter()) + "\tsipmID: " + str(digit.channelID().sipmId()) + "\tsipm channel: " + str(digit.channelID().sipmCell()) + "\n")
  #f.write("\n")

  #hits = evt['/Event/MC/FT/Hits']#.containedObjects()
  #print("found " + str(len(hits)) + " hits")
  #for x in xrange(len(hits)):
  #  hit = hits.containedObject(x)
  #  f.write("particleID: " + str(hit.mcParticle().particleID().pid()) + "\norigin vertex: (" + str(hit.mcParticle().originVertex().position().X()) +
  #  ", " + str(hit.mcParticle().originVertex().position().Y()) + " ," + str(hit.mcParticle().originVertex().position().Z()) + ")\n")
  #  if(hit.mcParticle().mother()):
  #    f.write("mother: " + str(hit.mcParticle().mother().particleID().pid()) + "\n")
  #    if(hit.mcParticle().mother().mother()):
  #      f.write("mothers mother: " + str(hit.mcParticle().mother().mother().particleID().pid()) + "\n")
  #  f.write("hit at (" + str(hit.midPoint().X()) + ", " + str(hit.midPoint().Y()) + ", " + str(hit.midPoint().Z()) + ")\n" )
  #  f.write("started at (" + str(hit.entry().X()) + ", " + str(hit.entry().Y()) + ", " + str(hit.entry().Z()) + ")\n" )
  #  f.write("left at (" + str(hit.exit().X()) + ", " + str(hit.exit().Y()) + ", " + str(hit.exit().Z()) + ")\n" )
  #  f.write("\n")
  #f.write("\n\n")

  for t in outputTrees:
    t.Fill()
  resetSipmVals(sipmValPtr)


  #n+=1
  #if n>50: break

outputFile.cd()
for t in outputTrees:
  t.Write()
outputFile.Close()
