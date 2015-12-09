# SetupProject gauss v46r7p2
#

import sys
import inspect
import os
from math import tan, radians

from Gauss.Configuration import *
#from Gaudi.Configuration import *
from Configurables import Gauss, LHCbApp
import GaudiKernel.SystemOfUnits as units

local_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
sys.path.append(local_dir)

#from common import set_tags

from Configurables import LHCbApp, CondDB

def execute(pos="c", angle=0):
  importOptions("$APPCONFIGOPTS/Gauss/Beam7000GeV-md100-nu7.6-HorExtAngle.py")

  importOptions("$LBPYTHIA8ROOT/options/Pythia8.py")
  importOptions("$APPCONFIGOPTS/Gauss/G4PL_FTFP_BERT_EmNoCuts.py")
  importOptions("$APPCONFIGOPTS/Conditions/Upgrade.py")
  importOptions("$APPCONFIGOPTS/Persistency/Compression-ZLIB-1.py")

  #importOptions("$APPCONFIGOPTS/Gauss/Gauss-Upgrade-Baseline-20131029.py")
  # FTv5
  importOptions('$APPCONFIGOPTS/Gauss/Gauss-Upgrade-Baseline-20150522.py')
  

  outpath = "testbeam_simulation_position_" + pos  + '_at_' + str(angle) + 'deg'

  Gauss().DataType = "Upgrade"

  #LHCbApp().DDDBtag = "dddb-20150424"
  #LHCbApp().CondDBtag = "sim-20140204-vc-md100"

  #LHCbApp().DDDBtag = "dddb-20150424"
  #LHCbApp().CondDBtag = "sim-20140204-vc-md100"

  # FTv5 from Luigi
  LHCbApp().DDDBtag = "dddb-20150424"
  LHCbApp().CondDBtag = "sim-20140204-vc-md100"
  #DDDBConf().DbRoot = "/home/ttekampe/SciFi/FTv5/DDDB_FTv5_20150424_s20140204_lhcbv38r6/lhcb.xml"

  #work around for bug in DB
  CondDB().LoadCALIBDB = 'HLT1'
  CondDB().addLayer(dbFile = "DDDB_FTv5_20150424_s20140204_lhcbv38r6.db", dbName="DDDB" )

  importOptions('$LBPGUNSROOT/options/PGuns.py')
  from Configurables import ParticleGun
  #ParticleGun().EventType = 52210010

  # Set momentum
  from Configurables import MaterialEval
  ParticleGun().addTool(MaterialEval, name="MaterialEval")
  ParticleGun().ParticleGunTool = "MaterialEval"

# test beam position jargon
#position a: 225.5 cm (near mirror) ~5 cm distance from mirror
#position b: 125.5 cm
#position c: 30.5 cm (near sipm) ~ 5 cm distance from sipm
#default y table position: 72.4 cm


  moduleWidth = 552.4 + 3 # 3 = modul gap
  z_orig = 7834. # 7620
  z_target = 9439.
  x_orig = 4. * moduleWidth + 65.3 # centre of the innermost fibre mat of the second module from left when looking into beam direction (neglected half a gap)
  #y_orig = 2417.5
  if pos == "a":
    y_orig = 50 # 5 cm from mirror
  elif pos == "c":
    y_orig = 2417.5 - 50. # 5 cm from SiPM
  elif pos.isdigit():
    y_orig = float(pos)
  else:
    exit()

  ParticleGun().MaterialEval.Xorig = x_orig
  ParticleGun().MaterialEval.Yorig = y_orig
  #ParticleGun().MaterialEval.Zorig = 7620
  ParticleGun().MaterialEval.Zorig = z_orig
  ParticleGun().MaterialEval.ModP = 150000 #150GeV

  ParticleGun().MaterialEval.ZPlane = z_target
  ParticleGun().MaterialEval.Xmin = x_orig - 1.7 + (z_target - z_orig) / tan( radians(90 - angle) )
  ParticleGun().MaterialEval.Xmax = x_orig + 1.7 + (z_target - z_orig) / tan( radians(90 - angle) )
  ParticleGun().MaterialEval.Ymin = y_orig - 1.7
  ParticleGun().MaterialEval.Ymax = y_orig + 1.7
  ParticleGun().MaterialEval.PdgCode = 211

  # Set min and max number of particles to produce in an event
  from Configurables import FlatNParticles
  ParticleGun().addTool(FlatNParticles, name="FlatNParticles")
  ParticleGun().NumberOfParticlesTool = "FlatNParticles"
  ParticleGun().FlatNParticles.MinNParticles = 1
  ParticleGun().FlatNParticles.MaxNParticles = 1

  GaussGen = GenInit("GaussGen")
  GaussGen.FirstEventNumber = 1
  GaussGen.RunNumber = 1082

  LHCbApp().EvtMax = 10

  HistogramPersistencySvc().OutputFile = outpath+'-GaussHistos.root'

  OutputStream("GaussTape").Output = "DATAFILE='PFN:%s.sim' TYP='POOL_ROOTTREE' OPT='RECREATE'"%outpath

#execute()
