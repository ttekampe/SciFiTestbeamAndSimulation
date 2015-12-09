from Gaudi.Configuration import importOptions
from Configurables import LHCbApp
from Gauss.Configuration import *

import sys
import inspect
import os
from math import tan, radians

local_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
sys.path.append(local_dir)

#information on colliding beams: momenta, crossing angles
#importOptions('$APPCONFIGOPTS/Gauss/Beam7000GeV-md100-nu7.6-HorExtAngle.py')

#event type
#importOptions('$DECFILESROOT/options/13104012.py')    #Bs phi phi


def execute(pos="c", angle=0):

    #MC generator to use - all upgrade studies are done with PYTHIA 8
    importOptions('$LBPYTHIA8ROOT/options/Pythia8.py')

    #to enable hadronic physics in GEANT4
    importOptions('$APPCONFIGOPTS/Gauss/G4PL_FTFP_BERT_EmNoCuts.py')

    #Upgrade conditions
    importOptions('$APPCONFIGOPTS/Conditions/Upgrade.py')

    importOptions("$APPCONFIGOPTS/Gauss/Beam7000GeV-md100-nu7.6-HorExtAngle.py")
    importOptions("$APPCONFIGOPTS/Persistency/Compression-ZLIB-1.py")


    ################
    #user options
    sub_id = 0

    #number of events to generate
    LHCbApp().EvtMax = 10000

    #geometry options
    GeoV5 = True

    ################

    #to enable spillover and 25 ns bunch crossing
    #EnableSpillover = True

    #if EnableSpillover :
    #    importOptions('$APPCONFIGOPTS/Gauss/EnableSpillover-25ns.py')
    #else :
    #    #without spillover, I manually specify the bunch crossing
    #    from GaudiKernel import SystemOfUnits
    #    GenInit("GaussGen").BunchSpacing = 25 * SystemOfUnits.ns

    ################


    #generation seeds are controlled by the event number and the run number
    GaussGen = GenInit("GaussGen")
    GaussGen.FirstEventNumber = (sub_id * LHCbApp().EvtMax) + 1
    GaussGen.RunNumber        = (sub_id * 1000) 

        
    #detector conditions

    from Configurables import DDDBConf

    if GeoV5 :
        #customized geometry
        from Configurables import Gauss, CondDB
        
     
        # V5 geo based on dddb-20150424, June 2015
        importOptions('$APPCONFIGOPTS/Gauss/Gauss-Upgrade-Baseline-20150522.py')
        
        #from Configurables import Gauss, CondDB
        #CondDB().Upgrade = True
        #Gauss().DetectorGeo  = { "Detectors": ['VP', 'UT', 'FT', 'Rich1Pmt', 'Rich2Pmt', 'Ecal', 'Hcal', 'Muon', 'Magnet' ] }
        #Gauss().DetectorSim  = { "Detectors": ['VP', 'UT', 'FT', 'Rich1Pmt', 'Rich2Pmt', 'Ecal', 'Hcal', 'Muon', 'Magnet' ] }
        #Gauss().DetectorMoni = { "Detectors": ['VP', 'UT', 'FT', 'Rich1Pmt', 'Rich2Pmt', 'Ecal', 'Hcal', 'Muon', 'Magnet' ] }

        #Gauss().DetectorGeo  = { "Detectors": ['FT'] }
        #Gauss().DetectorSim  = { "Detectors": ['FT'] }
        #Gauss().DetectorMoni = { "Detectors": ['FT'] }

        #Gauss().DataType = "Upgrade"


        LHCbApp().DDDBtag    = "dddb-20150424"
        LHCbApp().CondDBtag  = "sim-20140204-vc-md100"
        
        #xml files
        #local interactive DB root file
        #DDDBConf().DbRoot = "/afs/cern.ch/user/d/delbuono/cmtuser/DDDB_FTv5_20150424_s20140204_lhcbv38r6/lhcb.xml"
        #ganga (sandbox) non local db file
        CondDB().addLayer(dbFile = "DDDB_FTv5_20150424_s20140204_lhcbv38r6.db", dbName="DDDB" ) #if loaded in the ganga script it will be automatically loaded by Gauss
           
    #########################################################################

    # This is usually not needed, but right now there is a bug
    # which tries to search caliboff.db and fails
    from Configurables import CondDB
    CondDB().LoadCALIBDB = 'HLT1'

    #########################################################################
    #Options = "Gauss_"


    #Options = "BsPhiPhi"+"_V5_"+"Spillover"
    #OutputStream("GaussTape").Output = "DATAFILE='PFN:"+Options+".sim' TYP='POOL_ROOTTREE' OPT='RECREATE'"
    outpath = "testbeam_simulation_position_" + pos  + '_at_' + str(angle) + 'deg'
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


    #moduleWidth = 552.4 + 3 # 3 = modul gap
    moduleWidth = 529.0 + 3 # 3 = modul gap
    z_orig = 7820. # 7620
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


    LHCbApp().EvtMax = 10

    HistogramPersistencySvc().OutputFile = outpath+'-GaussHistos.root'

    OutputStream("GaussTape").Output = "DATAFILE='PFN:%s.sim' TYP='POOL_ROOTTREE' OPT='RECREATE'"%outpath
