from Ganga.GPI import *
import sys
import inspect
import os



attScanPositions = ["305"
                    ,"355"
                    ,"455"
                    ,"555"
                    ,"655"
                    ,"755"
                    ,"855"
                    ,"955"
                    ,"1055"
                    ,"1155"
                    ,"1255"
                    ,"1355"
                    ,"1455"
                    ,"1555"
                    ,"1655"
                    ,"1755"
                    ,"1855"
                    ,"1955"
                    ,"2055"
                    ,"2155"
                    ,"2255"
                    ]

# angle between beam and fibre mat in x
angle = 0 #int(sys.argv[1])

local_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

#pos = "a"

print local_dir

for pos in attScanPositions:
    j = Job(application=Gauss(
            #version="v46r7p2",
            version="v48r2",
            optsfile=local_dir + "/gauss-pgun-job.py",
            extraopts="execute(\"%s\", %s)"%(pos, angle)
            #user_release_area=local_dir + "/../cmt",
            ))

    j.name = "SciFi-Pos" + pos + "_angle_" + str(angle)
    j.outputfiles = [LocalFile("*.sim")]

    j.backend=Dirac()
    #j.backend=Local()
    #j.backend.settings['BannedSites'] = ['LCG.RRCKI.ru']


    #events = 10000
    #eventsperjob = 500

    events = 10000
    eventsperjob = 1000

    j.splitter = GaussSplitter(numberOfJobs=int(round(events*1.00000001/eventsperjob)),
                               eventsPerJob=eventsperjob)


    j.postprocessors.append(LHCbFileMerger(files = ['testbeam_simulation_position_' + pos + '_at_' + str(angle) + 'deg.sim'],ignorefailed=True,overwrite=True))


    j.prepare()
    j.submit()
