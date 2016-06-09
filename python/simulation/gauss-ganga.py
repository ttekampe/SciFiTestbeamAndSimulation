from Ganga.GPI import *
import sys
import inspect
import os



#attScanPositions = ["305"
#                    ,"355"
#                    ,"455"
#                    ,"555"
#                    ,"655"
#                    ,"755"
#                    ,"855"
#                    ,"955"
#                    ,"1055"
#                    ,"1155"
#                    ,"1255"
#                    ,"1355"
#                    ,"1455"
#                    ,"1555"
#                    ,"1655"
#                    ,"1755"
#                    ,"1855"
#                    ,"1955"
#                    ,"2055"
#                    ,"2155"
#                    ,"2255"
#                    ]

attScanPositions = ["a", "c"]

# angle between beam and fibre mat in x
angles = [0, 10, 20, 30] #int(sys.argv[1])

local_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

#pos = "c"

print local_dir

for pos in attScanPositions:
	for angle in angles:
		j = Job(application=Gauss(
		        #version="v46r7p2",
		        version="v48r2",
		        #version="v48r3",
		        optsfile=local_dir + "/gauss-job-V5_renato.py",
		        extraopts="execute(\"%s\", %s)"%(pos, angle),
		        user_release_area="/home/ttekampe/cmtuser"
		        #user_release_area=local_dir + "/../cmt",
		        ))
		j.name = "SF-Pos_" + pos + "_at_" + str(angle)
		j.outputfiles = [LocalFile("*.sim")]

		#j.inputfiles=[LocalFile("/home/ttekampe/SciFi/FTv5/DDDB_FTv5_20150424_s20140204_lhcbv38r6.db")]
		j.inputfiles =      ['/home/ttekampe/SciFi/FTv5/DDDB_FTv5_20150424_s20140204_lhcbv38r6.db']

		j.backend=Dirac()
		#j.backend=Local()
		#j.backend.settings['BannedSites'] = ['LCG.RRCKI.ru']


		    #events = 10000
		    #eventsperjob = 500

		events = 10000
		eventsperjob = 1000
		j.splitter = GaussSplitter(numberOfJobs=int(round(events*1.00000001/eventsperjob)),
		                           eventsPerJob=eventsperjob)
		j.postprocessors.append(LHCbFileMerger(files = ['testbeam_simulation_position_' + pos + '_at_' + str(angle) + 'deg.sim'],overwrite=True))
		j.prepare()
		j.submit()
