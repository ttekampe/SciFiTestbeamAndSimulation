from Ganga.GPI import *
import sys
import inspect
import os


# Stereo angle to use
stereo = 5 #int(sys.argv[1])

local_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

pos = "c"

print local_dir

j = Job(application=Gauss(
        #version="v46r7p2",
        version="v48r2",
        optsfile=local_dir + "/gauss-pgun-job.py",
        extraopts="execute(\"%s\", %s)"%(pos, stereo)
        #user_release_area=local_dir + "/../cmt",
        ))

j.name = "SciFi-stereo-%i"%stereo
j.outputfiles = [LocalFile("*.sim")]

j.backend=Dirac()
#j.backend.settings['BannedSites'] = ['LCG.RRCKI.ru']


events = 10000
eventsperjob = 500


j.splitter = GaussSplitter(numberOfJobs=int(round(events*1.00000001/eventsperjob)),
                           eventsPerJob=eventsperjob)


j.postprocessors.append(LHCbFileMerger(files = ['testbeam_simulation_position_' + pos + '.sim','stdout'],ignorefailed=True,overwrite=True))


j.prepare()
j.submit()
