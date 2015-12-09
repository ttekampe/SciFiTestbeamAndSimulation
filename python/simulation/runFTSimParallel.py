import argparse, subprocess, sys, time

def nFinishedJobs(processes):
    nFinished = 0
    for process in processes:
        #print(process)
        if process.poll() is not None:
            nFinished += 1
    print(str(nFinished) + " jobs are finished so far")
    return nFinished

parser = argparse.ArgumentParser(description='Plot cluster properties from data and simulation.')
parser.add_argument('-f', '--file', type=str, nargs='+')
parser.add_argument('-n', '--nCpu', type=int)
parser.add_argument('-t', '--tag', type=str, default="")

cfg = parser.parse_args()

script = sys.path[0] + '/runFTSimulation.py'
#script = sys.path[0] + '/someTestScript.py'
devnull = open('/dev/null', 'w')

processes = []
try:
    for f in cfg.file:
        while len(processes) >= (cfg.nCpu + nFinishedJobs(processes) ):
            #print("len(processes) is " + str(len(processes)))
            print("Waiting 10 seconds until next test for free cpu")
            time.sleep(10)

        print("Starting process for file " + f)

        processes.append( subprocess.Popen(
        [
            "python", script
            ,'-f' , f
            ,'-t' , cfg.tag
        ]
        ,stdout=devnull
        ,stderr=devnull
        ))


    for process in processes:
        process.wait()

except KeyboardInterrupt:
    for process in processes:
        if not process.poll():
            process.kill()
    exit()

print("finished " + str(len(processes)) + " jobs")
