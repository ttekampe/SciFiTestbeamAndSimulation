import argparse, subprocess, sys


parser = argparse.ArgumentParser(description='Plot cluster properties from data and simulation.')
parser.add_argument('-f', '--file', type=str, nargs='+')
parser.add_argument('-n', '--nCpu', type=int)

cfg = parser.parse_args()


processes = []
try:
    for f in cfg.file:
        if len(processes) < cfg.nCpu:
            print("Starting process for file " + f)
            #commands = '''source /home/ttekampe/.bashrc
#python '''
            #commands += sys.path[0] + "/runFTSimulation.py -f" + f
            #process = subprocess.Popen('/bin/bash', stdin=subprocess.PIPE, stdout=subprocess.PIPE)
            #out, err = process.communicate(commands)
            #processes.append((process, out, err))
            script = sys.path[0] + '/runFTSimulation.py'
            processes.append( subprocess.Popen([
                "python", script
                ,'-f' , f
            ]))
        else:
            for process in processes:
                process[0].wait()

except KeyboardInterrupt:
    for process in processes:
        process[0].kill()
    exit()
