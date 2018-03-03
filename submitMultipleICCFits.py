#!/usr/bin/python
import sys
import numpy as np
import os
import socket

hostName = socket.gethostname()
if 'sns' in hostName:
    print 'Running on an SNS machine!'
    mantidPythonPath="/SNS/users/ntv/mantid/mantid/release/bin"
elif 'pcwf' in hostName: 
    mantidPythonPath="/home/ntv/mantid/mantid/bin"
else:
    mantidPythonPath=None
    
#numRuns = 181; numJobs = 25
#numRuns=11; numJobs=6
numRuns=9; numJobs=3

jobsList = np.array_split(range(numRuns), numJobs)

for jobs in jobsList:
    command = 'nohup python doICCFit.py -r '+ str(jobs)[1:-1]+ ' -p '+ mantidPythonPath + ' >/dev/null 2>&1  &'
    print command
    #os.system(command)
