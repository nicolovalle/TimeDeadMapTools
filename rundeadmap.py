#!/usr/bin/env python

import sys
import re
import os
import time
import inspect
from datetime import datetime
import matplotlib.pyplot as plt
import statistics
import requests
import json
import subprocess

INFO = 'INFO'
WARNING = 'WARNING'
ERROR = 'ERROR'
FATAL = 'FATAL'
DEBUG = 'DEBUG'

# Usage: ./rundeadmap.py <run_number>

logfile = 'log.log'

#_________________________________________________________________________________
def LOG(severity, *message):
    global logfile
    filename = str(inspect.stack()[1][1])
    filename = filename.split('/')[-1]
    funcname = str(inspect.stack()[1][3])

    tt = datetime.now()
    tstamp = "%s-%s-%s-%s:%s:%s"%(str(tt.year)[-2:],str(tt.month).zfill(2),str(tt.day).zfill(2),str(tt.hour).zfill(2),str(tt.minute).zfill(2),str(tt.second).zfill(2))
    print("[%s][%s][%s]"%(tstamp,severity,filename+':'+funcname),*message)

    writestring = ''
    for s in message:
        writestring = writestring + ' ' + str(s)
    with open(logfile,'a') as f:
        f.write("[%s][%s][%s] "%(tstamp,severity,filename+':'+funcname)+writestring+'\n')

#________________________________________________________________________________
def Exit(severitylog = INFO):
    print(severitylog,'EXITING')
    exit()

#________________________________________________________________________________
def execute(command, log = True, severitylog = INFO):
    if log:
        LOG(severitylog,"Exec -- ",command)
    os.system(command)

#________________________________________________________________________________
def querylogbook(run):
    LOG(INFO,'Querying bookkeeping for run',run)
    req = requests.get('https://ali-bookkeeping.cern.ch/api/runs?filter[runNumbers]=%s&page[offset]=0'%(str(run)),verify=False)
    data = json.loads(req.text)['data']
    if len(data) != 1:
        LOG(ERROR,'Data from bookeeping of wrong size:',len(data))
        Exit(FATAL)
    else:
        for run in data:
            period = run['lhcPeriod']
            list_det = run['detectors'].split(',')
        IsITS = 'ITS' in list_det
        IsMFT = 'MFT' in list_det
        secDuration = int(run['runDuration'])/1000       
        LOG(INFO,'Returning period',period,' duration',secDuration,' Detectors: ITS',IsITS,', MFT',IsMFT)
        return str(period), secDuration, IsITS, IsMFT
        
    
    

#_______________________________ MAIN ___________________________________________

if __name__ == "__main__":

    LOG(INFO,"Starting ",sys.argv[0],sys.argv[1])

    logfile = str(sys.argv[1])+'.log'
    
    with open(logfile,'w') as f:
        f.write(logfile+'\n')
    LOG(INFO,"Starting ",sys.argv[0],sys.argv[1])
        
    try:
        run = int(sys.argv[1])
    except:
        LOG(FATAL, "Inputs not recognized.")
        Exit()
        
    period, secDuration, doITS, doMFT = querylogbook(run)

    year = '20'+period[3:5]
    run = str(run)
    LOG(INFO,"Processing run",run,"from period",period,"year",year,'. Duration(min) = ',secDuration/60)

    if os.path.exists('./output/'+run):
        LOG(INFO,"Directory ./output/"+run,"already exists. Deleting it.")
        execute("rm -fr ./output/"+run, False)
         

    targetdir = './output/'+run
    execute('mkdir -p '+targetdir)

    execute('alien_find /alice/data/'+year+'/'+period+'/'+run+'/raw/ o2_ctf* > '+targetdir+'/full_ctf_list.dat')

    with open(targetdir+'/full_ctf_list.dat') as f:
        tflist = f.readlines()

    if len(tflist) < 3:
        LOG(FATAL,"CTF list is too short. Run will not be processed. Exiting")
        Exit(FATAL)

    LOG(INFO,"Found:",len(tflist),"CTFs")


    try:
        epnlist = [ re.search('epn[0-9]{3}', l).group(0) for l in tflist]
    except:
        LOG(FATAL,"Failing in extracting EPN list")
        Exit(FATAL)

    targetEPN = max(set(epnlist), key=epnlist.count)

    LOG(INFO,"Choosing epn",targetEPN,"producing",epnlist.count(targetEPN),"CTFs, and adding first and last file.")

    ctflist = targetdir+'/alien_ctf_'+targetEPN+'.dat'
    execute('head -n 1 '+targetdir+'/full_ctf_list.dat | grep -v '+targetEPN+' | sed "s_/alice/_alien:///alice/_" > '+ctflist)
    execute('grep '+targetEPN+' '+targetdir+'/full_ctf_list.dat | sed "s_/alice/_alien:///alice/_" >> '+ctflist)
    execute('tail -n 1 '+targetdir+'/full_ctf_list.dat | grep -v '+targetEPN+' | sed "s_/alice/_alien:///alice/_" >> '+ctflist)

    with open(ctflist, 'r') as ctfile:
        line_count = sum(1 for line in ctfile)
        LOG(INFO,'Number of files going to be processed:',line_count)

    wflog = targetdir+'/o2-deadmapbuilder.log'
    wferr = targetdir+'/o2-deadmapbuilder.err'
    LOG(INFO,"Executing workflow. std out, std err:",wflog,",",wferr)

    tflen = 128 if int(year)<2023 else 32

    LOG(INFO,'Using TF length =',tflen,'orbits')
    
    ITSMFTcommand = 'o2-ctf-reader-workflow -b --ctf-input '+ctflist+' --remote-regex "^alien:///alice/data/.+" --copy-cmd no-copy --onlyDet ITS,MFT --shm-segment-size 40000000000 | o2-itsmft-deadmap-builder-workflow --local-output --output-dir '+targetdir+' --source clusters --tf-sampling 1 --tf-length '+str(tflen)+' --shm-segment-size 4000000000 -b | o2-itsmft-deadmap-builder-workflow --runmft --local-output --output-dir '+targetdir+' --source clusters --skip-static-map --tf-sampling 1 --tf-length '+str(tflen)+' --shm-segment-size 4000000000 -b --run > '+wflog+' 2>'+wferr

    ITScommand = 'o2-ctf-reader-workflow -b --ctf-input '+ctflist+' --remote-regex "^alien:///alice/data/.+" --copy-cmd no-copy --onlyDet ITS --shm-segment-size 40000000000 | o2-itsmft-deadmap-builder-workflow --local-output --output-dir '+targetdir+' --source clusters --tf-sampling 1 --tf-length '+str(tflen)+' --shm-segment-size 4000000000 -b --run > '+wflog+' 2>'+wferr
 
    t__start = time.time()
    if doITS and doMFT:
        execute(ITSMFTcommand)
    elif doITS:
        exectute(ITScommand)
    else:
        LOG(FATAL,'ITS is not in the run, doing nothing')
        Exit()
    t__stop = time.time()

    LOG(INFO,"Process took",round(t__stop - t__start,2),"sec.",round((t__stop-t__start)/60.,2),"min.")

    LOG(INFO,"Checking orbit uniformity")

    with open(wferr) as f:
        errlines = f.readlines()
        if len(errlines) > 0:
            LOG(ERROR,'There are errors in the std err')
            
    with open(wflog) as f:
        loglines = f.readlines()

    orbits = []
    for ll in loglines:
        if 'TF received. First orbit' in ll and 'deadmap-builder_its' in ll:
            orbits.append(int(re.search('First orbit [0-9]+',ll).group(0).replace('First orbit','')))

    if len(orbits) == 0:
        LOG(FATAL,"No good workflow output, exiting")
        Exit(FATAL)

    for ll in loglines:
        if 'ERROR' in ll or 'Error' in ll:
            LOG(ERROR,'There are ERRORS in wf stdout')
            break

    orbits.sort()

    orbgap = [orbits[i+1] - orbits[i] for i in range(len(orbits)-1)]
    step = [i+1 for i in range(len(orbits)-1)]

    plt.plot(orbgap,'o-r')
    plt.savefig(targetdir+'/orbits.png')

    OrbitRangeSec = (orbits[-1]-orbits[0])*89.e-6
    LOG(INFO,"Number of TFs:",len(orbits)," Orbit range",orbits[0],":",orbits[-1]," corresponding to ",OrbitRangeSec/60,"min")
    if abs(OrbitRangeSec - secDuration)/secDuration > 0.05:
        LOG(ERROR,'Big difference between range in the map and run duration')
    sevcheck = INFO
    sigmaorb = statistics.pstdev(orbgap)
    if max(orbgap) > 320000:
        sevcheck = ERROR
    elif max(orbgap) > 32000:
        sevcheck = WARNING
    LOG(sevcheck,"Max orbit gap:",max(orbgap),"Min orbit gap:",min(orbgap),"Std orbit gap:",sigmaorb)
    LOG(INFO,targetdir+'/orbits.png created.')
    LOG(INFO,'Run',run,'completed. Running QA on object. Setting timout = 60 sec')

    try:
        
        execute('mkdir '+targetdir+'/ITSQA/')
        rootcommand = ['root', '-b', 'DeadMapQA.C("'+targetdir+'/its_time_deadmap.root",'+str(run)+',"'+targetdir+'/ITSQA/")']
        process = subprocess.Popen(rootcommand, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        stdout, stderr = process.communicate(timeout=60)

        with open(targetdir+'/ITSQA/root.log','w') as fqa:
            fqa.write('========\n stderr \n========')
            fqa.write(stderr.decode())
            fqa.write('========\n stdout \n========')
            fqa.write(stdout.decode())

        with open(targetdir+'/ITSQA/DeadMapQA.log') as fqa:
            loglines = fqa.readlines()
            for ll in loglines:
                if 'ERROR' in ll:
                    LOG(ERROR,'There are errors in the object QA')
                    break
                if 'WARNING' in ll:
                    LOG(WARNING,'There are warnings in the object QA')
                   
    except Exception as e:

        LOG(ERROR,'QA produced exceptions:',e)  

    LOG(INFO,'rundeadmap.py reached the end.')
    Exit()