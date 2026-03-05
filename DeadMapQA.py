import json
import numpy as np
import traceback
from functools import partial
import time
import os
import sys
import uproot

import MakeCanvas
from mylogger import *



tree_input = 'DeadMapTREE.root'
log_file = 'QApy.log'
traceback_file = 'exc.err'

NominalGap = 380*32
UnanchorableThreshold = 330000
TriggerRampSec = 10
zoom_threshold = 0.08

N_CHIPS = 24120
N_CHIPS_IB = 432
N_CHIPS_OB = N_CHIPS - N_CHIPS_IB
N_LANES = 3816
N_LANES_IB = N_CHIPS_IB
N_LANES_OB = N_LANES - N_LANES_IB
N_LANES_ML = 864 # L3,4
N_STAVES = 192
N_STAVES_IB = 48
vNStaves = [12, 16, 20, 24, 30, 42, 48]
vStaveBound = [0, 12, 28, 48, 72, 102, 144, 192]
vLaneBound = [0, 108, 252, 432, 816, 1296, 2472, 3816]
vNLanesPerStave = [9, 9, 9, 16, 16, 28, 28]
vNChipsPerLane = [1, 1, 1, 7, 7, 7, 7]

chipsPerStave = np.array([9 if s < N_STAVES_IB else 112 if s < N_STAVES_IB+24+30 else 196 for s in range(N_STAVES)])

LHCOrbitNS = 88924.6

QAcheck = {}
QAFLAG = 'UNKNOWN' # Updated when printing with the worst score
GLO_RUN = 0
FORUN = -1 # first orbit run
LORUN = -1 # last orbit run

NA = -111

log_verbosity = 2 # 1: no DEBUG, 2: also DEBUG

#_________________________________________________
def LOG(severity, *message):
    logger = Logger(log_file)
    logger.set_highlight_keyword(True)
    logger.set_verbosity(log_verbosity)
    sp = f'[{GLO_RUN}]' if GLO_RUN > 0 else ''
    logger.log(severity,sp,*message)
def FLOG(severity, *message):
    logger = Logger(log_file)
    logger.set_print_on_terminal(False)
    logger.set_verbosity(log_verbosity)
    logger.log(severity,*message)
    
#________________________________________________
def Traceback(severity, *message):
    LOG(severity,*message)
    LOG(severity,f'...full error stack in {traceback_file}')
    logger = Logger(traceback_file)
    logger.set_print_on_terminal(False)
    logger.set_verbosity(9999)
    logger.log(severity,*message)
    with open(traceback_file,'a') as f:
        f.write(traceback.format_exc())
        f.write('\n'+'-'*50+'\n')
        
    
    
#_________________________________________________
#def ChipToLane(chipid):
#    if chipid < N_LANES_IB:
#        return chipid
#    else:
#        return N_LANES_IB + (chipid - N_LANES_IB) // 7

#_________________________________________________
def Mapping(dummy='dummy',chip='na',lane='na'): # use either chip or lane

    if dummy != 'dummy':
        LOG(FATAL,f'Invalid use. dummy={dummy}, chip={chip}, lane={lane}. Exiting to avoid troubles')
        exit()

    if chip != 'na':
        if int(chip) < N_LANES_IB:
            lane_ = int(chip)
        else:
            lane_ = N_LANES_IB + (int(chip) - N_LANES_IB) // 7
    elif lane != 'na':
        lane_ = int(lane)
    else:
        LOG(FATAL,f'Ivalid use of Mapping: dummy={dummy}, chip={chip}, lane={lane}. Exiting')
        exit()

    layer = 0
    laneinlayer = lane_
    for i in range(1,7):
        if lane_ >= vLaneBound[i]:
            layer = i
            laneinlayer = lane_ - vLaneBound[i]

    staveinlayer, laneinstave = divmod(laneinlayer, vNLanesPerStave[layer])
    stave = 0
    for l in range(7):
        stave += int(l<layer)*vNStaves[l] + int(l==layer)*staveinlayer

    #return lane, stave, staveinlayer, layer
    return lane_, stave, layer, staveinlayer, laneinstave
    

#_________________________________________________
def NChipsPerLane(lane):
    if isinstance(lane,np.ndarray):
        return np.where(lane < N_CHIPS_IB, 1, 7)
    else:
        return 1 if lane < N_CHIPS_IB else 7
            
#________________________________________________
def NDead(A,layers='all',element='chip'):
    """
    layers: IB, OB, all, layers (only with element = chip)
    element : chip, lane (=fully dead lanes), lanewchip (=lane with at leas one dead chip)
    """
    if layers == 'layers' and element != 'chip':
        LOG(FATAL,f'Invalid use of NDead. Exiting')
        exit()
        
    if layers == 'IB':
        l1,l2 = (0,N_LANES_IB)
    elif layers == 'OB':
        l1,l2 = (N_LANES_IB, N_LANES)
    elif layers == 'all':
        l1,l2 = (0,N_LANES)

    if element == 'chip':
        if layers == 'layers':
            return [A[vLaneBound[i]:vLaneBound[i+1]].sum() for i in range(7)]
        else:
            return np.sum(A[l1:l2])
    elif element == 'lane':
        return np.sum(A[l1:l2] == NChipsPerLane(np.arange(l1,l2))) 
    elif element == 'lanewchip':
        return np.sum(A[l1:l2] > 0)
        
#______________________________________________
def TimeRollingAverage(x,y,window_size=300):
    rolling_avg = np.convolve(np.array(y), np.ones(window_size)/window_size, mode='valid')
    down_x = np.array(x)[window_size-1::window_size] # Taking every 300th point after the window size
    down_y = rolling_avg[::window_size]  # Take the corresponding rolling averages
    return down_x, down_y

#______________________________________________
def index_clusterizer(steps, full_keys, padding_sec=60):  # passing list interesting steps and full list of orbits. Returning list (a,b) where a and b are the first and last index of each cluster

    if not steps or len(full_keys) < 2:
        LOG(INFO,f'Returning 0 clusters')
        return []

    padding_orb = int(padding_sec / (LHCOrbitNS * 1.e-9))
    LOG(INFO,f"Looking for clusters of high dead fraction. Input: {len(steps)} indices. Merging up to {padding_orb} orbits.")

    unique_steps = sorted(set(steps))
    clusters_nopadding = []
    current_cluster = [unique_steps[0],]

    for s in unique_steps:
        if full_keys[s] - full_keys[current_cluster[-1]] <= 2*padding_orb:
            current_cluster.append(s)
        else:
            clusters_nopadding.append((current_cluster[0], current_cluster[-1]))
            current_cluster = [s,]

    clusters_nopadding.append((current_cluster[0], current_cluster[-1])) # adding last

    clusters = []
    
    # extending each cluster up to padding
    LOG(INFO,f"Extending intervals")
    for A,B in clusters_nopadding:
        a,b = (A,B)
        while True:
            a -= 1
            if a < 1 or (full_keys[A] - full_keys[a]) > padding_orb:
                a += 1
                break
        while True:
            b += 1
            if b >= len(full_keys)-1 or (full_keys[b] - full_keys[B]) > padding_orb:
                b -= 1
                break
        clusters.append((a,b))

    if len(clusters) > 9:
        LOG(ERROR,f"Found {len(clusters)} > 9 clusters. This is not acceptable. Returning no clusters")
        return []

    LOG(INFO,f"Returning {len(clusters)} clusters")
    return clusters
        
 
        
#______________________________________________
def LogQAchecks(checks):
    global QAFLAG
    score = {'GOOD': 1, 'UNKNOWN': 0, 'MEDIUM': -1, 'BAD': -2, 'FATAL': -3}
    worst_score = 999
    if 'Default object' in checks and checks['Default object'] == 'FATAL':
        checks.pop('Map decoded',None)
    for check, val in checks.items():
        LOG(INFO,f'QA CHECK - {check}: {val}')
        if score[val] < worst_score:
            worst_score = score[val]
            QAFLAG = val
    
#_______________________________________________
def process_vector(kv): # function to parallelize over the orbits
    
    k, vec = kv    
    n_dead_l = np.zeros(N_LANES)
    n_dead_s = np.zeros(N_STAVES)
    chip_flg = np.zeros(N_CHIPS) # 1 if dead
    for C in vec:
        c = int(C)
        lan, sta, _, _, _ = Mapping(chip=c)
        n_dead_l[lan] += 1
        n_dead_s[sta] += 1
        chip_flg[c] += 1
    return k, n_dead_l, n_dead_s, chip_flg, np.sum(chip_flg)
    

#________ MAIN ________________________________
def main(doGraphics = True):

    global QAFLAG
    global QAcheck
    global GLO_RUN
    
    LOG(INFO,f'Start. Importing data from {tree_input}')

    now_ = time.time()

    with uproot.open(tree_input) as f:
        t_static = f["t_static"].arrays(library="np")
        staticchipmap1 = list(t_static["static"][0])
        run = int(t_static["run"][0])
        LOG(INFO,f'Run number: {run}')
        GLO_RUN = run
        rctstart = int(t_static["rctstart"][0])
        rctstop  = int(t_static["rctstop"][0])
        orbitreset = int(t_static["orbitreset"][0])
        version = int(t_static["version"][0])
        isdefault = bool(t_static["isdefault"][0])
        fatalcheck = list(t_static["fatal"][0])
        exp_norb = int(t_static["nkeys"][0])

        t_dynamic = f["t_dynamic"].arrays(library="np")
        nwords = list(t_dynamic['nwords'])
        rawkeys = list(t_dynamic['key'])


    if orbitreset > 0 and rctstart > 0 and rctstop > 0:
        firstorbitrun = int((rctstart - orbitreset) / (LHCOrbitNS * 1.e-6))  # millisec / millisec
        lastorbitrun = int((rctstop - orbitreset) / (LHCOrbitNS * 1.e-6))  # millisec / millisec
    else:
        firstorbitrun = lastorbitrun = -1
        
        

        
    lanemap = {}
    stavemap = {} # stavemap[orbit] = number of chips dead in the stave at that step
    ndeadchips = [] # ndeadchips[i] = total number of dead chips at step i --> to be implemented
    critical_steps = [] # list of step indices where either OB or IB dead time is above zoom_threshold
    counter_chip_by_chip = np.zeros(N_CHIPS)
   
    OrbitResetChecked = True

    parallelize = True

    exp_eta = int(exp_norb/1000)

    lanemap_invalidKeys = {} # same as lanemap but only for invalid keys

    if parallelize:
        import concurrent.futures
        LOG(INFO,f'CPU count = {os.cpu_count()}. Expected size {exp_norb}, {exp_eta} seconds to import it.')

        with concurrent.futures.ProcessPoolExecutor() as tor:
            #futures = tor.map(process_vector, raw_data.items())
            futures = tor.map(process_vector, list(zip(rawkeys, t_dynamic["deadchips"])))
            stcc = 0
            for k, v1, v2, v3, n3 in futures:
                
                if stcc % (1 + len(rawkeys) // 4) == 0:
                    LOG(INFO,f'{stcc} / {len(rawkeys)}...')
                    
                if k is None:
                    LOG(ERROR,f'Found None orbit. Skipped')
                    continue
                if orbitreset > 0 and firstorbitrun-UnanchorableThreshold < k < lastorbitrun+UnanchorableThreshold:
                    lanemap[int(k)] = v1
                    stavemap[int(k)] = v2
                    counter_chip_by_chip = counter_chip_by_chip + v3
                elif orbitreset <= 0:
                    lanemap[int(k)] = v1
                    stavemap[int(k)] = v2
                    counter_chip_by_chip = counter_chip_by_chip + v3
                    OrbitResetChecked = False
                else:
                    lanemap_invalidKeys[int(k)] = v1

                stcc += 1
                    

    else: # do not parallize
        exit()
        # MISSING IMPLEMENTATION OF ZERO ORBIT
        for k, values in raw_data.items():
            if isinstance(k,int) or k.isdigit():
                n_dead_l = np.zeros(N_LANES)
                n_dead_s = np.zeros(N_STAVES)
                for c in values:
                    lan, sta, _, _, _ = Mapping(chip=c)
                    n_dead_l[lan] += 1
                    n_dead_s[sta] += 1
                    counter_chip_by_chip[c] += 1
                lanemap[int(k)] = n_dead_l
                stavemap[int(k)] = n_dead_s
                
    
    lanemap = dict(sorted(lanemap.items()))

    # Neet to recompute the keys becuase invalid orbits has been removed
    if set(lanemap.keys()) | set(lanemap_invalidKeys.keys()) == set(rawkeys):
        pass
    else:
        LOG(FATAL,f'Error in building the list of keys. Exiting')
        exit()

    ninvalid = len(lanemap_invalidKeys)
    if OrbitResetChecked:
        LOG(INFO if ninvalid == 0 else WARNING,f'There are {ninvalid} invalid orbits')
    else:
        LOG(ERROR,f'Orbit reset not available. Map range checks will not be effective')
        
    keys = list(lanemap.keys())

    staticchipmap2 = np.where(counter_chip_by_chip == len(keys))[0].tolist() # when OB single chips are saved, this should be equal to statichipmap

    if staticchipmap1 == staticchipmap2:
        LOG(INFO,f'Static maps 1 and 2 are identical')
        staticchipmap = staticchipmap2
    elif len(staticchipmap1) == 0 and len(staticchipmap2) > 0:
        LOG(INFO,f'The static map is empty. Dead chips computed from evolving map')
        staticchipmap = staticchipmap2
    elif len(staticchipmap1) > 0 and len(staticchipmap2) == 0:
        LOG(WARNING,f'There are no chips which are always dead in the evolving map!')
        staticchipmap = staticchipmap1
    else:
        LOG(WARNING,f'Static maps 1 and 2 are both filled but different. Checking if all the dead chips in map 2 belong to a dead lane in the middle of the map.')
        for cc in staticchipmap2:
            lan, _, lay, _, _ = Mapping(chip=cc)
            midkey = keys[len(keys) // 2]
            if lanemap[midkey][lan] != vNChipsPerLane[lay]:
                LOG(FATAL,f'Found event chip {cc} lane {lan} layer {lay}')
                break
        else:
            LOG(INFO,f'Check passed')
        staticchipmap = staticchipmap1
        #exit()

    for check in fatalcheck:
        QAcheck[str(check)] = 'FATAL'
    
    if isdefault:
        QAcheck['Default object'] = 'FATAL'

    if len(keys) == 0:
        if not any(ccc == 'FATAL' for ccc in QAcheck.values()):
            QAcheck['Map size'] = 'FATAL'
        LogQAchecks(QAcheck)
        LOG(WARNING,f'No orbit keys found. Returning FATAL without further actions for this run {run}')
        return 'FATAL' 

    elapsed_time = time.time()-now_
    with open('track_time.dat','a') as tt:
        tt.write(f'{run} {len(keys)} {elapsed_time}\n')

    staticlanemap = np.zeros(N_LANES)
    for c in staticchipmap:
        staticlanemap[Mapping(chip=c)[0]] += 1
 

    minorbit = min(keys)
    maxorbit = max(keys)

    maprange = (maxorbit - minorbit) * LHCOrbitNS * 1.e-9
    rctduration = (rctstop-rctstart) / 1000 if rctstart > 0 else -1

    offsetstart =  minorbit - firstorbitrun
    offsetstartsec = offsetstart * LHCOrbitNS * 1.e-9
    offsetend = maxorbit - lastorbitrun
    offsetendsec = offsetend * LHCOrbitNS * 1.e-9

    LOG(INFO,f'Evolving map size: {len(keys)}')
    LOG(INFO,f'Map range {minorbit} to {maxorbit}, in seconds: {maprange}')
    LOG(INFO,f'Run duration from RCT object, in seconds: {rctduration}')
    LOG(INFO,f'Orbit at run start {hex(firstorbitrun)}')
    LOG(INFO,f'Delta first orbit map-run {offsetstart} = {offsetstartsec} sec')
    LOG(INFO,f'Orbit at run stop {hex(lastorbitrun)}')
    LOG(INFO,f'Delta last orbit map-run {offsetend} = {offsetendsec} sec')


    FullyDeadIB = int(NDead(staticlanemap, 'IB', 'chip'))
    FullyDeadOB = int(NDead(staticlanemap, 'OB', 'chip'))
    LanesWithFullyDeadOB = int(NDead(staticlanemap, 'IB', 'lanewchip'))

    # --- loop over the orbits --------
    gaps = []
    ngap_overnominal = 0
    unAnchorable = 0
    TimeStampFromStart = []
    DeadFractionIB = []
    DeadFractionOB = []
    DeadFractionLay = [[] for _ in range(7)]
    RecoveryRateLay = [[] for _ in range(7)]
  
    
    WorstIBN = -1
    WorstIBStep = 0
    WorstOBN = -1
    WorstOBStep = 0

    LaneDeadTime = np.zeros(N_LANES)
    LaneDeadTimeNoRamp = np.zeros(N_LANES)
    StaveDeadTimeNoRamp = np.zeros(N_STAVES)
    StaveRecoveryPerHour = np.zeros(N_STAVES)
    WorstIBLaneDeadFraction = np.zeros(N_LANES)
    WorstOBLaneDeadFraction = np.zeros(N_LANES)
    LastDeadFraction = np.zeros(N_LANES)

    LanesWithSingleChip = set()

    lane_range = np.arange(N_LANES)

    nRecoIB = nRecoOB = 0

    SecForTriggerRamp = TriggerRampSec
    for orr in keys:
        if (orr - minorbit) * LHCOrbitNS * 1.e-9 > TriggerRampSec:
            SecForTriggerRamp = (orr - minorbit) * LHCOrbitNS * 1.e-9
            LOG(INFO,f'Trigger ramp  duration initially set to {TriggerRampSec} moved to {round(SecForTriggerRamp,2)} sec')
            break
        
    for i in range(len(keys)):

        if i % (1 + len(keys) // 4) == 0:
            LOG(INFO,f'{i} / {len(keys)}...')

        currentorbit = keys[i]
        currentmap = lanemap[currentorbit] # np array with size N_LANES, of number of dead chips per lane

        if i < len(keys)-1:
            gaps.append(keys[i+1] - currentorbit)
            if gaps[-1] > NominalGap:
                ngap_overnominal += 1
            if gaps[-1] > UnanchorableThreshold:
                unAnchorable += (gaps[-1] - UnanchorableThreshold)

        TimeStampFromStart.append( (currentorbit - minorbit) * LHCOrbitNS * 1.e-9)
        deltaTsec = 0 if i >= len(keys)-1 else (keys[i+1] - currentorbit) * LHCOrbitNS * 1.e-9

        #if TimeStampFromStart[-1] > TriggerRampSec and SecForTriggerRamp < 0:
        #    SecForTriggerRamp = TimeStampFromStart[-1]

        # Fraction of dead chips per Barrel
        IBdead = NDead(currentmap,'IB','chip')
        OBdead = NDead(currentmap,'OB','chip')
        DeadFractionIB.append(IBdead / N_CHIPS_IB)
        DeadFractionOB.append(OBdead / N_CHIPS_OB)

        # Fraction of dead chips per Stave
        LayDead = NDead(currentmap,'layers','chip') #LayDead[4] = number of dead chips in L4
                
        for ilay in range(7): 
            DeadFractionLay[ilay].append(LayDead[ilay] / (vNStaves[ilay] * vNLanesPerStave[ilay] * vNChipsPerLane[ilay]))
        

        # Build array of indices with large dead time (>=zoom_threshold)
        if TimeStampFromStart[-1] >= SecForTriggerRamp:
            if DeadFractionIB[-1] > zoom_threshold or DeadFractionOB[-1] > zoom_threshold:
                critical_steps.append(i)
                
                
        # Fill set of lanes with signle chips        
        for lan,n in enumerate(currentmap[N_LANES_IB:], start=N_LANES_IB):
            if 0 < n < NChipsPerLane(lan):
                nsing = n - staticlanemap[lan]
                if nsing > 0:
                    LanesWithSingleChip.add(lan)
                if nsing < 0:
                    LOG(FATAL,'Unexpected number of dead chips in time-evolving vs static map')
                    QAcheck['Other'] = 'FATAL'

        if i < len(keys)-1:
            ddtime = currentmap[lane_range]*(keys[i+1] - currentorbit)/NChipsPerLane(lane_range)
            LaneDeadTime[lane_range] += ddtime # to be normalized by time
            if TimeStampFromStart[-1] >= SecForTriggerRamp and SecForTriggerRamp > 0:
                LaneDeadTimeNoRamp[lane_range] += ddtime

        if IBdead > WorstIBN and TimeStampFromStart[-1] >= SecForTriggerRamp and SecForTriggerRamp >= 0 and currentorbit != maxorbit:
            WorstIBN = IBdead
            WorstIBStep = i 
        
        if OBdead > WorstOBN and TimeStampFromStart[-1] >= SecForTriggerRamp and SecForTriggerRamp >= 0 and currentorbit != maxorbit:
            WorstOBN = OBdead
            WorstOBStep = i

        if currentorbit == maxorbit:
            LastDeadFraction = currentmap[lane_range] / NChipsPerLane(lane_range)

        if i < len(keys)-1 and TimeStampFromStart[-1] >= SecForTriggerRamp and SecForTriggerRamp >= 0:
            isRecoed = (stavemap[keys[i]] == chipsPerStave) & (stavemap[keys[i+1]] < chipsPerStave)
            StaveRecoveryPerHour += isRecoed  # to be normalized by number of hours
            nRecoIB += np.sum(isRecoed[:N_STAVES_IB])
            nRecoOB += np.sum(isRecoed[N_STAVES_IB:])
            for ilay in range(7):
                RecoveryRateLay[ilay].append(np.sum(isRecoed[vStaveBound[ilay]:vStaveBound[ilay+1]]) / deltaTsec)
        else:
            for ilay in range(7):
                RecoveryRateLay[ilay].append(0)  
            
    # -- end loop over orbits

    WorstIBLaneDeadFraction = lanemap[keys[WorstIBStep]][lane_range] / NChipsPerLane(lane_range)
    WorstOBLaneDeadFraction = lanemap[keys[WorstOBStep]][lane_range] / NChipsPerLane(lane_range)


    if len(keys) > 1:
        LaneDeadTime /= (maxorbit-minorbit)
        if TimeStampFromStart[-1] >= SecForTriggerRamp and SecForTriggerRamp >= 0:
            LaneDeadTimeNoRamp /= (maxorbit-minorbit-SecForTriggerRamp*1e9/LHCOrbitNS)
        else:
            LaneDeadTimeNoRamp[:] = NA
        unAnchorableFrac = unAnchorable / (maxorbit - minorbit)
        recoIBperH = nRecoIB / (maprange / 3600)
        recoOBperH = nRecoOB / (maprange / 3600)
        StaveRecoveryPerHour /= (maprange / 3600)
    else:
        LaneDeadTime[:] = NA
        unAnchorableFrac = NA
        LaneDeadTimeNoRamp[:] = NA
        recoIBperH = NA
        recoOBperH = NA
        StaveRecoveryPerHour[:] = NA

    lane_to_stave = np.array([Mapping(lane=l)[1] for l in range(N_LANES)])
    StaveDeadTimeNoRamp = np.bincount(lane_to_stave, weights=LaneDeadTimeNoRamp) / np.bincount(lane_to_stave)

   
    AvgDeadTimeIB = np.mean(LaneDeadTimeNoRamp[:N_LANES_IB])
    AvgDeadTimeOB = np.mean(LaneDeadTimeNoRamp[N_LANES_IB:])
    
    LOG(INFO,f'Lanes with single dead chips: {len(LanesWithSingleChip)}')

    ### LOG for debug and other studies
    # Lanes with single chips
    LOG(DEBUG,f'LWSC run {run} duration {rctduration:.1f} lanes {" ".join(str(x) for x in sorted(LanesWithSingleChip))}')
    if len(LanesWithSingleChip) < 50:
        lnames = ''
        for l_ in sorted(LanesWithSingleChip):
            _, _, la_, st_, ll_ = Mapping(lane=l_)
            lnames += f'L{la_}_{st_}_{ll_} '
        LOG(DEBUG,f'LNWSC {lnames}')
    # Lanes almost completely dead
    lacd = [ilane for ilane in range(N_LANES) if LaneDeadTimeNoRamp[ilane] >= 0.95]
    LOG(DEBUG,f'LACD run {run} duration {rctduration:.1f} lanes {" ".join(str(x) for x in sorted(lacd))}')
    if len(lacd) < 100:
        lnames = ''
        for l_ in sorted(lacd):
            _, _, la_, st_, ll_ = Mapping(lane=l_)
            lnames += f'L{la_}_{st_}_{ll_} '
        LOG(DEBUG,f'LNACD {lnames}')
   
    
    
    wind_size = 1200 if maprange > 12*60*60 else 900 if maprange > 8*60*60 else 600 if maprange > 5*60*60 else 300

    LOG(INFO,f'Using {wind_size} sec as window for the rolling average')

    DeadFrac_rolling_IB_x, DeadFrac_rolling_IB_y = TimeRollingAverage(TimeStampFromStart,DeadFractionIB,window_size=wind_size)
    DeadFrac_rolling_OB_x, DeadFrac_rolling_OB_y = TimeRollingAverage(TimeStampFromStart,DeadFractionOB,window_size=wind_size)
    
    CriticalStepsClusters = index_clusterizer(critical_steps, keys)
    clusterizer_summary = ''
    if len(critical_steps) > 0 and len(CriticalStepsClusters) > 0:
        clusterizer_summary = f'{len(CriticalStepsClusters)} regions with > {100*zoom_threshold}% dead time'
    if len(critical_steps) > 0 and len(CriticalStepsClusters) == 0:
        clusterizer_summary = f'WARN: too many regions with > {100*zoom_threshold}% dead time?'
       
    # Performing the checks
    
    ## Avg dead time IB and OB
    QAcheck['Avg dead time IB'] = 'GOOD' if AvgDeadTimeIB < 0.03 else 'MEDIUM' if AvgDeadTimeIB < 0.1 else 'BAD'
    QAcheck['Avg dead time OB'] = 'GOOD' if AvgDeadTimeOB < 0.05 else 'MEDIUM' if AvgDeadTimeOB < 0.1 else 'BAD'

    ## Single chips
    singfrac = len(LanesWithSingleChip) / N_LANES_OB
    QAcheck['Single chips'] = 'GOOD' if singfrac < 0.02 else 'MEDIUM' if singfrac < 0.05 else 'BAD'

    ## Fully dead IB and OB
    QAcheck['Fully dead IB'] = 'GOOD' if FullyDeadIB < 9 else 'MEDIUM' if FullyDeadIB < 0.1*N_CHIPS_IB else 'BAD'
    QAcheck['Fully dead OB'] = 'GOOD' if LanesWithFullyDeadOB < 68 else 'BAD'

    ## Deafault -> set at the beginning

    ## Map size:
    if isdefault:
        QAcheck['Map size'] = 'GOOD' if len(keys)==0 and len(staticchipmap)==0 else 'FATAL'
    else:
        QAcheck['Map size'] = 'GOOD'
        if len(staticchipmap) == 0:
            QAcheck['Map size'] = 'BAD'
        if len(keys) < 2:
            QAcheck['Map size'] = 'FATAL'

    ## Invalid orbit
    if OrbitResetChecked and ninvalid == 0:
        QAcheck['Invalid orbit'] = 'GOOD'
    elif not OrbitResetChecked:
        QAcheck['Invalid orbit'] = 'UNKNOWN'
    else:
        ncio = [kio for kio, v in lanemap_invalidKeys.items() if NDead(v,'all','chip') < N_CHIPS]
        if ncio:
            QAcheck['Invalid orbit'] = 'MEDIUM'
            if len(ncio) > 10:
                LOG(WARNING,f'More than 10 orbits are invalid with alive chips')
            else:
                for nci in ncio:
                    LOG(WARNING,f'Invalid orbit {hex(nci)} has {NDead(lanemap_invalidKeys[nci],"all","chip")} alive chips')
        else:
            QAcheck['Invalid orbit'] = 'MEDIUM'
    

    ## Orbit gaps
    QAcheck['Orbit gaps'] = 'GOOD'
    n_above = sum([g > NominalGap for g in gaps])
    n_above2 = sum([g > 2*NominalGap for g in gaps])
    if n_above2 > 0 or n_above > 2:
        QAcheck['Orbit gaps'] = 'MEDIUM'
    if unAnchorable > 0 or n_above > max(2, 0.25 * len(keys)):
        QAcheck['Orbit gaps'] = 'BAD'

    ## Orbit range
    if not OrbitResetChecked:
        QAcheck['Orbit range'] = 'UNKNOWN'
    else:
        if offsetstart < -0.001 or offsetend > 0.001:  # zero, modulo rounding errors
            QAcheck['Orbit range'] = 'MEDIUM'
        elif abs(offsetstart) < 3*NominalGap and abs(offsetend) < 3*NominalGap:
            QAcheck['Orbit range'] = 'GOOD'
        elif abs(offsetstart) < UnanchorableThreshold and abs(offsetend) < UnanchorableThreshold:
            QAcheck['Orbit range'] = 'MEDIUM'
        else:
            QAcheck['Orbit range'] = 'BAD'
            

    if QAcheck['Orbit range'] == 'GOOD' and abs(offsetstartsec) > 5:
        QAcheck['Orbit range'] = 'BAD'

    ## Unanchorable fraction
    QAcheck['Un-anchorable fraction'] = 'GOOD' if unAnchorableFrac < 0.02 else 'MEDIUM' if unAnchorableFrac < 0.05 else 'BAD'
    
    LOG(INFO,f'Average IB dead time (no trg ramp): {AvgDeadTimeIB if AvgDeadTimeIB != NA else "n/a"}')
    LOG(INFO,f'Average OB dead time (no trg ramp): {AvgDeadTimeOB if AvgDeadTimeOB != NA else "n/a" }')
    LOG(INFO,f'Stave recoveries (IB/OB): {nRecoIB}/{nRecoOB}')
    LOG(INFO,f'Stave recovery rate (IB/OB) (1/h): {recoIBperH if recoIBperH != NA else "n/a"}/{recoOBperH if recoOBperH != NA else "n/a"}')
    LOG(INFO,f'Nominal gap is {NominalGap}. Number of steps above: {n_above}. Number of steps above 2xnominal: {n_above2}')
    LOG(INFO,f'Unanchorable orbits: {unAnchorable} corresponding to {unAnchorableFrac} of the run duration')

    
    if not np.array_equal(DeadFrac_rolling_IB_x, DeadFrac_rolling_OB_x):
        LOG(WARNING,f'The time steps of the rolling average for IB and OB differ')

       
    Text1 =  f'k#Orbit keys: {len(keys)}'
    if ninvalid > 0:
        Text1 += f' + {ninvalid} neglected'
    Text1 += '#k#' + ','.join(hex(o) for o in keys[:3]) + '...#k#...' + ','.join(hex(o) for o in keys[-3:])
    
    if ninvalid > 0:
        Text1 += f'#r#{ninvalid} keys neglected:'
        if ninvalid < 4:
            Text1 += '#r#' + ','.join(hex(o) for o in lanemap_invalidKeys.keys())
        else:
            Text1 += '#r#too many to print'

    Text1 += f'#w#empty#k#Run start/stop: {hex(firstorbitrun)}, {hex(lastorbitrun)}'
    
    Text1 += f'#w#empty#k#RCT run duration (s): {rctduration:.1f}#k#MAP duration (s): {maprange:.1f}'
    Text1 += f'#k#Offsets (s): {round(offsetstartsec,3)}, {round(offsetendsec,3)}'
    Text1 += f'#w#empty#k#Dead chips (IB+OB): {FullyDeadIB} + {FullyDeadOB}'
    Text1 += f'#k#OB lanes w/ single dead chips: {len(LanesWithSingleChip)}'
    if clusterizer_summary:
        Text1 += f'#w#empty#k#{clusterizer_summary}'

    Text2 = f'b#Run {run}#w#empty'
    colorcode = {'GOOD': 'g', 'MEDIUM': 'orange', 'BAD': 'r', 'FATAL': 'm'}
    for check, val in QAcheck.items():
        try:
            col = colorcode[val]
        except:
            col = 'k'
        Text2 += f'#{col}#{check}: {val}'

    # Logging the QA checks. This will also set the global quality returned by main()
    LogQAchecks(QAcheck)

    if doGraphics:
        LOG(INFO,'Passing results to graphic functions')

        try:
            MakeCanvas.make_canvas1(
                lane_dead_time = LaneDeadTimeNoRamp.tolist(),
                stave_dead_time = StaveDeadTimeNoRamp.tolist(),
                number_of_fully_dead = staticlanemap.tolist(),
                stave_recovery_rate = StaveRecoveryPerHour.tolist(),
                gaps = gaps,
                dead_fraction = [{'both':list(range(len(keys)))}, {'IB':DeadFractionIB, 'OB':DeadFractionOB}],
                dead_fraction_rolling = [{'IB':(DeadFrac_rolling_IB_x/60).tolist(), 'OB':(DeadFrac_rolling_IB_x/60).tolist()}, {'IB':DeadFrac_rolling_IB_y.tolist(), 'OB':DeadFrac_rolling_OB_y.tolist()}],
                words = [[],nwords], # first is dummy for the number of dead chips step by step... to be implemented
                WorstOBstep = WorstOBStep,
                WorstIBstep = WorstIBStep,
                worst_ob = WorstOBLaneDeadFraction.tolist(),
                worst_ib = WorstIBLaneDeadFraction.tolist(),
                #last = LastDeadFraction.tolist(),
                last = [int(i in LanesWithSingleChip) for i in range(N_LANES)],
                text1 = Text1,
                text2 = Text2
            )
        except Exception as e:
            Traceback(ERROR,f'Exception canvas 1: {e}')

        if QAcheck['Map size'] == 'FATAL':
            LOG(WARNING,f'Map contains {len(keys)} key. Returning FATAL without further actions for this run {run}')
            return 'FATAL'    
    
        lanemap_fraction = {}
        for k,vec in lanemap.items():
            ff = np.zeros(N_LANES)
            ff[lane_range] = vec[lane_range] / NChipsPerLane(lane_range)
            lanemap_fraction[k] = ff
            
        try:
            MakeCanvas.make_canvas2(
                lanemap = lanemap_fraction,
                idx = [[0,N_LANES_IB], [N_LANES_IB, N_LANES_IB+N_LANES_ML], [N_LANES_IB+N_LANES_ML,N_LANES]],
                offset_sec = 0,
                spec1= [f' - run {run}',]*3, 
                run=run
                )
        except Exception as e:
            Traceback(ERROR,f'Exception canvas 2: {e}')
    
        try:
            MakeCanvas.make_canvas4(
                lane_dead_time = LaneDeadTimeNoRamp.tolist(),
                stave_dead_time = StaveDeadTimeNoRamp.tolist(),
                run=run
                )
        except Exception as e:
            Traceback(ERROR,f'Exception canvas 4: {e}')

        try:
            MakeCanvas.make_canvas5(
                dead0 = list(TimeRollingAverage(TimeStampFromStart,DeadFractionLay[0],window_size=wind_size)),
                dead1 = list(TimeRollingAverage(TimeStampFromStart,DeadFractionLay[1],window_size=wind_size)),
                dead2 = list(TimeRollingAverage(TimeStampFromStart,DeadFractionLay[2],window_size=wind_size)),
                dead3 = list(TimeRollingAverage(TimeStampFromStart,DeadFractionLay[3],window_size=wind_size)),
                dead4 = list(TimeRollingAverage(TimeStampFromStart,DeadFractionLay[4],window_size=wind_size)),
                dead5 = list(TimeRollingAverage(TimeStampFromStart,DeadFractionLay[5],window_size=wind_size)),
                dead6 = list(TimeRollingAverage(TimeStampFromStart,DeadFractionLay[6],window_size=wind_size)),
                reco0 = list(TimeRollingAverage(TimeStampFromStart,RecoveryRateLay[0],window_size=wind_size)),
                reco1 = list(TimeRollingAverage(TimeStampFromStart,RecoveryRateLay[1],window_size=wind_size)),
                reco2 = list(TimeRollingAverage(TimeStampFromStart,RecoveryRateLay[2],window_size=wind_size)),
                reco3 = list(TimeRollingAverage(TimeStampFromStart,RecoveryRateLay[3],window_size=wind_size)),
                reco4 = list(TimeRollingAverage(TimeStampFromStart,RecoveryRateLay[4],window_size=wind_size)),
                reco5 = list(TimeRollingAverage(TimeStampFromStart,RecoveryRateLay[5],window_size=wind_size)),
                reco6 = list(TimeRollingAverage(TimeStampFromStart,RecoveryRateLay[6],window_size=wind_size)),
                run = str(GLO_RUN)
                )
        except Exception as e:
            Traceback(ERROR,f'Exception canvas 5: {e}')
               

        # Making several canvas2
        zoom_index = 0
        for A,B in CriticalStepsClusters: # A,B are the first and last index of the zoom window
            
            zoom_index += 1

            try:
                lanemap_fraction_zoom = {keys[i]: lanemap_fraction[keys[i]] for i in range(A,B+1)}
                maxIB = 100*max(DeadFractionIB[i] for i in range(A,B+1) if TimeStampFromStart[i] > SecForTriggerRamp)
                maxOB = 100*max(DeadFractionOB[i] for i in range(A,B+1) if TimeStampFromStart[i] > SecForTriggerRamp)
                iat = [i for i in range(A,B+1) if max(DeadFractionIB[i], DeadFractionOB[i]) > zoom_threshold and TimeStampFromStart[i] > SecForTriggerRamp]
                first_orb = keys[iat[0]]
                last_orb = keys[iat[-1]]
                first_sec = TimeStampFromStart[iat[0]]
                last_sec = TimeStampFromStart[iat[-1]]
                firstIB = 100*DeadFractionIB[iat[0]]
                firstOB = 100*DeadFractionOB[iat[0]]
                lastIB = 100*DeadFractionIB[iat[-1]]
                lastOB = 100*DeadFractionOB[iat[-1]]
                maxIB = 100*max(DeadFractionIB[i] for i in iat)
                maxOB = 100*max(DeadFractionOB[i] for i in iat)           

                c2spec0 = f', with > {100*zoom_threshold}% missing. Zoom #{zoom_index}: orb {hex(first_orb)} = {first_sec:.1f}s ({firstIB:.1f}; {firstOB:.1f})% to orb {hex(last_orb)} = {last_sec:.1f}s ({lastIB:.1f}; {lastOB:.1f})%. Max dead fraction ({maxIB:.1f}; {maxOB:.1f})%'
                LOG(INFO,f'Zoom details {c2spec0}')
                c2spec2 = f'zoom{zoom_index}'
            
                MakeCanvas.make_canvas2(
                    lanemap = lanemap_fraction_zoom,
                    idx = [[0,N_LANES_IB], [N_LANES_IB, N_LANES_IB+N_LANES_ML], [N_LANES_IB+N_LANES_ML,N_LANES]],
                    offset_sec = TimeStampFromStart[A],
                    run = run,
                    spec0 = c2spec0,
                    spec1= [f' - run {run}',]*3, 
                    spec2= c2spec2
                    )
            except Exception as e:
                Traceback(ERROR,f'Exception while creating zoom canvas n. {zoom_index}: {e}')

        
                    

    LOG(INFO,f'Dumping statistics on ITSstat.json')
    # stave dead time
    j_stave_dead_time = dict(enumerate(StaveDeadTimeNoRamp))
    # problematic lanes
    problanes = [ilane for ilane in range(N_LANES) if LaneDeadTimeNoRamp[ilane] >= 0.70]
    if len(problanes) < 50:
        j_problematic_lanes = {}
        for l_ in problanes:
            _, _, la_, st_, ll_ = Mapping(lane=l_)
            j_problematic_lanes[f'L{la_}_{st_}_{ll_}'] = LaneDeadTimeNoRamp[l_]
    else:
        j_problematic_lanes = "too_many"
    # disappearing single chips
    if len(LanesWithSingleChip) < 50:
        j_lanes_single_chips = ''
        for l_ in sorted(LanesWithSingleChip):
            _, _, la_, st_, ll_ = Mapping(lane=l_)
            j_lanes_single_chips += f'L{la_}_{st_}_{ll_},'
        j_lanes_single_chips = j_lanes_single_chips[:-1]
    else:
        j_lanes_single_chips = "too_many"
        
            
     
    JJ = {
        'run': GLO_RUN,
        'rct_duration': rctduration,
        'map_flag': QAFLAG,
        'stave_dead_time': j_stave_dead_time,
        'lanes_dead_time_070cut': j_problematic_lanes,
        'lanes_with_disappearing_single_chips': j_lanes_single_chips
    }

    with open(f"ITSstat.json", "w") as j_f:
        json.dump(JJ, j_f, indent=4)
                
    LOG(INFO,f'Returning worst quality {QAFLAG}')
    return QAFLAG

###########  
if __name__ == "__main__":

    Usage = f"""
       {sys.argv[0]} [no-graphics]
       or
       {sys.argv[0]} input.root [no-graphics]
       default input is {tree_input}
    """

    if '-h' in sys.argv or '--help' in sys.argv:
        print(Usage)
        exit()

    if len(sys.argv) > 1:
        if '.json' in sys.argv[1]:
            tree_input = str(sys.argv[1])

    LOG(INFO,f'Running QA on file {tree_input}')
    nographics = 'no-graphics' in sys.argv
    main(not nographics)
    with open("QAHANDSHAKE", "w") as f:
        f.write("ok")
    exit()
