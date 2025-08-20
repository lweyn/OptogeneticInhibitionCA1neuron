import numpy as np
import scipy.signal as sig

#______________________________________HELP FUNCTIONS____________________________________________________
def getSpiketimes(t,v,threshold=-10):
    t=np.array(t)
    v = np.array(v)
    idx_spikes = sig.find_peaks(np.array(v), height = threshold)[0]
    spiketimes = t[idx_spikes]
    spikes = v[idx_spikes]
    return spikes,spiketimes, idx_spikes

def getFiringRate(spiketimes,t_start,t_end, t_transient = 0):
    t_start = t_start+t_transient
    return np.where((spiketimes > t_start) & (spiketimes < t_end))[0].size/(t_end-t_start)*1000

def CalctInhibitMax(spiketimes,t_start,t_end):
    if len(np.where((spiketimes > t_start) & (spiketimes < t_end))[0]) >1:
        return np.max(np.diff(spiketimes[np.where((spiketimes > t_start) & (spiketimes < t_end))[0]]))
    else: return None 

def CalctInhibit(spiketimes,t_optostart):
    if len(np.where((spiketimes > t_optostart))[0]) > 1:
        return spiketimes[np.where((spiketimes > t_optostart))[0]][0]-t_optostart
    else: return None 

def calcAPheight(spikes, spiketimes, t_start, t_end):
    id = np.where((spiketimes>t_start) & (spiketimes<t_end))[0]
    if len(id>0): return np.mean(np.array(spikes[id]))
    else: return 0

def calcAPbase(t, v, spiketimes, idx_spikes, t_start, t_end):
    spikes_id = idx_spikes[np.array(np.where((spiketimes>t_start) & (spiketimes<t_end))[0])]
    not_id=[np.arange(spikes_id[i]-30, spikes_id[i]+30) for i in range(len(spikes_id))]
    id = np.array(np.where((t>t_start) & (t<t_end))[0])
    new_id = []
    for i in id:
        if i not in np.array(not_id): new_id.append(i)
    return np.mean(v[new_id])

def calcAPpeaktopeak(v, spiketimes, idx_spikes, t_start, t_end):
    id = idx_spikes[np.array(np.where((spiketimes>t_start) & (spiketimes<t_end))[0])]
    id_min = [np.argmin(v[id[i]:id[i+1]])+id[i] for i in range(0, len(id)-1)]
    if len(id>0): return np.mean(v[id[:-1]]-v[id_min])
    else: return 0 

def ionMov(c):
    dc = np.diff(c)
    iPos = dc[np.where(dc>=0)[0]]
    iNeg = dc[np.where(dc< 0)[0]]
    return (np.sum(iPos), np.sum(iNeg))

def getStartEndTimes(input, name):
    if name == 'init': return input.stimopt.Estimparams.delay, input.stimopt.Ostimparams.delay
    elif name == 'end': return input.stimopt.Ostimparams.delay + input.stimopt.Ostimparams.dur, input.stimopt.Estimparams.delay + input.stimopt.Estimparams.dur
    elif name == 'opt': return input.stimopt.Ostimparams.delay, input.stimopt.Ostimparams.delay + input.stimopt.Ostimparams.dur
    else: assert("feature error: unknown time window")

#______________________________________SPECIFIC FEATURES____________________________________________________

def OstimPower(input, results, loc):
    return input.stimopt.Ostimparams.dc*input.stimopt.Ostimparams.amp

def FR(input, results, loc):
    FRs = {}
    for name in ['init', 'end', 'opt']:
        t_start, t_end = getStartEndTimes(input, name)
        FRs[name] = {'all': getFiringRate(results.APs[loc][1], t_start, t_end),
                    'tr': getFiringRate(results.APs[loc][1], t_start, t_start+input.analysesopt.t_tr),
                    'st': getFiringRate(results.APs[loc][1], t_start + input.analysesopt.t_tr, t_end)}
    return FRs

def APheight(input, results, loc):
    AP_heights = {}
    spikes, spiketimes, _ = results.APs[loc]
    for name in ['init', 'end', 'opt']:
        t_start, t_end = getStartEndTimes(input, name)
        AP_heights[name] = {'all': calcAPheight(spikes, spiketimes, t_start, t_end),
                            'tr': calcAPheight(spikes, spiketimes, t_start, t_start+input.analysesopt.t_tr),
                            'st': calcAPheight(spikes, spiketimes, t_start+input.analysesopt.t_tr, t_end)}
    return AP_heights

def APbase(input, results, loc):
    AP_bases = {}
    _, spiketimes, idx_spikes = results.APs[loc]
    for name in ['init', 'end', 'opt']:
        t_start, t_end = getStartEndTimes(input, name)
        AP_bases[name] = {'all': calcAPbase(results.t, results.recordings.v[loc], spiketimes, idx_spikes, t_start, t_end),
                          'tr': calcAPbase(results.t, results.recordings.v[loc], spiketimes, idx_spikes, t_start, t_start+input.analysesopt.t_tr),
                          'st': calcAPbase(results.t, results.recordings.v[loc], spiketimes, idx_spikes, t_start+input.analysesopt.t_tr, t_end)}
    return AP_bases

def APpeaktopeak(input, results, loc):
    AP_p2ps = {}
    _, spiketimes, idx_spikes = results.APs[loc]
    for name in ['init', 'end', 'opt']:
        t_start, t_end = getStartEndTimes(input, name)
        AP_p2ps[name] = {'all': calcAPpeaktopeak(results.recordings.v[loc], spiketimes, idx_spikes, t_start, t_end),
                          'tr': calcAPpeaktopeak(results.recordings.v[loc], spiketimes, idx_spikes, t_start, t_start+input.analysesopt.t_tr),
                          'st': calcAPpeaktopeak(results.recordings.v[loc], spiketimes, idx_spikes, t_start+input.analysesopt.t_tr, t_end)}
    return AP_p2ps

def tInhibit(input, results, loc):
    return {'all': CalctInhibit(results.APs[loc][1], input.stimopt.Ostimparams.delay),
            'tr': CalctInhibit(results.APs[loc][1], input.stimopt.Ostimparams.delay+input.analysesopt.t_tr)}

def iOpto(input, results, loc):
    Opto = results.recordings.iopto[loc]
    iPos = Opto[np.where(Opto>=0)[0]]
    iNeg = Opto[np.where(Opto<0)[0]]
    return {'end': Opto[np.where(results.t<input.stimopt.Ostimparams.dur+input.stimopt.Ostimparams.delay)[0][-1]],
            'min': np.min(Opto),
            'max': np.max(Opto),
            'totP': np.sum(iPos),
            'totN': np.sum(iNeg)}

def gOpto(input, results, loc):
    Opto = results.recordings.gopto[loc]
    iPos = Opto[np.where(Opto>=0)[0]]
    iNeg = Opto[np.where(Opto<0)[0]]
    return {'end': Opto[np.where(results.t<input.stimopt.Ostimparams.dur+input.stimopt.Ostimparams.delay)[0][-1]],
            'min': np.min(Opto),
            'max': np.max(Opto),
            'totP': np.sum(iPos),
            'totN': np.sum(iNeg)}

def clErev(input, results, loc):
    t_start_opt, t_stop_opt = getStartEndTimes(input, 'opt')
    id_opt = np.where((results.t>t_start_opt) & (results.t<t_stop_opt))[0]
    return {'initEstim':results.recordings.ecl[loc][np.where(results.t<input.stimopt.Estimparams.delay)[0][-1]],
            'initOstim': results.recordings.ecl[loc][np.where(results.t<input.stimopt.Ostimparams.delay)[0][-1]],
            'endEstim': results.recordings.ecl[loc][np.where(results.t<input.stimopt.Estimparams.delay+input.stimopt.Estimparams.dur)[0][-1]],
            'endOstim': results.recordings.ecl[loc][np.where(results.t<input.stimopt.Ostimparams.delay+input.stimopt.Ostimparams.dur)[0][-1]],
            'optMax': np.max(results.recordings.ecl[loc][id_opt]),
            'optMin': np.min(results.recordings.ecl[loc][id_opt])}

def kErev(input, results, loc):
    t_start_opt, t_stop_opt = getStartEndTimes(input, 'opt')
    id_opt = np.where((results.t>t_start_opt) & (results.t<t_stop_opt))[0]
    return {'initEstim':results.recordings.ek[loc][np.where(results.t<input.stimopt.Estimparams.delay)[0][-1]],
            'initOstim': results.recordings.ek[loc][np.where(results.t<input.stimopt.Ostimparams.delay)[0][-1]],
            'endEstim': results.recordings.ek[loc][np.where(results.t<input.stimopt.Estimparams.delay+input.stimopt.Estimparams.dur)[0][-1]],
            'endOstim': results.recordings.ek[loc][np.where(results.t<input.stimopt.Ostimparams.delay+input.stimopt.Ostimparams.dur)[0][-1]],
            'optMax': np.max(results.recordings.ek[loc][id_opt]),
            'optMin': np.min(results.recordings.ek[loc][id_opt])}

def clo(input, results, loc):
    t_start_opt, t_stop_opt = getStartEndTimes(input, 'opt')
    id_opt = np.where((results.t>t_start_opt) & (results.t<t_stop_opt))[0]
    return {'initEstim':results.recordings.clo[loc][np.where(results.t<input.stimopt.Estimparams.delay)[0][-1]],
            'initOstim': results.recordings.clo[loc][np.where(results.t<input.stimopt.Ostimparams.delay)[0][-1]],
            'endEstim': results.recordings.clo[loc][np.where(results.t<input.stimopt.Estimparams.delay+input.stimopt.Estimparams.dur)[0][-1]],
            'endOstim': results.recordings.clo[loc][np.where(results.t<input.stimopt.Ostimparams.delay+input.stimopt.Ostimparams.dur)[0][-1]],
            'optMax': np.max(results.recordings.clo[loc][id_opt]),
            'optMin': np.min(results.recordings.clo[loc][id_opt])}

def ko(input, results, loc):
    t_start_opt, t_stop_opt = getStartEndTimes(input, 'opt')
    id_opt = np.where((results.t>t_start_opt) & (results.t<t_stop_opt))[0]
    return {'initEstim':results.recordings.ko[loc][np.where(results.t<input.stimopt.Estimparams.delay)[0][-1]],
            'initOstim': results.recordings.ko[loc][np.where(results.t<input.stimopt.Ostimparams.delay)[0][-1]],
            'endEstim': results.recordings.ko[loc][np.where(results.t<input.stimopt.Estimparams.delay+input.stimopt.Estimparams.dur)[0][-1]],
            'endOstim': results.recordings.ko[loc][np.where(results.t<input.stimopt.Ostimparams.delay+input.stimopt.Ostimparams.dur)[0][-1]],
            'optMax': np.max(results.recordings.ko[loc][id_opt]),
            'optMin': np.min(results.recordings.ko[loc][id_opt])}

def cli(input, results, loc):
    t_start_opt, t_stop_opt = getStartEndTimes(input, 'opt')
    id_opt = np.where((results.t>t_start_opt) & (results.t<t_stop_opt))[0]
    return {'initEstim':results.recordings.cli[loc][np.where(results.t<input.stimopt.Estimparams.delay)[0][-1]],
            'initOstim': results.recordings.cli[loc][np.where(results.t<input.stimopt.Ostimparams.delay)[0][-1]],
            'endEstim': results.recordings.cli[loc][np.where(results.t<input.stimopt.Estimparams.delay+input.stimopt.Estimparams.dur)[0][-1]],
            'endOstim': results.recordings.cli[loc][np.where(results.t<input.stimopt.Ostimparams.delay+input.stimopt.Ostimparams.dur)[0][-1]],
            'optMax': np.max(results.recordings.cli[loc][id_opt]),
            'optMin': np.min(results.recordings.cli[loc][id_opt])}

def ki(input, results, loc):
    t_start_opt, t_stop_opt = getStartEndTimes(input, 'opt')
    id_opt = np.where((results.t>t_start_opt) & (results.t<t_stop_opt))[0]
    return {'initEstim':results.recordings.ki[loc][np.where(results.t<input.stimopt.Estimparams.delay)[0][-1]],
            'initOstim': results.recordings.ki[loc][np.where(results.t<input.stimopt.Ostimparams.delay)[0][-1]],
            'endEstim': results.recordings.ki[loc][np.where(results.t<input.stimopt.Estimparams.delay+input.stimopt.Estimparams.dur)[0][-1]],
            'endOstim': results.recordings.ki[loc][np.where(results.t<input.stimopt.Ostimparams.delay+input.stimopt.Ostimparams.dur)[0][-1]],
            'optMax': np.max(results.recordings.ki[loc][id_opt]),
            'optMin': np.min(results.recordings.ki[loc][id_opt])}

def cliMov(input, results, loc):
    c = results.recordings.cli[loc]
    cMovDict = {'tot': ionMov(c[100:])}
    for name in ['init', 'end', 'opt']:
        t_start, t_stop = getStartEndTimes(input, name)
        cMovDict[name] = ionMov(c[np.where((results.t>t_start) & (results.t<t_stop))[0]])
    return cMovDict

def cloMov(input, results, loc):
    c = results.recordings.clo[loc]
    cMovDict = {'tot': ionMov(c[100:])}
    for name in ['init', 'end', 'opt']:
        t_start, t_stop = getStartEndTimes(input, name)
        cMovDict[name] = ionMov(c[np.where((results.t>t_start) & (results.t<t_stop))[0]])
    return cMovDict

def koMov(input, results, loc):
    c = results.recordings.ko[loc]
    cMovDict = {'tot': ionMov(c[100:])}
    for name in ['init', 'end', 'opt']:
        t_start, t_stop = getStartEndTimes(input, name)
        cMovDict[name] = ionMov(c[np.where((results.t>t_start) & (results.t<t_stop))[0]])
    return cMovDict

def kiMov(input, results, loc):
    c = results.recordings.ki[loc]
    cMovDict = {'tot': ionMov(c[100:])}
    for name in ['init', 'end', 'opt']:
        t_start, t_stop = getStartEndTimes(input, name)
        cMovDict[name] = ionMov(c[np.where((results.t>t_start) & (results.t<t_stop))[0]])
    return cMovDict
