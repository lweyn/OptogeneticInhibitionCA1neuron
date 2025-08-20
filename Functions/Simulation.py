from neuron import h, load_mechanisms, rxd

from Functions.utils import Dict, MyEncoder

from Functions.CellSetup import setupCell

import Functions.features as Features
import matplotlib.pyplot as plt
import numpy as np

import json
from datetime import datetime
from time import time

from Functions.CellSetup import setupCell

h.load_file('stdrun.hoc') #initialize NEURON
load_mechanisms('./Model/Mods_Opto/', warn_if_already_loaded=True)
load_mechanisms('./Model/Mods_Gentiletti/', warn_if_already_loaded=True)

def runSimulation(input):

    #to detect a problem with deletion of neurons
    for sec in h.allsec():
        print(sec)

    #____________________________________________SETUP_____________________________________________ 
    StartTime = time()

    #set integration time step of neuron simulator
    h.dt = input.dt 
    #to avoid double recordings, if allsec is in loclist, then loclist should only contain allsec
    for x in input.recordDict:         
        if 'allsec' in input.recordDict[x]: input.recordDict[x] = ['allsec'] 

    #------------------------------------------ CELL SETUP -----------------------------------------
    # This is necessary because finitialize moves parts of the cell sometimes and it can then not be restored moveSomaToOrigin
    OStim, optlocs, estims, cell, cl, k, na, cyt, ecs  = setupCell(input)
    
    #_______________________________________SETUP RECORDINGS________________________________________
    Traces = Dict()

    #------- time signal -------
    Traces.t = h.Vector().record(h._ref_t, input.dt)

    #------- electrical stimulation recording ------
    #injected current input
    Traces.iStim = h.Vector().record(estims[0]._ref_i, input.dt)

    #------- opsin recordings -------
    Traces.Iopto = h.Vector().record(OStim._ref_Io, input.dt)

    #------- Other recordings via recordDict ------- 
    KeyToPointer = {'v': 'v', 
                    'iopto': 'iopsin'+ '_' + input.cellsopt.opsin_options.opsinmech,
                    'Oopto':   'o' + '_' + input.cellsopt.opsin_options.opsinmech,
                    'DAopto': 'da' + '_' + input.cellsopt.opsin_options.opsinmech,
                    'gopto': 'g' + '_' + input.cellsopt.opsin_options.opsinmech,
                    'icl': 'icl', 'ecl': 'ecl', 'cli': 'cli', 'clo': 'clo',
                    'ik': 'ik', 'ek': 'ek', 'ki': 'ki', 'ko': 'ko',
                    'ina': 'ina', 'ena': 'ena', 'nai': 'nai', 'nao': 'nao',
                    'ik_kcc2': 'ik_' + input.cellsopt.KCC2type, 'ik_nakpump': 'ik_nakpump',
                    'tauO': 'tauo_'+ input.cellsopt.opsin_options.opsinmech, 'tauR': 'tauda_' + input.cellsopt.opsin_options.opsinmech, 'Rinf': 'dainf_'+ input.cellsopt.opsin_options.opsinmech, 'Oinf': 'oinf_' + input.cellsopt.opsin_options.opsinmech}
    
    #create new dictionary with sections as locations instead of region strings
    NewRecordDict = Dict()
    for x, loclist in input.recordDict.items(): NewRecordDict[x] = [seg for locstr in loclist for seg in getattr(cell, locstr)]

    #checks to avoid errors regarding non-existing pointers, traces for location that get deleted here will be set to zero numpys when results dict is assembles
    AddToResults = {'v': [], 'iopto': [], 'Iopto': [], 'Oopto': [], 'DAopto': [], 'gopto': [], 'cli': [], 'clo': [], 'ecl': [], 'icl': [], 'ki': [], 'ko': [], 'ek': [], 'ik': [], 'nai': [], 'nao': [], 'ena': [], 'ina': [], 'ik_kcc2': [], 'ik_nakpump': [],'tauO': [], 'tauR': [], 'Rinf': [], 'Oinf': []}
    for x in ['iopto', 'Iopto', 'Oopto', 'DAopto', 'gopto', 'cli', 'clo', 'ecl', 'icl', 'ik_kcc2', 'ik_nakpump']:
        if ('cl' not in x) or ('cl' not in input.cellsopt.ion_options.iontypes and input.cellsopt.opsin_options.opsinlocations != 'all'):
            if x in NewRecordDict.keys():
                for loc in NewRecordDict[x].copy(): 
                    if loc not in optlocs: 
                        NewRecordDict[x].remove(loc)
                        AddToResults[x].append(loc)
    
    #create list recordingpointers for each trace with a pointer for each location in the list
    for x, loclist in NewRecordDict.items(): Traces[x] = [h.Vector().record(getattr(sec(0.5), '_ref_' + KeyToPointer[x]), input.dt) for sec in loclist]
    
    #_______________________________________RUN SIMULATION____________________________________________
    #initialize neuron simulation
    if input.stimopt.Estimparams.vClamp is False:
        h.cvode.active(True)
        h.cvode.atol(1e-4)

    h.celsius = 35
    h.v_init = input.v0

    h.init()
    

    #set initial concentrations of ecs 
    for ion in input.cellsopt.ion_options.iontypes:
        nodes_in_ecs, nodes_in_cyt = locals()[ion].nodes(ecs), locals()[ion].nodes(cyt)
        for node in nodes_in_ecs: node.value = input.cellsopt.ion_options.cextra0[ion]
        for node in nodes_in_cyt: node.value = input.cellsopt.ion_options.cintra0[ion]


    #initialization of ion concentrations messes with placement of neuron so soma movement to origin done at this point, after initialisation
    if input.cellsopt.init_options.movesomatoorigin: cell.moveSomaToOrigin()
    
    #run simulation
    h.continuerun(input.duration)

    #_________________________________ADD IMPORTANT THINGS TO RESULTS OUTPUT (NO POINTERS!)_______________________________
    results = Dict()

    results.t = np.array(Traces.t)

    results.stimulations  = Dict()
    for x in ['Iopto', 'iStim']: results.stimulations[x] = np.array(Traces[x])

    results.recordings  = Dict()
    for name, loclist in NewRecordDict.items(): results.recordings[name] = {str(loclist[i]): np.array(vector) for i, vector in enumerate(Traces[name])}
    #add traces that were not recorded because of possible pointer errors
    for name, loclist in AddToResults.items(): 
        for loc in loclist: results.recordings[name][str(loc)] = np.zeros(results.t.size) * np.nan

    #save which segments belong to which region, handy for plotting
    results.regionsDict = {'somatic' : [], 'axon' : [], 'alldend' : [], 'apicalTrunk' : [], 'apicalTuft' : [], 'apicalTrunk_ext': [], 'apicalTuft': [], 'apical_obliques': [], 'basaldend': []}
    for region in results.regionsDict: 
        for seg in getattr(cell, region): results.regionsDict[region].append(str(seg)) 

    #________________________________________EXTRACT FEATURES_______________________________________
    results.features = Dict()

    if input.analysesopt.FeaturesLocOverride is not None: 
        for name in input.analysesopt.featuresDict.keys(): input.analysesopt.featuresDict[name] = input.analysesopt.FeaturesLocOverride

    NewFeaturesDict = Dict()
    for x, loclist in input.analysesopt.featuresDict.items() : NewFeaturesDict[x] = [str(sec) for locstr in loclist for sec in getattr(cell, locstr)]
    
    results.APs = Dict()
    for sec in cell.allsec: results.APs[str(sec)] = Features.getSpiketimes(results.t, results.recordings.v[str(sec)],threshold=-10)

    if 'OstimPower' in input.analysesopt.featuresDict.keys(): results.features['OstimPower'] = getattr(Features, 'OstimPower')(input, results, loc)
    for name, loclist in NewFeaturesDict.items(): 
        if name != 'OstimPower': results.features[name] = {str(loc): getattr(Features, name)(input, results, loc) for loc in loclist}
    
    #_____________________________________________SAVE_______________________________________________
    
    if input.save_flag:

        #set unique filename using datetime
        filename = input.resultsFolder + '/' + datetime.now().strftime("%Y%m%d_%H%M%S") + '_' + str(input.id)

        #save input data
        inputData = {}
        inputData['settings'] = input

        with open(filename + '_input.json', 'w') as outfile: json.dump(inputData, outfile, indent=4, signif=input.signif, cls=MyEncoder)

        #save results in different subdictionaries
        Data = {}

        Data['runtime'] = time() - StartTime 

        if 't' in input.saveDict: Data['t'] = results.t

        if 'regionsDict' in input.saveDict: Data['regionsDict'] = results.regionsDict

        Data['stimulations']  = {}
        if 'iStim' in input.saveDict: Data['stimulations']['iStim'] = results.stimulations.iStim
        if 'Iopto' in input.saveDict: Data['stimulations']['Iopto'] = results.stimulations.Iopto

        Data['recordings']  = {}
        for name, loclist in NewRecordDict.items():
            if name in input.saveDict.keys():
                savelocs = [seg  for i in range(len(input.saveDict[name])) for seg in getattr(cell, input.saveDict[name][i])]
                Data['recordings'][name] = {}
                for i, vector in enumerate(Traces[name]):
                    if loclist[i] in savelocs: Data['recordings'][name][str(loclist[i])] = np.array(vector)
        
        Data['features'] = {}
        for name in results.features: Data['features'][name] = results.features[name]

        with open(filename + '_results.json', 'w') as outfile: json.dump(Data, outfile, indent=4, signif=input.signif,  cls=MyEncoder)

    del cell

    if input.returnResults: return results
    else: return {}



