import inspect
from datetime import datetime

import neuron

from Functions.utils import Dict, replaceDictODict

class simParams(object):
    def __init__(self, simParamsDict=None):

        self.recordDict = {'v': ['allsec'], 
                                'ecl': ['allsec'], 'icl': ['allsec'], 'cli': ['allsec'], 'clo': ['allsec'],
                                'ek': ['allsec'],  'ik':  ['allsec'],  'ki': ['allsec'], 'ko':  ['allsec'],
                                'ena': ['allsec'], 'ina': ['allsec'], 'nai': ['allsec'], 'nao': ['allsec'],
                                'iopto': ['allsec'], 'gopto': ['allsec'], 'DAopto': [], 'Oopto': [],
                                'tauO': [], 'tauR': [], 'Rinf': [], 'Oinf': [],
                                'ik_kcc2': ['somatic'], 'ik_nakpump': ['somatic']}
        self.axlims = False
        
        self.returnResults = True
        
        self.save_flag = True
        self.saveDict = {'t': None, 'regionsDict': None,'v': ['somatic']}
        self.resultsFolder = 'Results'
        self.signif = 4 
        self.id = 0 #necessary so files don't get the same name when multiprocessing used

        self.seed = datetime(2021, 12, 1, 13, 0).timestamp()

        self.duration = 2500 #simulation time in ms
        self.dt = 0.1 #integration step of simulation solver
        self.defaultThreshold = -10 #firing threshold in mV
        self.v0 = -65 #initial membrane potential

        self.cellsopt = cellsOptions()
        self.stimopt = stimOptions()
        self.analysesopt = analysesOptions()

        if simParamsDict:
            # overwrite defaults
            netParamsComponents = ['cellsopt', 'stimopt', 'analysesopt']
            for k, v in simParamsDict.items():
                if k in netParamsComponents:
                    recursiveDictUpdate(getattr(self, k), v)
                else:
                    if hasattr(self, k):
                        if isinstance(v, dict):
                            setattr(self, k, Dict(v))
                        else:
                            setattr(self, k, v)
                    else:
                        raise AttributeError(k)

    def todict(self, reduce=None, inplace=False):
        if not inplace:
            out = replaceDictODict(self.__dict__).copy()
            out = replaceNeuronSectionsandFunTostr(out).copy()
        else:
            out = replaceDictODict(self.__dict__)
            out = replaceNeuronSectionsandFunTostr(out)
        return out


class cellsOptions(Dict):
    """
    Class to hold options form Model/Cells
    """

    def __init__(self):
        super().__init__()
        # for other options see cells.Neurontemplates
        self.neurontemplate = 'CA1PYR_TK21' #neurontemplate subclass as defined in Cells
        self.insertKCC = True
        self.insertNaKpump = True
        self.KCC2type = 'kcc2'

        # initialization options
        # options of the neurontemplate class
        self.init_options = Dict()
        self.init_options.replace_axon = True
        self.init_options.morphologylocation = './Model/morphologies'
        self.init_options.movesomatoorigin = True

        # opsin options
        self.opsin_options = Dict()
        self.opsin_options.opsinlocations = 'all' #'all' or a sublist of ['somatic', 'alldend', 'apicalTrunk', 'apicalTrunk_ext', 'apicalTuft', 'apical_obliques', 'basaldend']
        self.opsin_options.opsinmech = 'Cl22OM' 

        self.opsin_options.gmax = 1e-3 # conductance of opsin channel, in S/cm2,
        self.opsin_options.Erev = -65.95 # value as calculated via fitting procedure
        self.opsin_options.EvarFlag = 1 # when this is set to 1, the reversal potential of the opsin channel will change depeding on the varying ion concentrations

        # ion concentrations options
        self.ion_options = Dict()
        self.ion_options.iontypes = ['cl',      'k',       'na']       
        self.ion_options.diffc =    {'cl':2.03, 'k': 1.96, 'na': 1.33} # diffusion coefficient per ion type

        self.ion_options.cextra0 =  {'cl':145,  'k': 4,    'na':145}   # initial extracellular concentration per ion type
        self.ion_options.cintra0 =  {'cl':7,    'k': 96,   'na': 10}   # initial intracellular concentration per ion type

        #parameters for simplified opsin model
        self.opsin_options.Iratio = 0.677 #0.9 #
        self.opsin_options.tauOn = 0.85 #0.017 #
        self.opsin_options.tauOff = 68.032 #140.014 #
        self.opsin_options.tauInact = 112.80518971118501 #236.82 #
        self.opsin_options.tauRecov = 5123.4942834484855 #4803 #
        self.opsin_options.Oinf_p1 = 2.473
        self.opsin_options.Oinf_p2 = 0.677

    def setParam(self, label, param, value):
        if label in self:
            d = self[label]
        else:
            return False

        d[param] = value

        return True

    def rename(self, old, new, label=None):
        return self.__rename__(old, new, label)


class stimOptions(Dict):
    def __init__(self):

        # optical stimulation parameters
        self.Ostimparams = Dict()
        self.Ostimparams.delay = 750 #start time of light pulse in ms
        self.Ostimparams.dur = 1000 #duration of optogenetic stimulation in ms
        self.Ostimparams.amp = 1000 #amplitude of optogenetic stimulation in W/m2
        self.Ostimparams.off = 0 #if different from 0, a repetititive light pulse is formed and this parameter defines the off duration in each cycle 
        self.Ostimparams.dc = 1 #if different from 0, a repetititive light pulse is formed and this parameter defines the on duration in each cycle 
        self.Ostimparams.prp = 0 #pulse repitition period

        # electrical stimulation parameters
        self.Estimparams = Dict()

        self.Estimparams.vClamp = False

        self.Estimparams.estimlocations = ['somatic'] #'all' or a sublist of ['somatic', 'alldend', 'apicalTrunk', 'apicalTrunk_ext', 'apicalTuft', 'apical_obliques', 'basaldend']
        self.Estimparams.method = ['currentInjection']

        self.Estimparams.delay = 250 # start time of electrical stimulation in ms
        self.Estimparams.dur = 2000 # duration of electrical stimulation in ms

        self.Estimparams.amp = 0.8 #amplitude of injected current in case of currentInjection

        self.Estimparams.netstimInterval = 30 # (mean) time between spikes of netstim spiketrain in ms
        self.Estimparams.netstimNoise = 0.1 # noise level of netstim spiketrain between 0 to 1


    def setParam(self, label, param, value):
        if label in self:
            d = self[label]
        else:
            return False

        d[param] = value

        return True

    def rename(self, old, new, label=None):
        return self.__rename__(old, new, label)


class analysesOptions(Dict):
    """
    Class to hold options form Model/Cells

    """

    def __init__(self):

        self.featuresDict = {   'FR':['allsec'], 
                                'APheight': ['allsec'], 'APbase': ['allsec'], 'APpeaktopeak': ['allsec'],
                                'tInhibit': ['allsec'],
                                'iOpto': ['allsec'], 'gOpto': ['allsec'],
                                'clErev': ['allsec'], 'clo': ['allsec'],  'cli': ['allsec'],
                                'kErev': ['allsec'], 'ko': ['allsec'], 'ki': ['allsec'], 
                                'cliMov': ['allsec'], 'cloMov': ['allsec'], 
                                'kiMov': ['allsec'], 'koMov': ['allsec'], 
                                'OstimPower': []}
        
        self.FeaturesLocOverride= None

        self.t_tr = 200

    def setParam(self, label, param, value):
        if label in self:
            d = self[label]
        else:
            return False

        d[param] = value

        return True

    def rename(self, old, new, label=None):
        return self.__rename__(old, new, label)


def replaceNeuronSectionsandFunTostr(obj):
    if type(obj) == list:
        for i, item in enumerate(obj):
            if type(item) in [list, dict]:
                item = replaceNeuronSectionsandFunTostr(item)
            elif type(item) == neuron.nrn.Section or type(item) == neuron.nrn.Segment:
                obj[i] = str(item)
            elif callable(item):
                return inspect.getsource(item)

    elif type(obj) == dict:
        for key, val in obj.items():
            if type(val) in [list, dict]:
                obj[key] = replaceNeuronSectionsandFunTostr(val)
            elif type(val) == neuron.nrn.Section or type(val) == neuron.nrn.Segment:
                obj[key] = str(val)
            elif callable(val):
                obj[key] = inspect.getsource(val)
    elif type(obj) == neuron.nrn.Section or type(obj) == neuron.nrn.Segment:
        return str(obj)
    elif callable(obj):
        return inspect.getsource(obj)

    return obj


def recursiveDictUpdate(obj, value, checkhasattr_flag=True):
    for k, v in value.items():
        # add in brackets parameter names that do not need to be checked typically dict with variable keys if dict with dict here else add on next line
        if isinstance(v, (dict, Dict)) and (not k in ['options']):
            if len(v) == 0:
                setattr(obj, k, Dict(v))
            else:
                if k in []:
                    recursiveDictUpdate(getattr(obj, k), Dict(
                        v), checkhasattr_flag=False)
                else:
                    recursiveDictUpdate(getattr(obj, k), Dict(
                        v), checkhasattr_flag=True)
        else:
            if v == 'removedToReduceMemory':
                v = None
            # elif isinstance(v,str) and any([x in v for x in ['= lambda','=lambda']]):
                # convert str(functions) back to functions/callables
            #    v = eval(v.split('=')[-1])
            if checkhasattr_flag:
                if hasattr(obj, k):
                    setattr(obj, k, v)
                else:
                    raise AttributeError(k)
            else:
                setattr(obj, k, v)
