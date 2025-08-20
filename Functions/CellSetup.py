
import Model.Cells as Cell
from neuron import h, rxd, load_mechanisms
from Functions.utils import Dict
from random import randint

def setupCell(input):
    OStim, optlocs, estims, cell, cl, k, na, cyt, ecs = None, None, None, None, None, None, None, None, None
    cell = getattr(Cell, input.cellsopt.neurontemplate)(movesomatoorigin = False, replace_axon = input.cellsopt.init_options.replace_axon, morphologylocation = input.cellsopt.init_options.morphologylocation)
    
    #---------------------------------------DEFINE DIFFUSION REGIONS AND IONS-------------------------
    # define cytosol. MUST specify nrn_region for concentrations to update
    if len(input.cellsopt.ion_options.iontypes)>0: cyt = rxd.Region(cell.allsec, name='cyt', nrn_region='i')

    #define extracellular diffusion space
    xlo, ylo, zlo, xhi, yhi, zhi = cell.extrema()
    padding = 50
    if len(input.cellsopt.ion_options.iontypes)>0: ecs = rxd.Extracellular( xlo - padding, ylo - padding, zlo - padding,
                                                                            xhi + padding, yhi + padding, zhi + padding,
                                                                            dx=(20, 20, 20), volume_fraction=0.2, tortuosity=1.6)

    # add different ions as species in cyt and ecs, connection automatically made to currents in mod files
    if 'cl' in input.cellsopt.ion_options.iontypes: cl = rxd.Species([cyt, ecs], d=input.cellsopt.ion_options.diffc['cl'], name='cl', charge=-1)
    if 'k'  in input.cellsopt.ion_options.iontypes: k  = rxd.Species([cyt, ecs], d=input.cellsopt.ion_options.diffc['k'], name='k', charge=1)
    if 'na' in input.cellsopt.ion_options.iontypes: na = rxd.Species([cyt, ecs], d=input.cellsopt.ion_options.diffc['na'], name='na', charge=1)

    #--------------------------------------- Optional ion channels/pumps -----------------------
    if input.cellsopt.insertKCC:
        for sec in cell.allsec: 
            if 'axon' not in str(sec): 
                sec.insert(input.cellsopt.KCC2type)
    if input.cellsopt.insertNaKpump: 
        for sec in cell.alldend: 
            sec.insert('nakpump')
            for seg in sec: setattr(getattr(seg, 'nakpump'), 'imax', 0.009)
        for sec in cell.somatic: 
            sec.insert('nakpump')
            for seg in sec: setattr(getattr(seg, 'nakpump'), 'imax', 0.014)

    #---------------------------------------INSERT OPTOGENETICS---------------------------------------
    # get opsin locations from list of strings
    if input.cellsopt.opsin_options.opsinlocations == 'all': optlocs = cell.allsec
    else: optlocs = [sec for locstr in input.cellsopt.opsin_options.opsinlocations for sec in getattr(cell, locstr)]
   
    #create optic stimulation input
    OStim = h.OPTO_pulse(getattr(cell, 'somatic')[0](0.5))
    OStim.osstart = input.stimopt.Ostimparams.delay
    OStim.osdur = input.stimopt.Ostimparams.dur
    OStim.Amp = input.stimopt.Ostimparams.amp #W/m2
    OStim.osdc = input.stimopt.Ostimparams.dc
    OStim.osprp = input.stimopt.Ostimparams.prp

    #insert opsin
    for l in optlocs:

        l.insert(input.cellsopt.opsin_options.opsinmech) #insert mechanism at location

        for seg in l:
            h.setpointer(OStim._ref_Io, 'Iopto', getattr(seg, input.cellsopt.opsin_options.opsinmech)) #connect optic pulse to every section

            #change opsin parameters
            setattr(getattr(seg, input.cellsopt.opsin_options.opsinmech), 'gmax', input.cellsopt.opsin_options.gmax)
            setattr(getattr(seg, input.cellsopt.opsin_options.opsinmech), 'Erev', input.cellsopt.opsin_options.Erev)
            setattr(getattr(seg, input.cellsopt.opsin_options.opsinmech), 'EvarFlag', input.cellsopt.opsin_options.EvarFlag)
            if 's' in input.cellsopt.opsin_options.opsinmech:
                setattr(getattr(seg, input.cellsopt.opsin_options.opsinmech), 'Iratio', input.cellsopt.opsin_options.Iratio)
                setattr(getattr(seg, input.cellsopt.opsin_options.opsinmech), 'tauOn' , input.cellsopt.opsin_options.tauOn)
                setattr(getattr(seg, input.cellsopt.opsin_options.opsinmech), 'tauOff' , input.cellsopt.opsin_options.tauOff)
                setattr(getattr(seg, input.cellsopt.opsin_options.opsinmech), 'tauInact', input.cellsopt.opsin_options.tauInact)
                setattr(getattr(seg, input.cellsopt.opsin_options.opsinmech), 'tauRecov', input.cellsopt.opsin_options.tauRecov)
                setattr(getattr(seg, input.cellsopt.opsin_options.opsinmech), 'Oinf_p1', input.cellsopt.opsin_options.Oinf_p1)
                setattr(getattr(seg, input.cellsopt.opsin_options.opsinmech), 'Oinf_p2', input.cellsopt.opsin_options.Oinf_p2)

    #---------------------------------------INSERT ELECTRICAL STIMULATION----------------------------------------

    # .......... create current injections in different estimlocs .........

    # get stimulation locations from list of strings
    estimlocs = [sec for locstr in input.stimopt.Estimparams.estimlocations for sec in getattr(cell, locstr)]

    input.stimopt.Estimparams.estim_amp_uAcm2 = input.stimopt.Estimparams.amp/((1e-5)*estimlocs[0](0.5).area())

    if 'currentInjection' not in input.stimopt.Estimparams.method: input.stimopt.Estimparams.amp = 0

    estims = []
    for estimloc in estimlocs:
        if not input.stimopt.Estimparams.vClamp:
            estims.append(h.IClamp(estimloc(0.5))) 
            #set currentInjection parameters
            estims[-1].dur = input.stimopt.Estimparams.dur
            estims[-1].delay = input.stimopt.Estimparams.delay
            estims[-1].amp = input.stimopt.Estimparams.amp
        else:
            estims.append(h.VClamp(estimloc(0.5)))
            estims[-1].dur[0], estims[-1].dur[1], estims[-1].dur[2] = input.stimopt.Estimparams.delay, input.stimopt.Estimparams.dur, input.duration
            estims[-1].amp[0], estims[-1].amp[1], estims[-1].amp[2] = 0, input.stimopt.Estimparams.amp, 0
    
    return OStim, optlocs, estims, cell, cl, k, na, cyt, ecs