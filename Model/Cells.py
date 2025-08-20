import mmap
import re
from math import pi, sqrt
import numpy as np
import matplotlib.pyplot as plt

from neuron import h, load_mechanisms

NeuronTemplates = ['CA1PYR_TK21']

class NeuronTemplate:
    def __init__(self, templatepath, templatename, replace_axon=True, morphologylocation='./Model/morphologies', ID=0, ty=0, col=0, phi=0, theta=0, psi=0, movesomatoorigin=True, **kwargs):
        self.templatepath = templatepath
        self.templatename = templatename
        self.morphologylocation = morphologylocation
        self.replace_axon = str(replace_axon).lower()
        self.celsius = 34
        self.ID = ID
        self.ty = ty
        self.col = col
        self.phi = phi
        self.theta = theta
        self.psi = psi
        self.movesomatoorigin = movesomatoorigin
        self.morphology = ''
        
    def load_template(self):
        h.load_file(self.templatepath)  # Load cell info
        self.template = getattr(h, self.templatename)(
            self.replace_axon, self.morphologylocation)  # initialize cell
        try:
            # try to create allsec list (used in many methods)
            self.allsec = []
            for x in self.template.all:
                self.allsec.append(x)
        except Exception as e:
            print(e)
    
    def move_attributes(self):
        # cell is loaded under template attribute -> move level up
        for x in self.template.__dir__():
            if x != '__class__':
                setattr(self, x, getattr(self.template, x))
        self.template = None

    def move_Cell(self, rt):
        # move Cell with translate vector rt
        for section in self.allsec:
            for i in range(section.n3d()):
                xyz = np.array(
                    [section.x3d(i), section.y3d(i), section.z3d(i)])
                xyz = xyz+rt
                h.pt3dchange(i, xyz[0], xyz[1], xyz[2],
                             section.diam3d(i), sec=section)

    def moveSomaToOrigin(self):  
        soma_pos = np.array([[self.soma[0].x3d(i), self.soma[0].y3d(
        i), self.soma[0].z3d(i)] for i in range(self.soma[0].n3d())])
        soma_pos = np.mean(soma_pos, axis=0)

        self.move_Cell(-soma_pos)

    def extrema(self):
        """Give the bounding box that contains the cell"""
        xlo = ylo = zlo = xhi = yhi = zhi = None
        for sec in self.allsec:
            n3d = sec.n3d()
            xs = [sec.x3d(i) for i in range(n3d)]
            ys = [sec.y3d(i) for i in range(n3d)]
            zs = [sec.z3d(i) for i in range(n3d)]
            my_xlo, my_ylo, my_zlo = min(xs), min(ys), min(zs)
            my_xhi, my_yhi, my_zhi = max(xs), max(ys), max(zs)
            if xlo is None:
                xlo, ylo, zlo = my_xlo, my_ylo, my_zlo
                xhi, yhi, zhi = my_xhi, my_yhi, my_zhi
            else:
                xlo, ylo, zlo = min(xlo, my_xlo), min(ylo, my_ylo), min(zlo, my_zlo)
                xhi, yhi, zhi = max(xhi, my_xhi), max(yhi, my_yhi), max(zhi, my_zhi)
        return (xlo, ylo, zlo, xhi, yhi, zhi)

    
    def __del__(self): 
        for sec in self.allsec: h.delete_section(sec=sec)

class CA1PYR_TK21(NeuronTemplate):
    '''
    CA1 pyramidal cell
    loads './Model/CA1PYR_TK21.hoc' with templatename = 'CA1PYR_TK21'
    '''

    def __init__(self, **kwargs):
        super().__init__(templatepath='./Model/CA1PYR_TK21.hoc',
                         templatename='CA1PYR_TK21', **kwargs)  # init of parent class
        load_mechanisms('./Model/Mods_Tomko/', warn_if_already_loaded=False)
        self.load_template()
        h.define_shape()
        self.move_attributes()
        self.make_lists()
        self.celltype = 'SP_PC'

    def __str__(self):
        try:
            return f'compartCell_{self.__class__.__name__}_{self.ID}'
        except:
            return 'compartCell%d' % self.ID

    def __repr__(self):
        return self.__str__()
    
    def make_lists(self):
        # create some specific list of section -> facilitates access
        self.alldend = []

        self.somatic = []
        somatic = ['soma']

        self.axon = []
        axon = ['axon']

        self.apicalTrunk = []
        apicalTrunk = ['radTmed', 'radTprox', 'radTdist']

        self.apicalTrunk_ext = []
        apicalTrunk_ext = []

        self.apicalTuft = []
        apicalTuft = ['lm_thick2', 'lm_medium2', 'lm_thin2', 'lm_thick1', 'lm_medium1', 'lm_thin1']

        self.apical_obliques = []
        obliques = ['rad_t1', 'rad_t2', 'rad_t3']

        self.basaldend = []
        basaldend = ['oriprox1', 'oridist1_1', 'oridist1_2', 'oriprox2', 'oridist2_1', 'oridist2_2']

        for x in self.allsec:
            if (not 'soma' in str(x)) and (not 'axon' in str(x)):
                self.alldend.append(x)
            if any([y in str(x) for y in somatic]):
                self.somatic.append(x)
            if any([y in str(x) for y in axon]):
                self.axon.append(x)
            if any([y in str(x) for y in apicalTrunk]):
                self.apicalTrunk.append(x)
            if any([y in str(x) for y in apicalTrunk_ext]):
                self.apicalTrunk_ext.append(x)
            if any([y in str(x) for y in apicalTuft]):
                self.apicalTuft.append(x)
            if any([y in str(x) for y in obliques]):
                self.apical_obliques.append(x)
            if any([y in str(x) for y in basaldend]):
                self.basaldend.append(x)
            if 'dend' in str(x):
                self.basaldend.append(x)
