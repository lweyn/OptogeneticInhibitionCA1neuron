import json
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
from _ctypes import PyObj_FromPtr


# select fraction of list with max value. Usefull for recording
def sellist(frac, lmax, dmax=10):
    mylist = np.random.permutation(lmax)
    listend = int(np.ceil(frac*len(mylist)))
    listend = min(listend, dmax)
    mylist = mylist[:listend]
    return list(mylist)


def points_in_cylinder(pt1, pt2, r, q):
    # https://stackoverflow.com/questions/47932955/how-to-check-if-a-3d-point-is-inside-a-cylinder
    vec = pt2 - pt1
    const = r * np.linalg.norm(vec)
    return np.dot(q - pt1, vec) >= 0 and np.dot(q - pt2, vec) <= 0 and np.linalg.norm(np.cross(q - pt1, vec)) <= const


def get_memoryusage(local_vars, print_flag=False):
    # local_vars = list(locals().items()) paste this line of code before use of function
    totalmem = 0
    for var, obj in local_vars:
        if print_flag:
            print(var, sys.getsizeof(obj))
        totalmem += sys.getsizeof(obj)
    print("Total memory in use = %5.2f" % totalmem)
    return totalmem


def get_size(obj, seen=None):
    """Recursively finds size of objects"""
    size = sys.getsizeof(obj)
    if seen is None:
        seen = set()
    obj_id = id(obj)
    if obj_id in seen:
        return 0
    # Important mark as seen *before* entering recursion to gracefully handle
    # self-referential objects
    seen.add(obj_id)
    if isinstance(obj, dict):
        size += sum([get_size(v, seen) for v in obj.values()])
        size += sum([get_size(k, seen) for k in obj.keys()])
    elif hasattr(obj, '__dict__'):
        size += get_size(obj.__dict__, seen)
    elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
        size += sum([get_size(i, seen) for i in obj])
    return size


def Inputlist_tobatchsets(inputFileNames, val, method='nrsets', saveFolder=None):
    '''
    Noticed not ideal use of resources when running sequential simulation based on list or self devided list.
    Function is to create batch simulaiton in order not one cpu to remain idle
    methods: nrsets => devide subdivede among nr sets
             nrinputs => create ceil(len(inputFileNames)/nrinputs) sets
    '''
    nsets = 1
    N = len(inputFileNames)
    if method == 'nrsets':
        nsets = val
    elif method == 'nrinputs':
        nsets = np.ceil(N/val)
    else:
        raise KeyError("Incorrect method")

    nsets = int(nsets)
    nipset = N/nsets  # n inputs per set
    Nipset = [np.floor(nipset).astype(int)]*nsets
    if np.floor(nipset) != nipset:
        rest = N-sum(Nipset)
        for i in range(rest):
            Nipset[i] += 1
    csum = [0]+list(np.cumsum(Nipset).astype(int))
    output = [inputFileNames[x:y] for x, y in zip(csum[:-1], csum[1:])]

    # test if correct
    Ntest = sum([len(x) for x in output])
    if Ntest != N:
        raise ValueError(f"test N={Ntest} is not equal to length input={N}")

    if saveFolder is not None:
        import os
        setnames = []
        path = os.path.join(saveFolder, 'batch')
        os.makedirs(os.path.join(saveFolder, 'batch'), exist_ok=True)
        for i, items in enumerate(output):
            filepath = os.path.join(path, f"set{i}.txt")
            setnames.append(filepath)
            with open(filepath, 'w') as f:
                for item in items:
                    f.write("%s\n" % item)
        setpath = os.path.join(path, "setnames.txt")
        with open(setpath, 'w') as f:
            for item in setnames:
                f.write("%s\n" % item)
    return output

# significant figures


def signif(x, p):
    x = np.asarray(x)
    x_positive = np.where(np.isfinite(x) & (x != 0), np.abs(x), 10**(p-1))
    mags = 10 ** (p - 1 - np.floor(np.log10(x_positive)))
    if abs(x) > 1:
        return np.round(x, p)
    else:
        return np.round(x * mags) / mags


def applysigniftoall(x, p):
    if isinstance(x, (float)):
        return signif(x, p)
    elif isinstance(x, (list, np.ndarray)):
        for ix in range(len(x)):
            if isinstance(x[ix], float):
                x[ix] = signif(x[ix], p)
            if isinstance(x[ix], (list, dict, np.ndarray)):
                x[ix] = applysigniftoall(x[ix], p)
    elif isinstance(x, dict):
        if len(x) > 1:
            for key, val in x.items():
                x[key] = applysigniftoall(val, p)
    return x


def Rx(theta):
    return np.array([[1, 0, 0],
                     [0, np.cos(theta), -np.sin(theta)],
                     [0, np.sin(theta), np.cos(theta)]])


def Ry(theta):
    return np.array([[np.cos(theta), 0, np.sin(theta)],
                     [0, 1, 0],
                     [-np.sin(theta), 0, np.cos(theta)]])


def Rz(theta):
    return np.array([[np.cos(theta), -np.sin(theta), 0],
                     [np.sin(theta), np.cos(theta), 0],
                     [0, 0, 1]])


def TaitBryan2rotMat(theta):
    """  Calculates Rotation Matrix given Tait-Bryan angles.
    theta = [phi,theta,psi] Tait-Bryan angles -> first roll then pitch then yawn (see plane starting at x-axis)"""
    R = np.dot(Rz(theta[2]), np.dot(Ry(theta[1]), Rx(theta[0])))
    return R


# Encoders for json dump

 # https://stackoverflow.com/questions/42710879/write-two-dimensional-list-to-json-file/42721412#42721412
 # https://stackoverflow.com/questions/13249415/how-to-implement-custom-indentation-when-pretty-printing-with-the-json-module

class NoIndent(object):
    """ Value wrapper. """

    def __init__(self, value):
        self.value = value


class MyEncoder(json.JSONEncoder):
    FORMAT_SPEC = '@@{}@@'
    regex = re.compile(FORMAT_SPEC.format(r'(\d+)'))

    def __init__(self, **kwargs):
        # Keyword arguments to ignore when encoding NoIndent wrapped values.
        ignore = {'cls', 'indent'}

        if 'signif' in kwargs.keys():
            self._signif = kwargs['signif']
            del kwargs['signif']
        else:
            self._signif = None
        # Save copy of any keyword argument values needed for use here.
        self._kwargs = {k: v for k, v in kwargs.items() if k not in ignore}
        super(MyEncoder, self).__init__(**kwargs)

    def default(self, obj):
        if isinstance(obj, np.ndarray) or isinstance(obj, NoIndent):
            return self.FORMAT_SPEC.format(id(obj))
        else:
            if isinstance(obj, float):
                if self._signif is not None:
                    return super(MyEncoder, self).default(float(signif(obj, self._signif)))
            return super(MyEncoder, self).default(obj)

    def iterencode(self, obj, **kwargs):
        format_spec = self.FORMAT_SPEC  # Local var to expedite access.

        # Replace any marked-up NoIndent wrapped values in the JSON repr
        # with the json.dumps() of the corresponding wrapped Python object.
        for encoded in super(MyEncoder, self).iterencode(obj, **kwargs):
            match = self.regex.search(encoded)
            if match:
                id = int(match.group(1))
                no_indent = PyObj_FromPtr(id)
                if isinstance(no_indent, np.ndarray):
                    if self._signif is not None:
                        no_indent = NoIndent(applysigniftoall(
                            no_indent, self._signif).tolist())
                    else:
                        no_indent = NoIndent(no_indent.tolist())
                json_repr = json.dumps(no_indent.value, **self._kwargs)
                # Replace the matched id string with json formatted representation
                # of the corresponding Python object.
                encoded = encoded.replace(
                    '"{}"'.format(format_spec.format(id)), json_repr)

            yield encoded


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

# ----------------------------------------------------------------------------
# Dict class (allows dot notation for dicts)
# modified from netpyne code
# ----------------------------------------------------------------------------


class Dict(dict):
    """
    Class for/to <short description of `netpyne.specs.dicts.Dict`>


    """

    __slots__ = []

    def __init__(*args, **kwargs):
        self = args[0]
        args = args[1:]
        if len(args) > 1:
            raise TypeError('expected at most 1 arguments, got %d' % len(args))
        if args:
            self.update(self.dotify(args[0]))
        if len(kwargs):
            self.update(self.dotify(kwargs))

    # only called if k not found in normal places
    def __getattr__(self, k):
        try:
            # Throws exception if not in prototype chain
            return object.__getattribute__(self, k)
        except AttributeError:
            try:
                return self[k]
            except KeyError:
                raise AttributeError(k)

    def __setattr__(self, k, v):
        try:
            # Throws exception if not in prototype chain
            object.__getattribute__(self, k)
        except AttributeError:
            try:
                self[k] = v
            except:
                raise AttributeError(k)
        else:
            object.__setattr__(self, k, v)

    def __delattr__(self, k):
        try:
            # Throws exception if not in prototype chain
            object.__getattribute__(self, k)
        except AttributeError:
            try:
                del self[k]
            except KeyError:
                raise AttributeError(k)
        else:
            object.__delattr__(self, k)

    def todict(self):
        return self.undotify(self)

    def fromdict(self, d):
        d = self.dotify(d)
        for k, v in d.items():
            self[k] = v

    def __repr__(self):
        keys = list(self.keys())
        args = ', '.join(['%s: %r' % (key, self[key]) for key in keys])
        return '{%s}' % (args)

    def dotify(self, x):
        # nested dict to nested Dict => at every level use dot notation possible
        if isinstance(x, dict):
            return Dict((k, self.dotify(v)) for k, v in x.items())
        elif isinstance(x, (list, tuple)):
            # initiate same class variable as type x [self.dotify(v) for v in x] is a generator object will be input argument
            return type(x)(self.dotify(v) for v in x)
        else:
            return x

    def undotify(self, x):
        if isinstance(x, dict):
            return dict((k, self.undotify(v)) for k, v in x.items())
        elif isinstance(x, (list, tuple)):
            return type(x)(self.undotify(v) for v in x)
        else:
            return x

    def __rename__(self, old, new, label=None):
        """
        old (string): old dict key
        new (string): new dict key
        label (list/tuple of strings): nested keys pointing to dict with key to be replaced;
            e.g. ('PYR', 'secs'); use None to replace root key; defaults to None

        returns: True if successful, False otherwse
        """

        obj = self
        if isinstance(label, (tuple, list)):
            for ip in range(len(label)):
                try:
                    obj = obj[label[ip]]
                except:
                    return False

        if old in obj:
            obj[new] = obj.pop(old)  # replace
            return True
        else:
            return False

    def rename(self, *args, **kwargs):
        self.__rename__(*args, **kwargs)

    # def __missing__(self, key):
    #    if key and not key.startswith('_ipython'):
    #        value = self[key] = Dict()
    #        return value

    def __getstate__(self):
        return self.todict()

    def __setstate__(self, d):
        self = self.fromdict(d)

# ------------------------------------------------------------------------------
# Replace Dict with dict
# modified from netpyne code
# ------------------------------------------------------------------------------


def replaceDictODict(obj):
    """
    Function for/to <short description of `netpyne.sim.utils.replaceDictODict`>

    Parameters
    ----------
    obj : <type>
        <Short description of obj>
        **Default:** *required*


    """

    if type(obj) == list:
        for item in obj:
            if type(item) == Dict or any(['Dict' in str(x) for x in item.__class__.__bases__]):
                item = item.todict()
            if type(item) in [list, dict]:
                replaceDictODict(item)

    elif type(obj) in [dict, Dict] or any(['Dict' in str(x) for x in obj.__class__.__bases__]):
        for key, val in obj.items():
            if type(val) == Dict or any(['Dict' in str(x) for x in val.__class__.__bases__]):
                obj[key] = val.todict()
            if type(val) in [list, dict]:
                replaceDictODict(val)

    return obj
