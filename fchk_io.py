# -*- coding: utf-8 -*-
"""
temporary functions to read data from fchk files
"""
import os
import re
from math import ceil
import numpy as np



NCOLS_FCHK_R = 5
NCOLS_FCHK_I = 6

def fchk_geom_parser(fname):
    """Function that parse the Gaussian Formated checkpoint file and get
    currect geometry

    Arguments:
        fname {str} -- the file name of the .fchk file
    """
    keys = {'nat': 'Number of atoms',
            'crd': 'Current cartesian coordinates',
            'ian': 'Atomic numbers',
            'atm': 'Real atomic weights',
            'cm5': 'CM5 Charges'
            }
    #name = os.path.split(fname)[-1][:-5]
    data = {}
    qtt = len(keys)
    with open(fname, 'r') as fopen:
        line = fopen.readline()
        while line:
            if line.startswith(keys['nat']):
                data['natoms'] = int(line.split()[-1])
                nmnum = data['natoms'] * 3 - 6
                qtt -= 1
            elif line.startswith(keys['crd']):
                nval = int(line.split()[-1])
                nline = ceil(nval/NCOLS_FCHK_R)
                tmp = []
                for _ in range(nline):
                    line = fopen.readline()
                    tmp.extend([float(x) for x in line.split()])
                tmp = np.array(tmp)
                data['crd'] = tmp.reshape(-1, 3)
                qtt -= 1
            elif line.startswith(keys['ian']):
                nval = int(line.split()[-1])
                nline = ceil(nval/NCOLS_FCHK_I)
                tmp = []
                for _ in range(nline):
                    line = fopen.readline()
                    tmp.extend([int(x) for x in line.split()])
                data['atnum'] = np.array(tmp)
                qtt -= 1
            elif line.startswith(keys['atm']):
                nval = int(line.split()[-1])
                nline = ceil(nval/NCOLS_FCHK_R)
                tmp = []
                for _ in range(nline):
                    line = fopen.readline()
                    tmp.extend([float(x) for x in line.split()])
                data['atmas'] = np.array(tmp)
                qtt -= 1
            elif line.startswith(keys['cm5']):
                nval = int(line.split()[-1])
                nline = ceil(nval/NCOLS_FCHK_R)
                tmp = []
                for _ in range(nline):
                    line = fopen.readline()
                    tmp.extend([float(x) for x in line.split()])
                data['cm5'] = np.array(tmp)
                qtt -= 1

            if qtt == 0:
                break
            line = fopen.readline()
    return data


