"""
Add the virtual sites Oxygen only
"""
import sys
import os
import argparse
import numpy as np
import cclib as cc
from mol_data import AT_RVDW
from fchk_io import fchk_geom_parser


B2A = 0.529177210903
PROGNAME = os.path.basename(sys.argv[0])
# Values taken from 
# https://www.frontiersin.org/articles/10.3389/fchem.2020.00584/full
# 
LP = [0.3067, 111.3, 90]

def get_atom_connect(atom, ian, coord):
    """ 
    Given return a list with the connectivity
    Args:
        atom (int): atom index (base 0)
        ian (np.array(N, dtype=int)): atomic numbers 
        coord (np.array((N,3), dtype=float)): atomic coordinates (Natoms x 3)

    Returns:
        (list[int]): the indices of connected atoms
    """
    crd = np.array(coord)
    atnum = crd.shape[0]
    # dist = np.zeros((atnum, atnum))
    connect = []
    for i in range(atnum):
        if i == atom:
            continue
        dist = np.sqrt((crd[i, 0] - crd[atom, 0])**2 +
                       (crd[i, 1] - crd[atom, 1])**2 +
                       (crd[i, 2] - crd[atom, 2])**2)
        vdw_dist = AT_RVDW[int(ian[i])] + AT_RVDW[int(ian[atom])]
        if dist < vdw_dist*B2A:
                connect.append(i)
    return connect

def read_xyz(fname):
    """read the molecular coordinates from a xyz file

    Args:
        fname (str): filename

    Returns:
        (list, list): returns a tuple of atomic numbers and cartesian coordinates
    """
    with open(fname, 'r') as fopen:
        line = fopen.readline()
        nline = int(line)
        line = fopen.readline()
        atm = []
        crd = []
        for _ in range(nline):
            line = fopen.readline()
            tokens = line.split()
            atm.append(tokens[0])
            crd.append([float(x) for x in tokens[1:]])
    return (atm, crd)

def add_atom(apos, bpos, cpos, intcrd):
    """[summary]

    Args:
        apos ([type]): [description]
        bpos ([type]): [description]
        cpos ([type]): [description]
        intcrd ([type]): [description]

    Returns:
        [type]: [description]
    """
    norm = lambda x: x/np.sqrt(np.dot(x, x))
    r_dis = intcrd[0]
    t_ang = np.deg2rad(intcrd[1])
    p_ang = np.deg2rad(intcrd[2])
    d_two = np.array([r_dis*np.cos(t_ang),
                      r_dis*np.cos(p_ang)*np.sin(t_ang),
                      r_dis*np.sin(p_ang)*np.sin(t_ang)])
    m_mat = np.zeros((3, 3))
    # BC
    m_mat[:, 0] = -norm(cpos-bpos)
    # AB
    ad_d = apos - bpos
    # ABxBC
    m_mat[:, 2] = norm(np.cross(ad_d, m_mat[:, 0]))
    # ABxBCxBC
    m_mat[:, 1] = norm(np.cross(m_mat[:, 2], m_mat[:, 0]))
    atm_pos = cpos + np.einsum('ij,j->i', m_mat, d_two)
    return atm_pos

def orth_vec(vertatm, pos1, pos2):
    norm = lambda x: x/np.sqrt(np.dot(x, x))
    vec1 = pos1-vertatm
    vec2 = pos2-vertatm
    return vertatm + norm(np.cross(vec1, vec2))

def build_parser():
    """
    Build options parser.
    """
    par = argparse.ArgumentParser(prog=PROGNAME,
                                  formatter_class=argparse.RawTextHelpFormatter)
    # MANDATORY ARGUMENTS
    txt = "Strucure file (es. xyz, gaussian log file)"
    par.add_argument('strucfile', help=txt)
    txt = "Oxygen index"
    par.add_argument('oind', type=int, help=txt)

    # OPTIONAL ARGUMENTS
    par.add_argument('-p', '--print', action='store_true',
                     help='Print molecule')
    return par

def main():
    parser = build_parser()
    opts = parser.parse_args()
    filename, file_extension = os.path.splitext(opts.strucfile)
    cm5_a = False
    if file_extension == ".fchk":
        data = fchk_geom_parser(opts.strucfile)
        xyz = np.array(data['crd'])*B2A
        anum = data['atnum']
        if 'cm5' in data:
            print('CM5 available added to XYZ')
            cm5_a = True
            cm5 = data['cm5']
    else:
        data = cc.ccopen(opts.strucfile).parse()
        xyz = data.atomcoords[-1]
        anum = data.atomnos
    ofile = filename + '_lp.xyz'
    pt = cc.parser.utils.PeriodicTable()
    oindex = opts.oind - 1
    if anum[oindex] != 8:
        print('Not oxigen')
        sys.exit()
    cno = get_atom_connect(oindex, anum, xyz)
    if len(cno) != 1:
        print("Not a carbonyl oxygen")
        sys.exit()
    cindex = cno[0]
    cnn = get_atom_connect(cindex, anum, xyz)
    cnn.remove(oindex)
    if len(cnn) != 2:
        print(cnn)
        print("Not an SP2 Carbon. Exit")
        sys.exit()
    # Ghost atom on the norm to the sp2 plane
    ghost = orth_vec(xyz[cindex, :], xyz[cnn[0], :], xyz[cnn[1], :])
    # Add the LP
    lps = []
    lps.append(add_atom(ghost, xyz[cindex, :], xyz[oindex, :], LP))
    LP[2] += 180.
    lps.append(add_atom(ghost, xyz[cindex, :], xyz[oindex, :], LP))
    xyz = np.vstack((xyz, np.array(lps)))
    anum = np.hstack((anum, np.zeros(2)))
    # anum
    line = '{a:4s}{b[0]:12.6f}{b[1]:12.6f}{b[2]:12.6f}\n'
    if cm5_a:
        vscharge = cm5[oindex] / 2
        cm5[oindex] = 0.0
        cm5 = np.hstack((cm5,vscharge, vscharge))
        line = '{a:4s}{b[0]:12.6f}{b[1]:12.6f}{b[2]:12.6f}{c:12.6f}\n'
    toprnt = ''
    for i, item in enumerate(anum):
        symb = pt.element[int(item)] if item else 'X'
        if cm5_a:
            toprnt += line.format(a=symb, b=xyz[i, :], c=cm5[i])
        else:
            toprnt += line.format(a=symb, b=xyz[i, :])
    if opts.print:
        print(toprnt)
    with open(ofile, 'w') as fopen:
        fopen.write('{:d}\n\n'.format(anum.shape[0]))
        fopen.write(toprnt)

if __name__ == "__main__":
    main()











