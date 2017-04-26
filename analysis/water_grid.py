#!/usr/bin/env python2.7

#----------------------------------------------------------------
# Water_grid: Simple script to calculate the oxygen water density
# along a molecular dynamics trajectory. The oxygen
# volume space is divided into voxels of specified size
#
# Author: Dr Gaetano Calabro and Dr David Mobley
#
# University of California, Irvine
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, see
# http://www.gnu.org/licenses/
#-----------------------------------------------------------------

import MDAnalysis
import argparse
import sys
import numpy as np



# create_map(opts)
#
# This function creates the oxygen water voxel map
#
# Parameters:
# -----------
# opts: Python Argparse object
#     The passed mandatory arguments and options

# Returns:
# -------
# cmap : Python dictionary
#    the dictionary containing the voxel indexes

def create_map(opts):

    univ = MDAnalysis.Universe(opts.parameter, opts.trajectory)

    frames = len(univ.trajectory)
    cmap = {}

    c_frame = 0

    oxygens = univ.select_atoms(opts.selection)
    o_num = oxygens.n_atoms

    print('\n')

    for ts in univ.trajectory:

        coords = oxygens.coordinates()
        c_ints = np.ceil(coords/opts.grid).astype(np.int64)

        check_wn = 0
        for c_int in c_ints:
            idx =  tuple(c_int)
            if not idx in cmap:
                cmap[idx]=1
            else:
                cmap[idx]+=1
            check_wn+=1

        if check_wn != o_num:
            raise ValueError("Error! Wrong oxygen atom numbers")

        ln =  "\rProcessing Frame: %d / %d" % (c_frame+1, frames)
        sys.stdout.write(ln)
        sys.stdout.flush()

        c_frame+=1

    print('\n')

    # Average
    cmap = {k:v/float(frames) for k,v in cmap.items()}

    return cmap



# write_density(cmap, opts)
#
# This function writes a pdb density file
# where coordinates are the voxel centers
# and the beta factor column contains
# the average number of oxygens found
# in each voxel along the trajectory
#
# Parameters:
# -----------
# cmap: Python dictionary
#     the dictionary containing the voxel indexes
#
# opts: Python Argparse object
#     The passed mandatory arguments and options

def write_density(cmap, opts):

    f = open(opts.output, 'w')

    count = 0
    at = 'H'

    for xi, yi, zi in cmap.keys():

        count+=1

        x=(opts.grid/2.0) * (2*xi - 1)
        y=(opts.grid/2.0) * (2*yi - 1)
        z=(opts.grid/2.0) * (2*zi - 1)

        bf = cmap[(xi,yi,zi)]

        ln = "ATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f  1.00%6.2f          %2s  " % (count%100000, at,
                                                                                               'DNS', 'A',
                                                                                               1,
                                                                                               x, y, z,
                                                                                               bf,
                                                                                               'H')

        f.write(ln+'\n')

    f.close()


# The main function

if ('__main__'==__name__):


    parser = argparse.ArgumentParser(description='Calculate the average number of oxygens water molecules in a grid space',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    parser.add_argument('parameter', type=str, metavar='FILE',
                        help='File containing the system topology')

    parser.add_argument('trajectory', type=str, metavar='FILE',
                        help='File containing the system trajectory')

    parser.add_argument('-s', '--selection', default='resname WAT and name O', type=str,
                        help='The MDAnalysis atom selection')


    parser.add_argument('-g', '--grid', default='1.0', type=float,
                        help='Grid size in Angstrom unit used to dived the space in voxels')


    parser.add_argument('-o', '--output', default='density.pdb', type=str, metavar='FILE',
                        help='File to output the calculate oxygen water density')

    opts = parser.parse_args()


    ln=40*'-'
    print(ln)
    print('Topology File: %s\nTrajectory File: %s\nOutput File: %s\n\n'% (opts.parameter, opts.trajectory, opts.output))
    print('Grid size: %4.2f A\nSelection: %s' % (opts.grid, opts.selection))

    cmap = create_map(opts)

    write_density(cmap, opts)
    print('Done...')
    print(ln)
