#!/usr/bin/env python3

#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

__maintainer__ = 'Rob Hoelzle'
__script_name__ = 'coverage_avg_var.py'
__version__ = '0.0.1'
__profiling__ = 'False'

import os
import sys
import argparse
import pandas as pd

###############################################################################
###############################################################################

## Function definitions

def running_vals(covg, oldAvg, oldVar, N):
    """
    Calcs running average and variance for current contig
    """
    newAvg = (oldAvg*(N-1) + covg)/N
    newVar = ((N-2)*oldVar + (covg-newAvg)*(covg-oldAvg))/(N-1)
    newVals = [newAvg, newVar]

    return newVals

def parse_contigs(file):
    """
    This function reads in the contig coverage file, parses positional
    coverage values for each unique contig to the running average and variance
    calculator, and stores the final average and variance as a dictionary entry
    """
    outfile = file.replace('coverage.tsv','avg_var_covg.tsv')
    i = 0 #first line indicator
    allVals = {} #define the empty dictionary
    newVals = [] #initialize avg and var list

    with open(file, "r") as data:
        print('Reading in raw data')
        for line in data:
            #breakpoint()
            entry = line.strip('\n').split('\t')
            if (int(entry[1]) == 1):
                if (i == 0): #no previous value to save for first contig
                    i = 1 #start positional counter
                    contig = entry[0] #capture contig name
                    newVals = [int(entry[2]), 0] #initialize avg and var list
                else:
                    allVals[contig] = newVals #store new dictionary entry
                    i = 1 #reset positional counter
                    contig = entry[0] #capture contig name
                    newVals = [int(entry[2]), 0] #initialize avg and var list
            else:
                i = i+1 #increase positional counter
                newCovg = int(entry[2]) #caputre positional coverage
                #calc new avg and var values
                newVals = running_vals(newCovg, newVals[0], newVals[1], i)

    print('Finalizing table')
    allVals[contig] = newVals #store final dictionary entry
    print('Writing output file')
    pd.DataFrame(allVals).to_csv(outfile, sep='\t', index = False)


###############################################################################
###############################################################################

## Main script

def main(tsv):
    parse_contigs(tsv)
    print('Success!')

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='coverage_avg_var.py')
    parser.add_argument('tsv_file'),

    args = parser.parse_args()

    if not os.path.exists(args.tsv_file):
        print("Input coverage file doesnt exist: {}".format(args.tsv_file))
        sys.exit(-1)


    main(args.tsv_file)
