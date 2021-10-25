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
__script_name__ = 'bamm_to_coverage_table.py'
__version__ = '0.0.0'
__profiling__ = 'True'

import os
import sys
import tempfile
import io
import re
import pandas as pd
import numpy as np
from datetime import datetime
import argparse

###############################################################################
###############################################################################

## Function definitions

def parse_file(file, covg_table, regex):
    """
    Calculate genome coverage, and parse MAG and sample name from each file
    """
    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),'Reading', file, 'coverage data')
    #read file to dataframe
    df = pd.read_csv(file, sep='\t', header=0)

    #regex match in <bam_file_name> column
    re_match = re.search(regex, df.columns[2])

    #pull sample name and mag name
    sample = re_match.group(0)
    mag = df.columns[2][0:re_match.start()-1]

    #calc contig total coverages
    df['total_covg'] = df['Length'] * df[df.columns[2]]

    #calc MAG size-normalised coverage
    covg = df['total_covg'].sum()/df['Length'].sum()

    #add MAG coverage to table (and MAG name if new)
    if mag in covg_table[0]:
        covg_table[0][mag].append(covg)
    else:
        covg_table[0][mag]=list()
        covg_table[0][mag].append(covg)

    #add sample name if new
    if sample not in covg_table[1]:
        covg_table[1].append(sample)

    return covg_table

###############################################################################
###############################################################################

## Main script

def main(tsvs, outfile):
    """
    This is what we're here to do
    """
    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),'Initialising coverage table')
    #initialize an empty dictionary and column/row lists
    mag_table = {}
    sample_names = []
    build_table = [mag_table, sample_names]

    #set sample name regex, looking for "AA##_S##"
    samp_name = re.compile('AA[0-9]{2}_S[0-9]{2}')

    for bam_out in tsvs:
        build_table = parse_file(bam_out, build_table, samp_name)

    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),'Generating and saving final coverage table')
    #convert dictionary to dataframe
    covg_df = pd.DataFrame(build_table[0], index=build_table[1])

    #write to tsv
    covg_df.to_csv(outfile, sep='\t', index=True)
    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),'Success!')

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='bamm_to_coverage_table.py')
    parser.add_argument('tsv_files', metavar='-b', nargs='+',
                        help='bamm parse output coverage file(s), tsv format, where each file is one genome and one sample'),
    parser.add_argument('outfile', metavar='-o',
                        help='Name of output coverage file, tsv format'),

    args = parser.parse_args()

    full_paths = [os.path.join(os.getcwd(), path) for path in args.tsv_files]
    files = list()
    for path in full_paths:
        if os.path.isfile(path):
            files.append(path)
        else:
            print("Input file does not exist: {}".format(path))
            sys.exit(-1)

    main(files, args.outfile)
