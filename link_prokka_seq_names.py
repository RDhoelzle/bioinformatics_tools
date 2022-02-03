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
__script_name__ = 'link_prokka_seq_names.py'
__version__ = '0.0.1'
__profiling__ = 'False'

import os
import sys
import tempfile
import io
import csv
from datetime import datetime
import argparse
from collections import defaultdict

###############################################################################
###############################################################################

## Main script

def main(tbl, out_file):
    """
    This is what we're here to do
    """
    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),'Building table')
    #define line identifiers
    next_entry = '>Feature'
    prok_id = 'locus_tag'

    #initiate output table
    out = defaultdict(list)

    #define initial vars
    proks = []
    i = 0
    k = 0

    with open(tbl, 'r') as file:
        for line in file:
            if next_entry in line:
                if (i == 0):
                    otus = line.strip('\n').split()[1]
                    i = i+1
                else:
                    out[otus].append(proks)
                    otus = line.strip('\n').split()[1]
                    proks = []
                    i = i+1

            elif prok_id in line:
                if (k < i):
                    proks = line.strip('\n').split()[1]
                    k = i

                elif (k == i):
                    proks = proks + ";" + line.strip('\n').split()[1]

    out[otus].append(proks)

    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),'Writing IDs to output file')
    #write out to tsv
    w = csv.writer(open(out_file, "w"), delimiter='\t')
    for key, val in out.items():
        w.writerow([key, val])

    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),'Success!')

###############################################################################
###############################################################################

## Arguments

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='link_prokka_seq_names.py')
    parser.add_argument('-i', dest='tbl_file', required=True,
                        help='Prokka table file (.tbl)'),
    parser.add_argument('-o', dest='out_file', required=True,
                        help='Name of output file (.tsv)')

    args = parser.parse_args()
    if not os.path.exists(args.tbl_file):
        print("Input list file does not exist: {}".format(args.tsv_file))
        sys.exit(-1)

    main(args.tbl_file, args.out_file)
