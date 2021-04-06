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
__script_name__ = 'subselect_fasta.py'
__version__ = '0.0.0'
__profiling__ = 'True'

import os
import sys
import tempfile
import io
from datetime import datetime
import argparse

###############################################################################
###############################################################################

## Main script

def main(tsv, fasta_in, fasta_out):
    """
    This is what we're here to do
    """
    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),'Searching sequence names in input fasta file')
    #initiate output fasta file
    out = []

    #set file read state
    #0: searching for entry
    #1: active entry
    state = 0

    #variable for next fasta entry
    next_entry = '>'

    #import tsv and iterate through accession numbers
    with open(tsv, 'r') as list:
        for entry in list:
            #search accession number in input fasta
            with open(fasta_in, 'r') as fasta:
                for line in fasta:
                    #searching for entry
                    if (state == 0):
                        if entry in line:
                            out.append(line)
                            state = 1
                    #active entry
                    elif (state == 1):
                        #check for start of next entry
                        if next_entry in line:
                            state = 0
                            break
                        #continue appending active entry
                        else:
                            out.append(line)

    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),'Writing subselected sequences to output fasta file')
    #write out to fasta_out
    with open(fasta_out, 'w') as writefile:
        for line in out:
            writefile.write(line)

    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),'Success!')

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='subselect_fasta.py')
    parser.add_argument('-t', dest='tsv_file', required=True,
                        help='List (tsv file) of sequence names to subselect'),
    parser.add_argument('-i', dest='i_fasta', required=True,
                        help='Sequenes to subselect from')
    parser.add_argument('-o', dest='o_fasta', required=True,
                        help='Name of output fasta file'),

    args = parser.parse_args()
    if not os.path.exists(args.tsv_list):
        print("Input list file does not exist: {}".format(args.tsv_file))
        sys.exit(-1)
    if not os.path.exists(args.i_fasta):
        print("Input fasta file does not exist: {}".format(args.db_name))
        sys.exit(-1)

    main(args.tsv_file, args.i_fasta, args.o_fasta)
