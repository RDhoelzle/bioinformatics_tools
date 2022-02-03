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
__script_name__ = 'merge_graftm_out.py'
__version__ = '0.0.0'
__profiling__ = 'True'

import pandas as pd
import os
import sys
import io
from datetime import datetime
import argparse

###############################################################################
###############################################################################

## Main script
def main(files, out_file):
    """
    This is what we're here to do
    """
    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"), "Propagating otu table")
    #list all of the subdirectories in the GraftM output file
    sub_files = os.listdir(files)

    #cycle through all samples
    for sample in sub_files:
        file=files + sample + '/combined_count_table.txt'
        #verify sub_file contains "combined_count_table.txt"
        if not os.path.exists(file):
            print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"), "WARNING: No GraftM results found in {}".format(sample))
        else:
            #read file to dataframe and drop ID column
            df = pd.read_csv(file, sep='\t', header=0)
            df.drop(['#ID'], axis=1, inplace=True)
            #Either create the otu_tab dataframe, or outer-join the current sample to it
            try:
                otu_tab
            except NameError:
                otu_tab = df
            else:
                otu_tab = pd.merge(otu_tab, df, on='ConsensusLineage', how='outer')

    #set index to 'ConsensusLineage'
    otu_tab.set_index('ConsensusLineage', inplace=True)
    otu_tab.fillna(0, inplace=True)

    out = files + out_file
    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"), "Saving final otu table to {}".format(out))
    otu_tab.to_csv(out, sep='\t')

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='run_analysis.py')
    parser.add_argument('-d', dest='directory', required=True,
                        help='[REQUIRED] GraftM results directory'),
    parser.add_argument('-o', dest='out_file', required=True,
                        help='[REQUIRED] Name of output otu tsv file')
    args = parser.parse_args()

    #add '/' to directory name if necessary
    if args.directory.endswith('/'):
        pass
    else:
        args.directory = args.directory + '/'

    #add '.tsv' to directory name if necessary
    if args.out_file.endswith('.tsv'):
        pass
    else:
        args.out_file = args.out_file + '.tsv'

    #verify GraftM results directory exists
    if not os.path.exists(args.directory):
        print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"), "ERROR: GraftM results directory doesnt exist: {}".format(args.directory))
        sys.exit(-1)

    main(args.directory, args.out_file)
