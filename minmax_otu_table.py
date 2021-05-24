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
__script_name__ = 'minmax_otu_table.py'
__version__ = '0.0.0'
__profiling__ = 'True'

import os
import sys
import tempfile
import gzip
import csv
import io
import re
import pdb
import pandas as pd
import numpy as np
from datetime import datetime
import argparse

###############################################################################
###############################################################################

## Main script

def main(otufile, minmax):
    #Load otufile to dataframe
    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),'Importing otu table')
    df=pd.read_csv(otufile, compression='gzip', error_bad_lines=False)

    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),'Filtering by minmax value')
    #drop otus with minmax < minmax
    for otu in df.columns[1:].tolist():
        if int(df[otu].max()) < int(minmax):
            df=df.drop(labels=otu,axis=1)

    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),'Dropping empty sample rows')
    #drop indices with sum = 0
    drop = []
    for row in df.index:
        if sum(df.iloc[row,1:].tolist()) == 0:
            drop.append(row)

    df=df.drop(labels=drop,axis=0)

    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),'Writing filtered otu table file')
    #write final filtered dataframe to csv
    writefile=otufile[:len(otufile)-7] + ".minmax" + str(minmax) + ".csv.gz"
    df.to_csv(writefile, index=False, compression='gzip')

    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),'Success!')

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='coverage_avg_var.py')
    parser.add_argument('-i', dest='csv_file', required=True,
                        help='OTU table to filter, must be .csv.gz')
    parser.add_argument('-m', dest='minmax', required=True,
                        help='MinMax value to filter by')

    args = parser.parse_args()

    if not os.path.exists(args.csv_file):
        print("Input OTU table does not exist: {}".format(args.csv_file))
        sys.exit(-1)
    #pdb.set_trace()
    if args.csv_file[-7:] != '.csv.gz':
        print("OTU table is not .csv.gz format")
        sys.exit(-1)

    main(args.csv_file, args.minmax)
