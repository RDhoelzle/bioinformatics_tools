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
__script_name__ = 'parse_enrichm_classify.py'
__version__ = '0.0.0'
__profiling__ = 'True'

import os
import sys
import argparse
import numpy as np
import pandas as pd
from datetime import datetime
from collections import defaultdict

###############################################################################
###############################################################################

## Function definitions

def make_df(file):
    """
    This function reads in the module completeness file, removes the column
    headers, calculates the # steps missing, then generates a dictionary
    """
    raw_data = defaultdict(list)

    with open(file, "r") as file1:
        print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),'Reading in raw data')
        for line in file1:
            #ignore header and blank lines
            if 'Genome_name' in line.replace('\t',''):
                continue
            elif not line.strip():
                continue

            #parse genome name
            genome = line.strip('\n').split('\t')[0]
            if genome not in raw_data['Genome']:
                raw_data['Genome'].append(genome)

            #parse module name
            module =  line.strip('\n').split('\t')[1] + ":" + line.strip('\n').split('\t')[2]
            mod_miss = module+'~~missing'
            mod_comp = module+'~~complete'

            #calc and append missing and complete values
            step =  int(line.strip('\n').split('\t')[3])
            need =  int(line.strip('\n').split('\t')[4])
            missing = need - step #calc number steps missing
            complete =  float(line.strip('\n').split('\t')[5])

            raw_data[mod_miss].append(missing)
            raw_data[mod_comp].append(complete)

    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),'Generating raw data frame')
    raw_df = pd.DataFrame(raw_data).set_index('Genome')

    return raw_df



###############################################################################
###############################################################################

## Main script

def main(tsv):
    raw_df = make_df(tsv)

    #Write data to tsvs
    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),'Writing output files')

    missing_cols = [col for col in raw_df.columns if '~~missing' in col]
    complete_cols = [col for col in raw_df.columns if '~~complete' in col]

    raw_df[missing_cols].to_csv('mod_steps_missing.tsv', sep='\t')
    raw_df[complete_cols].to_csv('mod_percent_complete.tsv', sep='\t')

    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),'Success!')

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='parse_enrichm_classify.py')
    parser.add_argument('tsv_file'),

    args = parser.parse_args()

    if not os.path.exists(args.tsv_file):
        print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),"Input analysis definition file doesnt exist: {}".format(args.tsv_file))
        sys.exit(-1)


    main(args.tsv_file)
