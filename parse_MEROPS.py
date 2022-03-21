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
__script_name__ = 'parse_MEROPS.py'
__version__ = '0.0.0'
__profiling__ = 'False'

import os
import sys
import argparse
import pandas as pd
import numpy as np
from datetime import datetime

###############################################################################
###############################################################################

## Main script

def main(file_list, outFile):
    """
    This is what where here to do
    """
    #Import ref files
    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),"Importing reference files")
    
    metad = '/safs-data02/dennislab/Public/public_dbs/merops/MERid_metadata.tsv'
    if not os.path.exists(metad):
        print("MERid metadata file does not exist: {}".format(metad))
        sys.exit(-1)
        
    peps = '/safs-data02/dennislab/Public/public_dbs/merops/PEPid_names.tsv'
    if not os.path.exists(metad):
        print("PEPid names file does not exist: {}".format(peps))
        sys.exit(-1)
    
    #MERids
    metaDict = {}
    with open(metad, "r") as file:
        for line in file:
            entry = line.strip('\n').split('\t')
            metaDict[entry[0]] = entry[1]
    
    #PEPids
    pepList = []
    with open(peps, "r") as file:
        next(file)
        for line in file:
            pepList.append(line.strip('\n').split('\t')[0])        
    
    #Read list of files
    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),"Reading input file names")
    fileList = []
    with open(file_list, "r") as file:
        for line in file:
            fileList.append(line.strip('\n').split('.merops')[0])

    #Generate blank dataframe
    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),"Initializing MEROPS frequency table")
    freqDF = pd.DataFrame(0, index=np.arange(len(pepList)), columns=fileList)
    freqDF.index = pepList
    freqDF.sort_index(inplace=True)

    #Add PERid counts for ID>=30% & EValue<=1e-05
    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),"Writing MEROPS frequency table")
    minID = 30.0
    maxEV = 1e-05
    for blast in fileList:
        with open(("".join((blast, ".merops"))), "r") as file:
            for line in file:
                entry = line.strip('\n').split('\t')
                if (float(entry[2]) >= minID) & (float(entry[10]) <= maxEV):
                    freqDF.loc[metaDict[entry[1]]][blast] += 1
                    
    #Write output frequency table
    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),"Writing output file: {}".format(outFile))
    freqDF.to_csv(outFile, sep='\t')
    
    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),"Finished!")

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='parse_MEROPS.py')
    parser.add_argument('-i', dest='file_list', required=True,
                        help='List of input blast.merops files (one per line, tsv format)'),
    parser.add_argument('-o', dest='outFile', required=True,
                        help='Output MEROPS frequency table name (tsv format)')
    args = parser.parse_args()
    
    if not os.path.exists(args.file_list):
        print("Input list file does not exist: {}".format(args.file_list))
        sys.exit(-1)
        
    if os.path.exists(args.outFile):
        print("Cannot overwrite output frequency table: {}".format(args.outFile))
        sys.exit(-1)

    main(args.file_list, args.outFile)
