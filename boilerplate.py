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
__script_name__ = 'Whats in a name'
__version__ = '0.0.0'
__profiling__ = 'True'

import os
import sys
import tempfile
import shutil
import io
import gzip
import subprocess
from datetime import datetime
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import biom
import skbio

from config import dbs
from common import logger, AnalysisException, PickingStrategy, qiime_export, fastq_to_fasta

import usearch_utils

GLOBAL_LOGGER = logger()

###############################################################################
###############################################################################

## Function definitions

def log(msg):
    print(msg)
    GLOBAL_LOGGER.write(msg + "\n")

def phelp():
    print("""
    some helpful info about this script

  Author: %s
  Version: %s
""" % (__maintainer__,__version__))

def some_fn( args ):

    """
    give me some def
    """

    """
    except some errors
    try:
        #something
    except:
        print('Some description of error:')
        raise
    """

    #do some stuff

    pass

###############################################################################
###############################################################################

## Main script

def main(tsv_fp, db_name, trim_length, threads, itsx_profile, picking_strategy=None, use_uclust=False):
    pass

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='run_analysis.py')
    parser.add_argument('-i', dest='tsv_file', required=True,
                        help='Input analysis definition file (tsv format)'),
    parser.add_argument('-d', dest='db_name', required=True,
                        help='The reference database to use for the QIIME analysis. Current available options: ' +
                        ", ".join(sorted(dbs.known_db_names)))
    parser.add_argument('-l', dest='trim_length', type=int,
                        help='Trim sequences to this read length in bp. If omitted the pipeline will stop ' +
                        'after sample merging/demultiplexing and allow the user to choose a trim length.' ),
    parser.add_argument('-t', dest='threads', type=int,
                        help='Threads to use', default=1),
    parser.add_argument('-I', dest='itsx_type', required=False, default='F',
                        help='ITSx profile (organism groups) to use for searching - ' +
                             'See ITSx User Guide for more info. Default: fungi. (Ignored for non-ITS analyses)')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--denovo', dest='picking_strategy', action='store_const', const=PickingStrategy.DENOVO,
                        help="Perform the pipeline using a denovo OTU pick with the uparse pipeline (previously implemented method)")
    group.add_argument('--open', dest='picking_strategy', action='store_const', const=PickingStrategy.OPEN_REF,
                       help="Perform the pipeline using an open reference OTU pick with the uparse pipeline for the denovo stage")
    group.add_argument('--uclust', dest='use_uclust', action='store_true', default=False,
                        help="Use uclust instead of uparse for OTU picking. (Invalid for ITS analyses)"),


    args = parser.parse_args()
    if not os.path.exists(args.tsv_file):
        print("Input analysis definition file doesnt exist: {}".format(args.tsv_file))
        sys.exit(-1)

    use_uclust = False
    if args.use_uclust:
        use_uclust = True

    picking_strategy=None
    if args.picking_strategy is not None:
        picking_strategy=args.picking_strategy

    main(args.tsv_file, args.db_name, args.trim_length, args.threads, args.itsx_type, picking_strategy=picking_strategy, use_uclust=use_uclust)
