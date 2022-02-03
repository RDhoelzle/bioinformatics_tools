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
__script_name__ = 'vcf_to_msa.py'
__version__ = '0.0.0'
__profiling__ = 'True'

import sys,os
import pandas as pd
from Bio import SeqIO
import vcf
import io
import gzip
import subprocess
from datetime import datetime
import argparse


###############################################################################
###############################################################################

## Function definitions

def fasta_alignment_from_vcf(vcf_file):
    """Get snp site alt bases as sequences from all samples in a vcf file"""

    import vcf
    from collections import defaultdict
    vcf_reader = vcf.Reader(open(vcf_file, 'rb'))
    def default():
        return []
    result = defaultdict(default)
    sites = []
    for record in vcf_reader:
        ref = record.REF
        result['ref'].append(record.REF)
        sites.append(record.POS)
        for sample in record.samples:
            name = sample.sample
            if sample.gt_bases != None:
                result[name].append(sample.gt_bases)
            else:
                result[name].append(record.REF)
    print (datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"), 'found %s sites' %len(sites))
    recs = []
    for sample in result:
        seq = ''.join(result[sample])
        seqrec = SeqRecord(Seq(seq),id=sample)
        recs.append(seqrec)

    return recs

###############################################################################
###############################################################################

## Main script

def main(vcf, msa):
    """
    This is what we're here to do
    """
    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),'Extracting multi sequence alignment')

    msa_file=fasta_alignment_from_vcf(vcf)

    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),'Writing multi sequence alignment to output fasta file')
    #write out to fasta_out
    SeqIO.write(msa_file, msa, "fasta")

    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),'Success!')

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='vcf_to_msa.py')
    parser.add_argument('-i', dest='vcf', required=True,
                        help='Input vcf file, must be .vcf'),
    parser.add_argument('-o', dest='msa', required=True,
                        help='Output msa file, must be .fasta')

    args = parser.parse_args()

    if not os.path.exists(args.vcf):
        print("Input vcf file doesnt exist: {}".format(args.vcf))
        sys.exit(-1)

    main(args.vcf, args.msa)
