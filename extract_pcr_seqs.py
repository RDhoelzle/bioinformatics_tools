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
__script_name__ = 'extract_pcr_seqs.py'
__version__ = '0.0.1'
__profiling__ = 'False'

## Libraries

import os
import sys
import math
import regex
import argparse
from datetime import datetime

###############################################################################
###############################################################################

## Function definitions
def format_primer(primer, score, rev):
    """
    Formats primer sequence into searchable string
    primer: primer sequence
    rev: True if reverse primer, False if forward
    """
    #Remove all whitespace and prime(') indicators
    primer = primer.replace(" ","")
    primer = primer.replace("'","")
    primer = primer.replace("’","")
    
    #Check allowable primer formats
    string1 = regex.compile(r'^5-[A|T|C|G|R|Y|M|K|S|W|H|B|V|D|N]+-3$')
    string2 = regex.compile(r'^[A|T|C|G|R|Y|M|K|S|W|H|B|V|D|N]+$')
    
    if not (string1.match(primer) or string2.match(primer)):
        if rev:
            print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),"Incorrect r\
                  everse primer sequence: {}".format(args.rev_pr))
            sys.exit(-1)
        else:
            print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),"Incorrect f\
                  orward primer sequence: {}".format(args.fwd_pr))
            sys.exit(-1)
    
    #Clean directional indicators
    if string1.match(primer):
        primer = primer.replace("5-","")
        primer = primer.replace("-3","")
    
    #Get primer length
    len_pr = len(primer)
    len_max = len_pr - math.ceil(len_pr*(score/100))

    #Sub alternative bases
    primer = regex.subf(r'R', '[G|A]', primer)
    primer = regex.subf(r'Y', '[T|C]', primer)
    primer = regex.subf(r'M', '[A|C]', primer)
    primer = regex.subf(r'K', '[G|T]', primer)
    primer = regex.subf(r'S', '[G|C]', primer)
    primer = regex.subf(r'W', '[A|T]', primer)
    primer = regex.subf(r'H', '[A|C|T]', primer)
    primer = regex.subf(r'B', '[G|T|C]', primer)
    primer = regex.subf(r'V', '[G|C|A]', primer)
    primer = regex.subf(r'D', '[G|A|T]', primer)
    primer = regex.subf(r'N', '[G|A|T|C]', primer)
        
    #Reverse compliment
    if rev:
        primer = primer[::-1]
        primer = regex.subf(r'A', '1', primer)
        primer = regex.subf(r'T', '2', primer)
        primer = regex.subf(r'C', '3', primer)
        primer = regex.subf(r'G', '4', primer)
        primer = regex.subf(r'\]', '5', primer)
        primer = regex.subf(r'\[', '6', primer)
        primer = regex.subf(r'1', 'T', primer)
        primer = regex.subf(r'2', 'A', primer)
        primer = regex.subf(r'3', 'G', primer)
        primer = regex.subf(r'4', 'C', primer)
        primer = regex.subf(r'5', '[', primer)
        primer = regex.subf(r'6', ']', primer)
    
    fuzz_ptn = regex.compile('(%s){e<=%s}' % (primer, len_max))

    return fuzz_ptn

###############################################################################
###############################################################################

## Main script
def main(i_fa, o_fa, fwd_pr, rev_pr, m_scr):
    """
    This is what we're here to do
    """
    #format primers and mismatches for fuzzy matching
    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),'Verifying primer formats')
    fwd_ptn = format_primer(fwd_pr, m_scr, False)    
    rev_ptn = format_primer(rev_pr, m_scr, True)
    
    #initiate result fasta
    out = []
    
    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),'Trimming %s' % i_fa)
    #search for pcr seqs and print to result fasta
    with open(i_fa, 'r') as fasta:
        for line in fasta:
            #searching for entry
            if line[0] == '>':
                out.append(line)
            elif not (fwd_ptn.search(line) and rev_ptn.search(line)):
                out.pop()
            else:
                line = regex.subf(fwd_ptn, '1', line, regex.BESTMATCH).split('1')[-1]
                line = regex.subf(rev_ptn, '2', line, regex.BESTMATCH).split('2')[0] + '\n'
                out.append(line)
                
    #write output file
    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),'Writing %s' % o_fa)
    with open(o_fa, 'w') as writefile:
        for line in out:
            writefile.write(line)
            
    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),'Finished')

###############################################################################
###############################################################################

## Argparse
if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='extract_pcr_seqs.py')
    
    #import fasta file(s)
    parser.add_argument('-i', dest='fa_in', required=True, nargs='+',
                        help='Fasta file to extract pcr sequences from (<file>\
                            .fasta'),
    #define output tag
    parser.add_argument('-o', dest='fa_out', required=True,
                        help='Fasta file to print extracted sequences to (<fil\
                            e>.fasta)')
    #forward primer sequence
    parser.add_argument('-f', dest='f_primer', required=True,
                        help='Forward primer sequence, formatted as "5’-AAA CT\
                            Y AAA KGA ATT GRC GG-3’", "AAA CTY AAA KGA ATT GRC\
                                GG", "AAACTYAAAKGAATTGRCGG". Allowable bases: \
                                    https://genome.ucsc.edu/goldenPath/help/iupac.html')
    #reverse primer sequence
    parser.add_argument('-r', dest='r_primer', required=True,
                        help='Reverse primer sequence, formatted as "5’-ACG GG\
                            C GGT GWG TRC-3’", "ACG GGC GGT GWG TRC", or "ACGG\
                                GCGGTGWGTRC". Allowable bases: \
                                    https://genome.ucsc.edu/goldenPath/help/iupac.html')
    #sequence match score
    parser.add_argument('-s', dest='aln_scr', required=False, default=90,
                        help='Minimum allowable alignment score for positively\
                            identifying primer sequences (Default = 90)')
    args = parser.parse_args()

    if not os.path.exists(args.fa_in):
        print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),"Input file does\
              not exist: {}".format(args.fa_in))
        sys.exit(-1)

    main(args.fa_in, args.fa_out, args.f_primer, args.r_primer, args.aln_scr)