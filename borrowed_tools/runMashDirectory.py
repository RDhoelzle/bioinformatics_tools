#!/usr/bin/env python
###############################################################################
#
# __runMashDirectory__.py 
#
###############################################################################
# #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or #
# (at your option) any later version. #
# #
# This program is distributed in the hope that it will be useful, #
# but WITHOUT ANY WARRANTY; without even the implied warranty of #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the #
# GNU General Public License for more details. #
# #
# You should have received a copy of the GNU General Public License #
# along with this program. If not, see <http://www.gnu.org/licenses/>. #
# #
###############################################################################

__author__ = "Joshua Daly"
__copyright__ = "Copyright 2017"
__credits__ = ["Joshua Daly"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "Joshua Daly"
__email__ = "joshua.daly@uqconnect.edu.au"
__status__ = "Development"

###############################################################################
# system imports
import argparse
import sys
from multiprocessing import Pool
from subprocess import Popen, PIPE
import subprocess
import re
import os
import fnmatch

# local imports

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class Mash(object):
    def __init__(self,
                 directory,
                 extension,
                 outfile,
                 threads):
        
        self.find = FIND()
        self.runMash(directory, extension, outfile, threads)
    
    def runMash(self,
                directory,
                extension,
                outfile,
                threads):
        # find genome fasta files
        genome_files = self.find.find('*%s' % extension, directory)
        cmds = []
        
        for i in range(len(genome_files)):
            for j in range(i+1, len(genome_files)):
                
                # pairwise distance calculations using mash
                cmds.append("mash dist %s %s >> %s" % (genome_files[i],
                                                       genome_files[j],
                                                       outfile))
        pool = Pool(processes=threads)
        print pool.map(runCommand, cmds)
        
class FIND(object):
    def find(self,
             pattern,
             path):
        result = []
        
        for root, dirs, files in os.walk(path):
            for name in files:
                if fnmatch.fnmatch(name, pattern):
                    result.append(os.path.join(root, name))
        return result

###############################################################################
###############################################################################
###############################################################################
###############################################################################

def runCommand(cmd):
    """Run a command and take care of stdout

expects 'cmd' to be a string like "foo -b ar"

returns (stdout, stderr)
"""
    args = cmd
    print cmd
    p = subprocess.Popen(args, shell=True, stdout=None, stderr=None) # shell=bash is not recommended. Only use when '>' must be in cmd.
    return p.communicate()

def doWork( args ):
    """ Main wrapper"""
    # run scripts
    Mash(args.directory,
         args.extension,
         args.outfile,
         args.threads)
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='PROG',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    #------------------------------
    
    parser.add_argument('directory',help='Directory containing genome fasta files')
    parser.add_argument('-x','--extension',default='fa',help='File extension')
    parser.add_argument('-o','--outfile',default='mash.txt',help='File extension')
    parser.add_argument('-t','--threads',default=5, type=int,help='Set number of processes to run')
    
    # parse the arguments
    args = parser.parse_args()
    
    # do what we came here to do
    doWork(args)
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################
