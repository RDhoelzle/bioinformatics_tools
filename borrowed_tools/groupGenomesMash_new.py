#!/usr/bin/env python
###############################################################################
#
# __groupGenomes__.py 
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
__copyright__ = "Copyright 2016"
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

# local imports

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class Data(object):
    def __init__(self, name):
        self.__name  = name
        self.__links = set()

    @property
    def name(self):
        return self.__name

    @property
    def links(self):
        return set(self.__links)

    def add_link(self, other):
        self.__links.add(other)
        other.__links.add(self)

class Taxonomy(object):
    def __init__(self,
                 taxonomy_file,
                 outfile,
                 taxon_rank):
        # data
        self.taxonomy_data = {}
        self.genome_groups = {}
        
        # run scripts
        self.wrapper(taxonomy_file, outfile, taxon_rank)
    
    def wrapper(self,
                taxonomy_file,
                outfile,
                taxon_rank):
        
        # read in taxonomy file, store data
        self.parseTaxonomyFile(taxonomy_file)
        
        # group genomes by specific taxonomy
        self.groupGenomesByTaxonomy(outfile, taxon_rank)
        
        # write groups to file
        self.writeGroupsToOutfile(outfile)
        
    def parseTaxonomyFile(self,
                          taxonomy_file):
        with open(taxonomy_file) as fh:
            for l in fh:
                tabs = l.rstrip().split("\t")
                
                if 'clc_assembled.scaffolds.metabat-bins-' in l:
                    genome   = tabs[0]
                    taxonomy = tabs[138]
                    
                    self.taxonomy_data[genome] = taxonomy
    
    def groupGenomesByTaxonomy(self,
                               outfile,
                               taxon_rank):
        # set to lower case
        taxon_rank = taxon_rank.lower()
        
        # loop through genomes, grouping by taxonomy
        for genome in self.taxonomy_data.keys():
            taxonomy = self.taxonomy_data[genome]
            
            semi_colon = taxonomy.split(";")
            _phylum = semi_colon[1]
            _class  = semi_colon[2]
            _order  = semi_colon[3]
            _family = semi_colon[4]
            _genus  = semi_colon[5]
            _species= semi_colon[6]
            
            if taxon_rank == 'species':
                if len(_species) > 4:
                    
                    # add genome as is
                    self.addGenome(genome, taxonomy)
                    
                else:
                    # add uid to taxonomy
                    self.makeTaxonomyUnique(taxonomy, genome)
            
            elif taxon_rank == 'genus':
                if len(_genus) > 4:
                    
                    # add genome as is
                    self.addGenome(genome, taxonomy)
                    
                else:
                    # add uid to taxonomy
                    self.makeTaxonomyUnique(taxonomy, genome)
            
            elif taxon_rank == 'family':
                if len(_family) > 4:
                    
                    # add genome as is
                    self.addGenome(genome, taxonomy)
                    
                else:
                    # add uid to taxonomy
                    self.makeTaxonomyUnique(taxonomy, genome)
            
            elif taxon_rank == 'order':
                if len(_order) > 4:
                    
                    # add genome as is
                    self.addGenome(genome, taxonomy)
                    
                else:
                    # add uid to taxonomy
                    self.makeTaxonomyUnique(taxonomy, genome)
            
            elif taxon_rank == 'class':
                if len(_class) > 4:
                    
                    # add genome as is
                    self.addGenome(genome, taxonomy)
                    
                else:
                    # add uid to taxonomy
                    self.makeTaxonomyUnique(taxonomy, genome)
            
            elif taxon_rank == 'phylum':
                if len(_phylum) > 4:
                    
                    # add genome as is
                    self.addGenome(genome, taxonomy)
                    
                else:
                    # add uid to taxonomy
                    self.makeTaxonomyUnique(taxonomy, genome)
    
    def writeGroupsToOutfile(self,
                             outfile):
        # open outfile
        of = open(outfile, 'w')
        
        group_number = 1
        
        # print header
        of.write("#group_number\tcount\ttaxonomy\tgenomes\n")
        
        for taxonomy in self.genome_groups.keys():
            genome_concatenated = self.genome_groups[taxonomy][0]
            
            try:
                for genome in self.genome_groups[taxonomy][1:]:
                    genome_concatenated = "%s,%s" % (genome_concatenated,genome)
            except IndexError:
                pass
                
            str_to_write = "\t".join(['Group %d:' % group_number,
                                      '%s' % str(len(genome_concatenated.split(","))),
                                      taxonomy,
                                      genome_concatenated])
            
            of.write("%s\n" % str_to_write)
            
            # increment group number
            group_number +=1 
    
    def makeTaxonomyUnique(self,
                           taxonomy,
                           genome
                           ):
        
        for uid in range(1000):
            
            # add uid to taxonomy name
            new_taxonomy = "%s_uid_%d" % (taxonomy, uid)
            
            # check if unique
            if new_taxonomy in self.genome_groups:
                pass
            
            else:
                # add genome
                self.addGenome(genome,new_taxonomy)
                break
    
    def addGenome(self,
                  genome,
                  taxonomy):
        try:
            self.genome_groups[taxonomy].append(genome)
        except KeyError:
            self.genome_groups[taxonomy] = [genome]
    
class AAI(object):
    def __init__(self,
                 comparem_aai_file,
                 outfile,
                 threshold):
        
        self.aai_data = {}
        
        self.singletons = {}
        
        # run scripts
        self.wrapper(comparem_aai_file,
                     outfile,
                     threshold)
    
    def wrapper(self,
                comparem_aai_file,
                outfile,
                threshold):
        
        self. parseComparemFile(comparem_aai_file,
                                threshold)
        
        self.dereplicateGenomes(outfile)
        
    def parseComparemFile(self,
                          comparem_aai_file,
                          threshold):
        
        with open(comparem_aai_file) as fh:
            
            line_number = 0
            
            for l in fh:
                
                tabs = l.rstrip().split("\t")
                
                if line_number > 0:
                    
                    genome_a            = tabs[0]
                    genes_a             = tabs[1]
                    genome_b            = tabs[2]
                    genes_b             = tabs[3]
                    orthologous_genes   = tabs[4]
                    mean_aai            = float(tabs[5])
                    std_aai             = tabs[6]
                    orthologous_fraction= tabs[7]
                    
                    if mean_aai >= threshold:
                        
                        self.addConnection(genome_a, genome_b)
                        
                        # add genomes to non-singelton list
                        self.addGenome(genome_a)
                        self.addGenome(genome_b)
                    else:
                        self.removeGenome(genome_a)
                        self.removeGenome(genome_b)
                        
                line_number += 1
        
    def addGenome(self,
                  genome):
        self.singletons[genome] = 1
    
    def removeGenome(self,
                     genome):
        try:
            if self.singletons[genome] == 1:
                pass
        except KeyError:
            self.singletons[genome] = 0
    
    def dereplicateGenomes(self,
                           outfile):
        node_dict = {}

        nodes = set()

        # build trees from blast data

        regions_to_nodes= {}
        
        for genome_a in self.aai_data.keys():
            
            if genome_a not in regions_to_nodes:
                D = Data(genome_a)
                regions_to_nodes[genome_a] = D
                nodes |= {D}
                
                for genome_b in self.aai_data[genome_a].keys():
                    
                    if genome_b not in regions_to_nodes:
                        D = Data(genome_b)
                        regions_to_nodes[genome_b] = D
                        nodes |= {D}

        for genome_a in self.aai_data.keys():
            
            for genome_b in self.aai_data[genome_a].keys():
                regions_to_nodes[genome_a].add_link(regions_to_nodes[genome_b])

        # Find all the connected components
        number = 1

        # initialise outfile
        outfile = open(outfile, 'w')
        
        # print groups
        for components in self.connected_components(nodes):
            names = sorted(node.name for node in components)
            names_count= len(names)
            names = ",".join(names)
            #outfile.write("Group %i:\t%d\t%s\n" % (number,names_count, names))
            outfile.write("Group %i:\t%d\t%s\n" % (number,names_count, names))
            #print "Group %i:\t%d\t%s" % (number,names_count, names)
            number += 1        
    
        # print singletons
        for genome in self.singletons.keys():
            if self.singletons[genome] == 0:
                outfile.write("Group %i:\t1\t%s\n" % (number,genome))
                number += 1 
    
    def addConnection(self, id1, id2):
        try:
            self.aai_data[id1][id2] = 1 
        except KeyError:
            self.aai_data[id1] = {id2:1}
        try:
            self.aai_data[id2][id1] = 1 
        except KeyError:
            self.aai_data[id2] = {id1:1}
    
    def connected_components(self,nodes):
        # List of connected components found. The order is random.
        result = []

        # Make a copy of the set, so we can modify it.
        nodes = set(nodes)

        nn = float(len(nodes))
        remains = nn

        # Iterate while we still have nodes to process.
        while nodes:

            # Get a random node and remove it from the global set.
            n = nodes.pop()

            # This set will contain the next group of nodes connected to each other.
            group = {n}

            # Build a queue with this node in it.
            queue = [n]

            # Iterate the queue.
            # When it's empty, we finished visiting a group of connected nodes.
            while queue:

                # Consume the next item from the queue.
                n = queue.pop(0)

                # Fetch the neighbors.
                neighbors = n.links

                # Remove the neighbors we already visited.
                neighbors.difference_update(group)

                # Remove the remaining nodes from the global set.
                remains -= float(len(neighbors))
                #print "%f completed" % (1 - remains/nn)
                nodes.difference_update(neighbors)

                # Add them to the group of connected nodes.
                group.update(neighbors)

                # Add them to the queue, so we visit them in the next iterations.
                queue.extend(neighbors)

            # Add the group to the list of groups.
            result.append(group)

        # Return the list of groups.
        return result

class MASH(object):
    def __init__(self,
                 mash_file,
                 checkm_file,
                 outfile,
                 threshold):
        # data
        self.mash_data   = {}
        self.singletons  = {}
        self.checkm_data = {}
        self.dereplicated_groups = {}
        
        # run scripts
        self.parseCheckmFile(checkm_file)
        self.parseMashFile(mash_file,threshold)
        self.dereplicateGenomes(outfile)
    
    def parseCheckmFile(self,
                        checkm_file):
        with open(checkm_file) as fh:
            line_number = 0
            for l in fh:
                tabs = l.rstrip().split("\t")
                if line_number >0:
                    genome          = "%s.fa" % tabs[0]
                    completion      = float(tabs[11])
                    contamination   = float(tabs[12])
                    
                    # add to checkm data dictionary
                    self.checkm_data[genome] = [completion, contamination]
                line_number += 1
        
    def parseMashFile(self,
                      mash_file,
                      threshold):
        
        with open(mash_file) as fh:
            line_number = 0
            for l in fh:
                tabs = l.rstrip().split("\t")
                if line_number > 0:
                    genome_a            = tabs[0].split("/")[-1].replace('_genomic.fna','')
                    genome_b            = tabs[1].split("/")[-1].replace('_genomic.fna','')
                    mash_dist           = (1 - float(tabs[2]))*100
                    p_value             = float(tabs[3])
                    matching_hashes     = tabs[4]
                    
                    if mash_dist >= threshold:
                        self.addConnection(genome_a, genome_b)
                        # add genomes to non-singelton list
                        self.addGenome(genome_a)
                        self.addGenome(genome_b)
                    else:
                        self.removeGenome(genome_a)
                        self.removeGenome(genome_b)
                line_number += 1
        
    def addGenome(self,
                  genome):
        self.singletons[genome] = 1
    
    def removeGenome(self,
                     genome):
        try:
            if self.singletons[genome] == 1:
                pass
        except KeyError:
            self.singletons[genome] = 0
    
    def dereplicateGenomes(self,
                           outfile):
        node_dict = {}

        nodes = set()

        # build trees from blast data

        regions_to_nodes= {}
        
        for genome_a in self.mash_data.keys():
            
            if genome_a not in regions_to_nodes:
                D = Data(genome_a)
                regions_to_nodes[genome_a] = D
                nodes |= {D}
                
                for genome_b in self.mash_data[genome_a].keys():
                    
                    if genome_b not in regions_to_nodes:
                        D = Data(genome_b)
                        regions_to_nodes[genome_b] = D
                        nodes |= {D}

        for genome_a in self.mash_data.keys():
            
            for genome_b in self.mash_data[genome_a].keys():
                regions_to_nodes[genome_a].add_link(regions_to_nodes[genome_b])

        # Find all the connected components
        number = 1

        # initialise outfile
        outfile = open(outfile, 'w')
        
        # print groups
        for components in self.connected_components(nodes):
            names = sorted(node.name for node in components)
            names_count= len(names)
            names = ",".join(names)
            genomes = names
            
            rep_genome = self.pickRepresentativeGenome(genomes, number)
            
            genomes_split = names.split(",")
            
            genomes_concatenated = []
            for genome in genomes_split:
                if genome != rep_genome:
                    genomes_concatenated.append(genome)
            
            #outfile.write("Group %i:\t%d\t%s\n" % (number,names_count, names))
            outfile.write("%s\t%s\n" % (rep_genome, ",".join(genomes_concatenated)))
            #print "Group %i:\t%d\t%s" % (number,names_count, names)
            number += 1        
    
        # print singletons
        for genome in self.singletons.keys():
            if self.singletons[genome] == 0:
                outfile.write("%s\t\n" % (genome))
    
    def pickRepresentativeGenome(self,
                                 genomes,
                                 group):
        genomes = genomes.split(",")
        
        for genome in genomes:
            if group not in self.dereplicated_groups:
                
                # add to dictionary
                self.dereplicated_groups[group] = genome
            
            else:
                # incumbent genome for group #
                genome_incumbent               = self.dereplicated_groups[group]
                genome_incumbent_contamination = self.checkm_data[genome_incumbent][1]
                genome_incumbent_completion    = self.checkm_data[genome_incumbent][0]  
                genome_incumbent_relative_completion = genome_incumbent_completion - genome_incumbent_contamination
                
                # new genome to compare to incumbent
                genome_contamination           = self.checkm_data[genome][1]
                genome_completion              = self.checkm_data[genome][0] 
                genome_relative_completion = genome_completion - genome_contamination
                
                if genome_relative_completion > genome_incumbent_relative_completion:
                    # replace incumbent with new genome
                    self.dereplicated_groups[group] = genome
        
        return self.dereplicated_groups[group]
        
    def addConnection(self, id1, id2):
        try:
            self.mash_data[id1][id2] = 1 
        except KeyError:
            self.mash_data[id1] = {id2:1}
        try:
            self.mash_data[id2][id1] = 1 
        except KeyError:
            self.mash_data[id2] = {id1:1}
    
    def connected_components(self,nodes):
        # List of connected components found. The order is random.
        result = []

        # Make a copy of the set, so we can modify it.
        nodes = set(nodes)

        nn = float(len(nodes))
        remains = nn

        # Iterate while we still have nodes to process.
        while nodes:

            # Get a random node and remove it from the global set.
            n = nodes.pop()

            # This set will contain the next group of nodes connected to each other.
            group = {n}

            # Build a queue with this node in it.
            queue = [n]

            # Iterate the queue.
            # When it's empty, we finished visiting a group of connected nodes.
            while queue:

                # Consume the next item from the queue.
                n = queue.pop(0)

                # Fetch the neighbors.
                neighbors = n.links

                # Remove the neighbors we already visited.
                neighbors.difference_update(group)

                # Remove the remaining nodes from the global set.
                remains -= float(len(neighbors))
                #print "%f completed" % (1 - remains/nn)
                nodes.difference_update(neighbors)

                # Add them to the group of connected nodes.
                group.update(neighbors)

                # Add them to the queue, so we visit them in the next iterations.
                queue.extend(neighbors)

            # Add the group to the list of groups.
            result.append(group)

        # Return the list of groups.
        return result

###############################################################################
###############################################################################
###############################################################################
###############################################################################
   
def doWork( args ):
    """ Main wrapper"""
    
    # run scripts
    if args.subparser_name == 'aai':
        AAI(args.comparem_aai_file,
            args.outfile,
            args.threshold)
    elif args.subparser_name == 'taxonomy':
        Taxonomy(args.taxonomy_file,
                 args.outfile,
                 args.taxon_rank)
    elif args.subparser_name == 'mash':
        MASH(args.mash_file,
             args.checkm_file,
             args.outfile,
             args.threshold)
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help="--", dest='subparser_name')
    
    #------------------------------
    
    aai_parser = subparsers.add_parser('aai',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                        help='Groups genome by aai',
                                        description='Groups genome by aai')
    aai_parser.add_argument('comparem_aai_file',help='File containing all pairwise compareM results ')
    aai_parser.add_argument('outfile',help='Output filename')
    aai_parser.add_argument('-threshold','--threshold',type=int,default=95,help='Set AAI threshold for clustering')
    
    #------------------------------
    
    mash_parser = subparsers.add_parser('mash',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                        help='Groups genome by mash',
                                        description='Groups genome by mash')
    mash_parser.add_argument('mash_file',help='File containing all pairwise mash results ')
    mash_parser.add_argument('checkm_file',help='Checkm tab-delimited file')
    mash_parser.add_argument('outfile',help='Output filename')
    mash_parser.add_argument('-threshold','--threshold',type=int,default=95,help='Set mash percentage threshold for clustering')
    
    #------------------------------
    
    aai_parser = subparsers.add_parser('taxonomy',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                        help='Groups genome by taxonomy',
                                        description='Groups genome by taxonomy')
    aai_parser.add_argument('taxonomy_file',help='File containing all genomes and their corresponding taxonomy')
    aai_parser.add_argument('outfile',help='Output filename')
    aai_parser.add_argument('-tr','--taxon_rank',default='Species',help='Group genome by specified taxonomy rank: Phylum, Class, Order, Family, Genus, Species (default)')
    
    # parse the arguments
    args = parser.parse_args()
    
    # do what we came here to do
    doWork(args)
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################