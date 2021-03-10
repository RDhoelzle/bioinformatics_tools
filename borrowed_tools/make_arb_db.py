from Bio import SeqIO
import IPython
import sys

records_path = sys.argv[1]
domains_path = sys.argv[2]
output_arbdb_path = sys.argv[3]
#output_gtdb_tax_path = sys.argv[3]
#output_ncbi_tax_path = sys.argv[4]

empty_tax_string = 'd__;p__;c__;o__;f__;g__;s__'
header = ['name', 'sequence_source', 'sequence_length', 'sequence_domains', 'sequence_ncbi_tax', 'sequence_gtdb_tax', 'annotation', 'phyla', 'sequence']
ncbi_taxonomy_path = '/srv/db/gtdb_reference_databases/gtdb_release75/taxonomy/ncbi_taxonomy.2016_04_27.tsv'
gtdb_taxonomy_path = '/srv/db/gtdb_reference_databases/gtdb_release75/taxonomy/gtdb_taxonomy.2016_04_27.tsv'

ncbi_taxonomy_dict = {x.strip().split('\t')[0]:x.strip().split('\t')[1]
                        for x in open(ncbi_taxonomy_path)}
gtdb_taxonomy_dict = {x.strip().split('\t')[0]:x.strip().split('\t')[1]
                        for x in open(gtdb_taxonomy_path)}

records_dict = SeqIO.to_dict(SeqIO.parse(records_path, "fasta"))

pfam_domain_dict = {}
for line in open(domains_path):
    sline = line.strip().split('\t')
    pfam_domain_dict[sline[0]] = sline[2]

output_dictionary = {}

for record_name, record_object in records_dict.items():
    
    record_entry = []
    
    if record_name.startswith("contig"):
        record_source = 'gtdb'
        record_ncbi_tax = ncbi_taxonomy_dict[record_name.split('~')[1]]
        record_gtdb_tax = gtdb_taxonomy_dict[record_name.split('~')[1]]
    else:
        record_source = 'SQISSprot'
        record_ncbi_tax = empty_tax_string
        record_gtdb_tax = empty_tax_string

    record_sequence_length = str(len(record_object.seq))
    
    try:
        record_domains = pfam_domain_dict[record_name]
    except:
        record_domains = 'none'
        
    record_entry.append(record_sequence_length)
    record_entry.append(record_domains)

    phyla = record_ncbi_tax.split(';')[1]

        
    annotation = '_'.join(record_object.description.split()[1:]).split('_[')[0]
    output_dictionary[record_name] = [record_source,
                                      record_sequence_length,
                                      record_domains,
                                      record_ncbi_tax,
                                      record_gtdb_tax,
                                      annotation,
                                      phyla,
                                      str(record_object.seq)]

# Write arb database
with open(output_arbdb_path, 'w') as out_io:
    out_io.write('\t'.join(header) + '\n')
    for record_name, entry in output_dictionary.items():
        output_line = [record_name] + entry

        out_io.write('\t'.join(output_line) + '\n')

            # Write tax
