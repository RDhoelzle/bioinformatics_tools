import sys

domtbl = sys.argv[1]

H="#"
DOM="domains"
LEN="length"

def extract_domains(domtbl):
    
    domain_dict = {}

    
    for line in open(domtbl):
        sline = line.strip().split()
        if line.startswith(H): continue
        
        s_name = sline[0]
        p_name = sline[3]
        
        s_length = float(sline[2])
        p_length = float(sline[5])
        
        p_from = int(sline[15])
        s_from = int(sline[19])
        
        p_to = int(sline[16])
        s_to = int(sline[20])
        
        bit = float(sline[13])
        
        if(len(range(p_from, p_to)) / p_length >= 0.75):
            if s_name in domain_dict:
    
                for old_dom_name, entry in domain_dict[s_name][DOM].items():
                    old_span = range(entry[0], entry[1])
                    new_span = set(range(s_from, s_to))
                    if(len(new_span.intersection(old_span))!=0):
                        if bit > entry[2]:
                            del domain_dict[s_name][DOM][old_dom_name]
                            domain_dict[s_name][DOM][p_name] = [s_from, s_to, bit]
                            
                    else:
                        domain_dict[s_name][DOM][p_name] = [s_from, s_to, bit]
            else:
                domain_dict[s_name] = {DOM: {p_name: [s_from, s_to, bit]},
                                       LEN: s_length}
    
    return domain_dict



domains = extract_domains(domtbl)           


for sequence_name, domain_info in domains.items():
    seq_order = []
    seq_dom = {}
    out_order = []
    for key, item in domain_info[DOM].items():
        start = item[0]
        seq_order.append(start)
        seq_dom[start] = key
    for pos in sorted(seq_order):
        out_order.append(seq_dom[pos])
    print "%s\t%i\t%s" % (sequence_name, domain_info[LEN], ','.join(out_order))
