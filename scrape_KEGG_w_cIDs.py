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
__script_name__ = 'scrape_KEGG_w_cIDs.py'
__version__ = '0.0.2'
__profiling__ = 'False'

import os
import re
import sys
import argparse
import requests
from bs4 import BeautifulSoup
from datetime import datetime

###############################################################################
###############################################################################

## Function definitions

#compound names
def scrape_names_list(cmpd_soup):
    """
    Builds list of compound names associated with input cID
    """
    name_list = []
    #names listed in first <div style="width:555px;overflow-x:auto;overflow-y:hidden">" entry
    names = cmpd_soup.find_all("div", style="width:555px;overflow-x:auto;overflow-y:hidden")
    names_str = re.split(r'>|\n|;<br/|<br/', str(names))

    for entry in names_str:
        if bool(entry):
            if "/div" in entry: #first "/div" iteration ends compound list
                break
            elif "div" not in entry: #compound list starts after last "div" iteration
                name_list.append(entry)

    if len(name_list) == 0:
        name_list.append("NA")

    return name_list

#mapIDs, mIDs, and rIDs
def scrape_metab_lists(cmpd_soup):
    """
    Builds a list of mapIDs and mIDs associated with input cID
    """
    map_list = []
    mod_list = []
    rxn_list = []
    metabs = cmpd_soup.find_all('a') #mapIDs and mIDs in nobr entries

    for entry in metabs:
        line = str(entry)
        if "pathway" in line: #mapIDs always in format "<a href="/pathway/map####+C####">map####</a>"
            map_list.append(line.split(">")[1].strip("</a"))

        if "module" in line: #mIDs always in format "<a href="/module/m####">m####</a>"
            mod_list.append(line.split(">")[1].strip("</a"))
            
        if re.search(r'entry/R\d{5}', line): #rIDs and ECs both have 'entry'
            rxn_list.append(line.split(">")[1].strip("</a"))

    if len(map_list) == 0:
        map_list.append("NA")
    if len(mod_list) == 0:
        mod_list.append("NA")

    return map_list, mod_list, rxn_list

#kIDs
def scrape_kid_list(rxn_list, link):
    """
    Builds dictionary of rID entries and their list of associated kIDs
    """
    rxn_dict = {rid:[] for rid in rxn_list} #blank dictionary

    for rxn in rxn_list:
        kid_list = []

        #get soup
        rxn_html = requests.get(link + rxn)
        rxn_soup = BeautifulSoup(rxn_html.content, "html.parser")
        kids = rxn_soup.find_all('nobr') #kIDs in nobr entries

        for entry in kids:
            line = str(entry)
            if "entry" in line: #kIDs always in format "<a href="/entry/K####">K####</a>"
                kID = line.split(">")[2].strip("</a")
                if kID.startswith("K"):
                    kid_list.append(kID)

        if len(kid_list) == 0:
            kid_list.append("NA")

        rxn_dict[rxn] = kid_list
        print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"), rxn, "kIDs scraped")

    return rxn_dict

###############################################################################
###############################################################################

## Main script

def main(cmpd_tsv, out_tsv):
    """
    This is what we're here to do
    """
    #set default link
    link = "https://www.kegg.jp/entry/"

    #populate compound list
    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),"Parsing cID list")
    cmpd_list = []
    with open(cmpd_tsv, 'r') as cmpds:
        for entry in cmpds:
            cmpd_list.append(entry.strip('\n'))

    #set default compound dictionary
    cmpd_dict = {cid:[] for cid in cmpd_list}

    #iterate through cID list
    for compound in cmpd_list:
        cid_dict = {"Name":[],"mapIDs":[],"mIDs":[],"rIDs":[]} #blank dict for cID

        #get soup
        print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),"Gathering HTML content:", link + compound)
        html = requests.get(link + compound)
        soup = BeautifulSoup(html.content, "html.parser")

        #Name entries
        cid_dict["Name"] = scrape_names_list(soup)
        print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"), compound, "names scraped")

        #mapID, mID, and rID entries
        cid_dict["mapIDs"], cid_dict["mIDs"], rxns = scrape_metab_lists(soup)
        print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"), compound, "mapIDs scraped")
        print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"), compound, "mIDs scraped")
        print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"), compound, "rIDs scraped")
        
        if len(rxns) == 0:
            cid_dict["rIDs"] = {"NA":"NA"}
        else:
            cid_dict["rIDs"] = scrape_kid_list(rxns, link)

        #add cID entry to cmpd_dict
        cmpd_dict[compound] = cid_dict
        print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"), compound, "completed")

    #write results
    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"), "Writing results file")
    with open(out_tsv, 'w+') as file:
        for key in cmpd_dict.keys():
            file.write(">"+key+":\n")
            file.write("\t#Name:\t"+'\t'.join(cmpd_dict[key]['Name'])+"\n")
            file.write("\t#mapIDs:\t"+'\t'.join(cmpd_dict[key]['mapIDs'])+"\n")
            file.write("\t#mIDs:\t"+'\t'.join(cmpd_dict[key]['mIDs'])+"\n")
            file.write("\t#rIDs:\n")
            for rxn in cmpd_dict[key]['rIDs'].keys():
                file.write("\t\t##"+rxn+":\t"+'\t'.join(cmpd_dict[key]['rIDs'][rxn])+"\n")

    print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"), "Success! cID attributes writen to", out_tsv)

    pass

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='scrape_KEGG_w_cIDs.py')
    parser.add_argument('-i', dest='cmpd_tsv', required=True,
                        help='Input list of cIDs to search (tsv format, one cID per line)'),
    parser.add_argument('-o', dest='out_tsv', required=True,
                        help='Name of cID attribute file (tsv format)')
    args = parser.parse_args()

    if not os.path.exists(args.cmpd_tsv):
        print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),"Input list of cIDs doesnt exist: {}".format(args.cmpd_tsv))
        sys.exit(-1)

    if os.path.exists(args.out_tsv):
        print(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]:"),"Cannot overwrite existing cID attribute file: {}".format(args.out_tsv))

    main(args.cmpd_tsv, args.out_tsv)
