#!/usr/bin/env python
# https://github.com/shenwei356/bio_scripts/
'''
fetch taxon information by species name or taxid.

Take home message:

    1. using cache to avoid repeatly search
    2. object of Entrez.read(Entrez.efetch()) could be treated as list,
       but it could not be rightly pickled. Using Json is also not OK.
       The right way is cache the xml text.

           search = Entrez.efetch(id=taxid, db="taxonomy", retmode="xml")
           # data = Entrez.read(search)
           ##  read and parse xml
           data_xml = search.read()
           data = list(Entrez.parse(StringIO(data_xml)))
    3. pickle file was fragile. a flag file could be used to detect whether
       data is rightly dumped.
    4. using multi-threads to accelerate fetching.

'''

from __future__ import print_function
import sys
import argparse
import os
import re
import shutil
import pickle
from StringIO import StringIO
from multiprocessing import Pool
from Bio import Entrez

parser = argparse.ArgumentParser(
    description=
    "fetch taxon information by species name or taxid. Cache used to avoid repeatly search",
    epilog="https://github.com/shenwei356/bio_scripts/")

parser.add_argument('infile', help='species name/taxid list')
parser.add_argument('-n',
                    '--by-name',
                    action='store_true',
                    help='search by species name')
parser.add_argument('-t',
                    '--threads',
                    type=int,
                    default=4,
                    help='threads number, default:4')

default_cache_path = os.path.join(
    os.path.expanduser("~"), '.taxon', 'taxon_map.pickle')
parser.add_argument(
    '-c',
    '--cache-file',
    type=str,
    default=default_cache_path,
    help='taxon_map cache file, default: {}'.format(default_cache_path))
parser.add_argument('-d',
                    '--delete-cache-file',
                    action='store_true',
                    help='delete cache file')

args = parser.parse_args()

# ================[ caching feteched data ]==================
cache = dict()

# a flag file to check if the pickle file is ok, its existance means not ok
flag_file = '{}.close-by-accident'.format(args.cache_file)

if args.delete_cache_file:
    if os.path.exists(args.cache_file):
        os.unlink(args.cache_file)
    if os.path.exists(flag_file):
        os.unlink(flag_file)

# read cache if available
if os.path.exists(args.cache_file):
    sys.stderr.write('[INFO] read taxon_map cache from file: {}\n'.format(
        args.cache_file))

    if not os.path.exists(flag_file):
        cache = pickle.load(open(args.cache_file, 'rb'))
    else:
        sys.stderr.write(
            '[INFO] it seems that last run failed. delete cache file.\n')
        os.unlink(flag_file)
    # cache = pickle.load(open(args.cache_file, 'rb'))
else:
    sys.stderr.write('[INFO] create new taxon_map cache file: {}\n'.format(
        args.cache_file))

    cache_dir = os.path.dirname(args.cache_file)
    if not os.path.exists(cache_dir):
        os.mkdir(cache_dir)

cache_fh = open(args.cache_file, 'wb')

open(flag_file, 'w').close()


# ================[ fetching method ]==================
def get_tax_id(species):
    species = species.replace(" ", "+").strip()

    search = Entrez.esearch(term=species, db="taxonomy", retmode="xml")
    record = Entrez.read(search)

    return record['IdList'][0]


# ================[ fetching method ]==================
def get_tax_data(taxid):
    if not re.search('^\d+$', taxid):
        sys.stderr.write(
            '[ERROR] do you use species name as query? you may use flag: -n\n')
        if os.path.exists(flag_file):
            os.unlink(flag_file)
        sys.exit(0)
    search = Entrez.efetch(id=taxid, db="taxonomy", retmode="xml")

    # return Entrez.read(search) # if not using pickle, this is enough

    # save xml for pickle
    data_xml = search.read()
    return list(Entrez.parse(StringIO(data_xml))), data_xml


# ================[ fetching and outputing ]==================
def fetch_taxon(query):
    if query in cache:
        sys.stderr.write('[INFO] cached query: {}\n'.format(query))

        taxon_data_xml = cache[query]
        data = list(Entrez.parse(StringIO(taxon_data_xml)))
    else:
        sys.stderr.write('[INFO] new query: {}\n'.format(query))

        if args.by_name:
            taxid = get_tax_id(query)
            data, taxon_data_xml = get_tax_data(taxid)

            cache[taxid] = taxon_data_xml
        else:
            data, taxon_data_xml = get_tax_data(query)

            cache[data[0]['ScientificName']] = taxon_data_xml

        # save xml for pickle
        cache[query] = taxon_data_xml

    # output
    lineage = data[0]['Lineage']
    division = data[0]['Division']
    taxid = data[0]['TaxId']

    CommonName = ''
    if 'OtherNames' in data[0] and 'GenbankCommonName' in data[0][
            'OtherNames']:
        CommonName = data[0]['OtherNames']['GenbankCommonName']

    ScientificName = data[0]['ScientificName']

    if args.by_name:
        print('\t'.join([taxid, query, division, CommonName, lineage]))
    else:
        print('\t'.join([query, ScientificName, division, CommonName, lineage
                         ]))

# ================[ read query list ]==================
Entrez.email = "tmp@gmail.com"

species_list = list()
with open(args.infile) as fh:
    for species in fh:
        species = species.rstrip().lstrip()
        if len(species) == 0:
            continue
        species_list.append(species)

# ================[ fetching with multiprocessing ]==================
pool = Pool(args.threads)
#pool.map(fetch_taxon, species_list)
map(fetch_taxon, species_list)

# ================[ caching ]==================
pickle.dump(cache, cache_fh, -1)
if os.path.exists(flag_file):
    os.unlink(flag_file)
