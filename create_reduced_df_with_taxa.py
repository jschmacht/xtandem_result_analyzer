from pathlib import Path
import argparse
import pandas as pd
from create_PSM_df import PSM_FDR
from ReadAccTaxon import ReadAccTaxon
from collections import defaultdict
from handling_acc_files import HelperMethod


def flatten_list(l):
    return [item for sublist in l for item in sublist]


def read_crap(file):
    crap = set()
    with open(file, "r") as f:
        for line in f:
            if line.startswith('>'):
                if '|' in line:
                    fields = [item for item in line.split('|')]
                    crap.add(fields[1])
                else:
                    crap.add(line[1:].strip())
    return crap

def get_accs_from_df(x_tandem_result_tsv, db_type, decoy_tag):
    acc_from_tsv = pd.read_csv(str(x_tandem_result_tsv), delimiter='\t')['Protein'].tolist()
    accs = get_accs_from_accs(acc_from_tsv, db_type, decoy_tag)
    return accs

def get_accs_from_accs(accs, db_type, decoy_tag):
    processed_accs = set()
    for acc in accs:
        if db_type == 'uniprot' or db_type == 'swissprot':
            try:
                if decoy_tag not in acc:
                    processed_accs.add(acc.split()[0].split('|')[1])
            # e.g. CRAP: KKA1_ECOLX
            except IndexError:
                processed_accs.add(acc.split()[0])
        elif db_type == 'ncbi':
            # ncbi: 'generic|AC:5988671|WP_080217951.1', 'generic|AC:2500558_REVERSED|EQV03804.1 ERA83742.1 WP_021536075.1-REVERSED'
            # problem: generic|AC:373747|pir||D85980 WP_000133047.1 AAG58304.1 -> maxsplit=2
            # CRAP accs starts with 'sp|' or Index Error
            try:
                # CRAP accs = sp|acc
                if decoy_tag not in acc and not acc.startswith('sp|'):
                    processed_accs.add(acc.split()[0].split('|', maxsplit=2)[2])
            except IndexError:
                accs.add(acc.split()[0])
        elif db_type == 'custom':
            if decoy_tag not in acc:
                processed_accs.add(acc.strip())
    return processed_accs


def get_acc2taxon_dict(db_path, db_type, accs=None):
    acc2tax_reader=ReadAccTaxon(db_path, db_type)
    if db_type == 'custom':
        acc_2_taxon_dict = acc2tax_reader.get_acc2taxonID_dict(db_path/'acc2tax_custom')
    elif db_type == 'uniprot' or db_type == 'swissprot':
        acc_2_taxon_dict = acc2tax_reader.read_acc2tax(accs)
    return acc_2_taxon_dict


def get_multispecies_accs(path_to_multispecies_acc, accs):
    """
    :param path_to_multispecies_acc: path to multispecies file
    :param accs: set of ncbi accs
    :return: multi_acc to list of accs dict
        """
    multiacc2accs_dict={}
    with open(path_to_multispecies_acc, 'r') as input:
        for line in input:
            fields = line.split()
            l = [acc in accs for acc in fields]
            if any(l):
                pos_of_acc = [i for i, x in enumerate(l) if x][0]
                multiacc2accs_dict[fields[pos_of_acc]] = fields[0:]
    return multiacc2accs_dict


def get_tsv_multspecies_acc_to_taxa_dict(multiacc2acc_dict, acc_2_taxon_dict):
    tsv_multi_acc_to_taxon_dict = defaultdict(set)
    for acc, acc_list in multiacc2acc_dict.items():
        for multi_acc in acc_list:
            if multi_acc in acc_2_taxon_dict.keys():
                tsv_multi_acc_to_taxon_dict[multi_acc].add(acc_2_taxon_dict[multi_acc])
    return tsv_multi_acc_to_taxon_dict


def remove_accs_with_unsupported_taxa(multi_acc_2_taxon_dict, taxa):
    accs_with_not_matching_taxa = []
    for acc, taxon in multi_acc_2_taxon_dict.items():
        if int(taxon) not in taxa:
            accs_with_not_matching_taxa.append(acc)
    for acc in accs_with_not_matching_taxa:
        del multi_acc_2_taxon_dict[acc]
    return multi_acc_2_taxon_dict


def get_ncbi_multiacc_to_accs_dict(ncbi_accs_from_file, db_path, db_type):
    path_to_multiaccs = '/home/jules/Documents/databases/databases_tax2proteome/multispecies_acc'
    # key = first acc in result tsv, value: all accs
    multacc_reader = ReadAccTaxon(db_path, db_type, path_to_multiaccs)
    multiacc2acc_dict = multacc_reader.read_multispecies_accs(ncbi_accs_from_file)
    return multiacc2acc_dict


def get_taxa_from_acc2taxid_files(db_path, db_type, all_accs, taxa=None):
    acc2tax_reader = ReadAccTaxon(db_path, db_type)
    # only acc, taxa in target taxa (multiaccs file bring a lot unspecific taxa
    acc_2_taxon_dict = acc2tax_reader.read_acc2tax(all_accs, taxa)
    return acc_2_taxon_dict


def create_acc_from_file_2_taxa_set_dict(ncbi_accs_from_file, multiacc2acc_dict, acc_2_taxon_dict):
    acc_in_tsv_2_taxa_set_dict = defaultdict(set)
    for acc in ncbi_accs_from_file:
        if acc in multiacc2acc_dict.keys():
            for m_acc in multiacc2acc_dict[acc]:
                try:
                    acc_in_tsv_2_taxa_set_dict[acc].add(int(acc_2_taxon_dict[m_acc]))
                    del acc_2_taxon_dict[m_acc]
                # acc not in acc_2_taxon_dict, because taxon was not in taxa_all_level
                except KeyError:
                    continue
        else:
            acc_in_tsv_2_taxa_set_dict[acc].add(int(acc_2_taxon_dict[acc]))
            del acc_2_taxon_dict[acc]
    return acc_in_tsv_2_taxa_set_dict


def get_ncbi_acc2taxon_dict(ncbi_accs_from_file, db_path, db_type, taxa=None):
    """
    :param x_tandem_result_tsv:
    :param db_path:
    :param db_type:
    :param taxa: list of all relevant taxa, all taxa in taxon graph below specified level
    :return: acc_in_tsv_2_taxa_set_dict dict {acc string: {set of taxa(int)}
    """
    # for ncbi find all multispecies accessions with file ncbi acc
    # key = Multiacc or accs assigned to multiaccs
    multiacc2acc_dict = get_ncbi_multiacc_to_accs_dict(ncbi_accs_from_file, db_path, db_type)
    all_accs = set(flatten_list(multiacc2acc_dict.values())).union(ncbi_accs_from_file)
    # acc_2_taxon_dict contains 0 values for some accessions (are in prot.acc...)
    acc_2_taxon_dict = get_taxa_from_acc2taxid_files(db_path, db_type, all_accs, taxa)
    all_accs.clear()
    # acc_in_tsv_2_taxa_set_dict contains in taxa_set items like '345' and '0'
    acc_in_tsv_2_taxa_set_dict = create_acc_from_file_2_taxa_set_dict(ncbi_accs_from_file, multiacc2acc_dict, acc_2_taxon_dict)
    multiacc2acc_dict.clear()
    return acc_in_tsv_2_taxa_set_dict


def main():
    """
    1. based on db_type: get accs from tsv (split acc string)
    2. get acc_2_taxon_dict, NCBI only if taxon of acc2taxon in taxa, added to dict -> change
    acc_2_taxon_dict = get_ncbi_acc2taxon_dict(accs, db_path, db_type, all_taxa_of_level_set): removed all_taxa_of_level_set
    3. create_PSM_dataframe: sorting by Hyperscore, all matches to one spectrum with same score in same row, determine decoy hits
         add taxID column based on acc2tax dict, taxa_level column based on taxID column
    so given taxIDs are never used (except for NCBI, but I removed this
    :return: write reduced_tsv with all sorted and relevant information
    """
    parser = argparse.ArgumentParser(description='Read xtandem output .tsv')
    parser.add_argument('-i', '--input', dest='input', default=None, help='xtandem tsv with columns:'
                                                                          '[spectra_file, Spectrum ID, Peptide, Protein Accession + Description, Hyperscore, Evalue')
    parser.add_argument('-p', '--path', dest='path', default=None, help='Path to folder with acc2tax files. '
                                                                        '(Uniprot: acc2tax_uniprot, NCBI: pdb/prot.accession2taxid, '
                                                                        'Custom: acc2tax_custom')
    parser.add_argument('-c', '--crap', dest='crap', default=None, help='Crap-File')
    parser.add_argument('-g', '--tax_graph', dest='tax_graph', help='Path to taxdump.tar.gz')
    parser.add_argument('-l', '--level', dest='level', choices=['subspecies', 'species', 'genus', 'family', 'order', 'superkingdom'],
                        help='Level of database')
    parser.add_argument('-d', '--database', dest='database', choices=['ncbi', 'uniprot', 'custom', 'swissprot'], default='uniprot',
                        help='Database format.')
    parser.add_argument('-s', '--taxonset', dest='taxonset', choices=['tanca', 'kleiner'], default=None,
                        help='Taxon dataset used for analysis.')
    parser.add_argument('-t', '--taxon', dest='taxon', nargs='+', action='append',
                        help='NCBI taxon ID/s for database extraction. Multiple taxonIDs seperated by space.')
    parser.add_argument('-z', '--add_db', dest='add_db', choices=['ncbi', 'uniprot', 'custom', 'swissprot'], default=None,
                        help='if databases are of mixed origin, only working for custom_db + swissprot/uniprot or ncbi entries')
    # parser.add_argument('-f', '--fdr', dest='fdr', type=float, default=0.01, help='FDR-rate, default  = 0.01')
    parser.add_argument('-y', '--decoy', dest='decoy', default='REVERSED', help='Decoy_tag.')
    parser.add_argument('-x', '--threads', dest='threads', type=int, action="store", help='Number of threads.')

    options = parser.parse_args()
    # not in Kleiner DB: 536: Chromobacterium violaceum, 1407502: Stenotrophomonas maltophilia SeITE02
    # in brackets: taxon ID used in bachelor thesis:
    # 83333 Escherichia coli K-12 (511145, Escherichia coli str. K-12 substr. MG1655), 1294143 : Paracoccus denitrificans (1302247)
    # 294 Pseudomonas fluorescens (1114970: Pseudomonas fluorescens F113), 216596 (1041145: Rhizobium leguminosarum bv. viciae VF39)
    # 1280 (93061: Staphylococcus aureus subsp. aureus NCTC 8325, 46170),
    # 0 = dummy ID unknown seq
    # new taxa based on Kleiner DB
    Kleiner_taxIDs_0 = [262724, 882, 176299, 228410, 44577, 926571, 323848, 12022, 1283336, 10754, 101570, 224308, 216596,
                      1004788, 266265, 266264, 99287, 1294143, 1149133, 3055, 1280, 1977402, 294, 83333,
                      536, 1407502]
    # taxa bachelor thesis with woring
    Kleiner_taxIDs_1 = [262724, 882, 176299, 228410, 44577, 926571, 323848, 12022, 1283336, 10754, 101570, 224308, 1041145,
                      1004788, 266265, 266264, 99287, 1302247, 1149133, 3055, 93061, 1977402, 1114970, 511145,
                      536, 1407502, 1294143, 1619948]
    # taxa bachelor thesis
    # 46170, 216596 not in Kleiner_taxIDs_1, not in here: 1619948(Pseudomonas sp 21)
    # no difference for result
    Kleiner_taxIDs = [536, 882, 44577, 228410, 323848, 46170, 93061, 224308, 99287, 511145, 176299, 216596, 1041145,
                      262724, 266264, 266265, 1004788, 1114970, 1149133, 1294143, 1407502, 3055, 1302247, 926571,
                      # virus
                      10754, 101570, 1283336, 12022, 1977402]


    Tanca_taxIDs = [747, 5535, 655183, 1579, 1255, 4932, 1465, 1351, 562]
    if options.taxonset == 'kleiner':
        taxonIDs = Kleiner_taxIDs
        levels = ['subspecies', 'species', 'genus', 'family', 'order', 'superkingdom']
    elif options.taxonset == 'tanca':
        taxonIDs = Tanca_taxIDs
        levels = ['species', 'genus', 'family', 'order', 'superkingdom']
    if options.taxon:
        taxonIDs.extend([taxID for taxonlist in options.taxon for taxID in taxonlist])

    db_path = Path(options.path)
    db_type = options.database
    path_to_x_tandem_result_tsv = Path(options.input)
    path_to_crap = Path(options.crap)
    decoy_tag = options.decoy
    taxon_graph = HelperMethod.load_taxa_graph(Path(options.tax_graph))

    # acc_2_taxon_dict for identification file identified accs
    if db_type == 'custom':
        accs = set()
    else:
        accs = get_accs_from_df(path_to_x_tandem_result_tsv, db_type, decoy_tag)

    if db_type == 'ncbi':
        all_taxa_of_level_set = set(flatten_list(taxon_graph.get_all_taxids(taxonIDs, options.level)))
        acc_2_taxon_dict = get_ncbi_acc2taxon_dict(accs, db_path, db_type, all_taxa_of_level_set)
        # acc_in_tsv_2_taxa_set_dict
        # acc_2_taxon_dict = get_ncbi_acc2taxon_dict(accs, db_path, db_type)
    else:
        acc_2_taxon_dict = get_acc2taxon_dict(db_path, db_type, accs)
    accs.clear()

    if options.add_db:
        # accs= full accs without split, if custom: get dict of all accs in db
        accs = get_accs_from_df(path_to_x_tandem_result_tsv, 'custom', decoy_tag)
        add_accs = {acc for acc in accs if acc not in acc_2_taxon_dict.keys()}
        processed_add_acc = get_accs_from_accs(add_accs, options.add_db, decoy_tag)
        processed_add_acc_2_taxon_dict = get_acc2taxon_dict(db_path, options.add_db, processed_add_acc)
        add_acc_2_taxon_dict = {}
        for processed_acc, taxon in processed_add_acc_2_taxon_dict.items():
            for acc in add_accs:
                if processed_acc in acc:
                    add_acc_2_taxon_dict[acc] = taxon
        acc_2_taxon_dict.update(add_acc_2_taxon_dict)
    for k, v in acc_2_taxon_dict.items():
        if v == 0 or v=='0':
            print('error in taxonGraph, taxon = root, accs:', k,v)

    psm = PSM_FDR(path_to_x_tandem_result_tsv, path_to_crap, decoy_tag)
    reduced_df = psm.create_PSM_dataframe(db_type, options.level, taxon_graph, acc_2_taxon_dict)
    print(f"writing data frame to {path_to_x_tandem_result_tsv.parent.joinpath(path_to_x_tandem_result_tsv.stem + '_new_reduced.tsv')}... ")
    reduced_df.to_csv(str(path_to_x_tandem_result_tsv.parent.joinpath(path_to_x_tandem_result_tsv.stem + '_new_reduced.tsv')), sep='\t')


if __name__ == '__main__':
    main()
