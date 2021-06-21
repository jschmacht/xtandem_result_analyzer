from pathlib import Path
import argparse
import pandas as pd
from PSM_FDR import PSM_FDR
from AccessionSearcher import AccessionSearcherNCBI
from ReadAccTaxon import ReadAccTaxon
from SearchAccessions import SearchAccessions
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


def get_acc2taxon_dict(x_tandem_result_tsv, db_path, db_type):
    acc2tax_reader=ReadAccTaxon(db_path, db_type)
    if db_type == 'custom':
        acc_2_taxon_dict = acc2tax_reader.get_acc2taxonID_dict(db_path/'acc2tax_custom')
    elif db_type == 'uniprot':
        # taxids_level = taxon_graph.get_all_taxids(taxonIDs, options.level)
        accs = {acc.split()[0].split('|')[1] for acc in pd.read_csv(str(x_tandem_result_tsv), delimiter='\t')['Protein'].tolist()}
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


def get_ncbi_acc2taxon_dict(x_tandem_result_tsv, db_path, db_type, taxa, decoy_tag):
    """
    :param x_tandem_result_tsv:
    :param db_path:
    :param db_type:
    :param taxa: list of all relevant taxa, all taxa in taxon graph below specified level
    :return:
    """
    # Protein: generic|AC:2905764_REVERSED|VVN60222.1-REVERSED, generic|AC:4682148|PMZ83684.1 PMZ91489.1 WP_054593944.1 PMZ34596.1 ALI00421.1 PMZ37242.1...
    # generic|AC:3582431_REVERSED|WP_032249104.1 KDU12656.1 KDT75424.1-REVERSED
    protein_accs_from_result_tsv = pd.read_csv(str(x_tandem_result_tsv), delimiter='\t')['Protein'].tolist()
    # problem: generic|AC:373747|pir||D85980 WP_000133047.1 AAG58304.1 -> maxsplit=2
    ncbi_accs_from_file = [acc.split('|', maxsplit=2)[2] for acc in protein_accs_from_result_tsv if decoy_tag not in acc]
    # for ncbi find all multispecies accessions with file ncbi acc
    path_to_multiaccs = '/home/jules/Documents/databases/databases_tax2proteome/multispecies_acc'
    multiaccs = set()
    single_accs = set()
    for i, acc in enumerate(ncbi_accs_from_file):
        if len(acc.split()) > 1:
            multiaccs.add(acc.split()[0].strip())
        else:
            single_accs.add(acc.strip())

    # key = first acc in result tsv, value: all accs
    multacc_reader = ReadAccTaxon(db_path, db_type, path_to_multiaccs)
    multiacc2acc_dict = multacc_reader.read_multispecies_accs(multiaccs)
    # key = Multiacc or accs assigned to multiaccs
    all_multiaccs = set(flatten_list(multiacc2acc_dict.values()))
    all_accs = all_multiaccs.union(single_accs)
    acc2tax_reader = ReadAccTaxon(db_path, db_type)
    acc_2_taxon_dict = acc2tax_reader.read_acc2tax(all_accs, taxa)
    print('acc_2_taxon_dict 1: ', len(acc_2_taxon_dict))
    acc_in_tsv_2_taxa_set_dict = defaultdict(set)
   # acc_2_taxon_dict = remove_accs_with_unsupported_taxa(acc_2_taxon_dict, taxa)
   # print('acc_2_taxon_dict with removed unsupported taxa: ', len(acc_2_taxon_dict))
    # taxon to set(taxon), in tsv_multi_acc_to_taxa_dict multiple taxa in set possible
    while single_accs:
        acc = single_accs.pop()
        try:
            acc_in_tsv_2_taxa_set_dict[acc].add(acc_2_taxon_dict[acc])
        # only from CRAP: sp|K22E_HUMAN| -> acc = ''
        except KeyError:
            continue
        del acc_2_taxon_dict[acc]
    while multiaccs:
        acc = multiaccs.pop()
        for multiacc in multiacc2acc_dict[acc]:
            try:
                acc_in_tsv_2_taxa_set_dict[acc].add(acc_2_taxon_dict[multiacc])
                del acc_2_taxon_dict[multiacc]
            except KeyError:
                continue
        del multiacc2acc_dict[acc]
  #  acc_in_tsv_2_taxon_set_dict = {acc:{taxon} for (acc,taxon) in acc_2_taxon_dict.items() if acc in single_accs}
  #  print('acc_in_tsv_2_taxon_dict ready1', len(acc_2_taxon_dict))
  #  tsv_multi_acc_to_taxa_dict = get_tsv_multspecies_acc_to_taxa_dict(multiacc2acc_dict, acc_2_taxon_dict)
    # single_accs.extend(flatten_list(multiacc2acc_dict.values()))
   # print('tsv_multi_acc_to_taxon_dict: ', len(tsv_multi_acc_to_taxa_dict))
   # acc_in_tsv_2_taxon_set_dict.update(tsv_multi_acc_to_taxa_dict)
    print('acc2tax_dict ready2', len(acc_in_tsv_2_taxa_set_dict))
    return acc_in_tsv_2_taxa_set_dict


def main():
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
    parser.add_argument('-d', '--database', dest='database', choices=['ncbi', 'uniprot', 'custom'], default='uniprot',
                        help='Database format.')
    parser.add_argument('-s', '--taxonset', dest='taxonset', choices=['tanca', 'kleiner'], default=None,
                        help='Taxon dataset used for analysis.')
    parser.add_argument('-f', '--fdr', dest='fdr', type=float, default=0.01, help='FDR-rate, default  = 0.01')
    parser.add_argument('-y', '--decoy', dest='decoy', default='REVERSED', help='Decoy_tag.')
    parser.add_argument('-x', '--threads', dest='threads', type=int, action="store", help='Number of threads.')

    options = parser.parse_args()
    # not in Kleiner DB: 536: Chromobacterium violaceum, 1407502: Stenotrophomonas maltophilia SeITE02
    # in brackets: taxon ID used in bachelor thesis:
    # 83333 Escherichia coli K-12 (511145, Escherichia coli str. K-12 substr. MG1655), 1294143 : Paracoccus denitrificans (1302247)
    # 294 Pseudomonas fluorescens (1114970: Pseudomonas fluorescens F113), 216596 (1041145: Rhizobium leguminosarum bv. viciae VF39)
    # 1280 (93061: Staphylococcus aureus subsp. aureus NCTC 8325),
    # 0 = dummy ID unknown seq
    # new taxa based on Kleiner DB
    Kleiner_taxIDs = [262724, 882, 176299, 228410, 44577, 926571, 323848, 12022, 1283336, 10754, 101570, 224308, 216596,
                      1004788, 266265, 266264, 99287, 1294143, 1149133, 3055, 1280, 1977402, 294, 83333,
                      536, 1407502]
    # taxa bachelor thesis
    Kleiner_taxIDs = [262724, 882, 176299, 228410, 44577, 926571, 323848, 12022, 1283336, 10754, 101570, 224308, 1041145,
                      1004788, 266265, 266264, 99287, 1302247, 1149133, 3055, 93061, 1977402, 1114970, 511145,
                      536, 1407502, 1294143,
                      1619948]

    Tanca_taxIDs = [747, 5535, 655183, 1579, 1255, 4932, 1465, 1351, 562]
    if options.taxonset == 'kleiner':
        taxonIDs = Kleiner_taxIDs
        levels = ['subspecies', 'species', 'genus', 'family', 'order', 'superkingdom']
    elif options.taxonset == 'tanca':
        taxonIDs = Tanca_taxIDs
        levels = ['species', 'genus', 'family', 'order', 'superkingdom']

    db_path = Path(options.path)
    db_type = options.database
    # for reduced_df True
    is_decoy_column_set = True
    x_tandem_result_tsv = options.input
    decoy_tag = options.decoy
    taxon_graph = HelperMethod.load_taxa_graph(Path(options.tax_graph))

    # acc_2_taxon_dict for identification file identified accs
    if db_type == 'ncbi':
        all_taxa_of_level_set = set(flatten_list(taxon_graph.get_all_taxids(taxonIDs, options.level)))
        acc_2_taxon_dict = get_ncbi_acc2taxon_dict(x_tandem_result_tsv, db_path, db_type, all_taxa_of_level_set, decoy_tag)
    else:
        acc_2_taxon_dict = get_acc2taxon_dict(x_tandem_result_tsv, db_path, db_type)

    psm = PSM_FDR(x_tandem_result_tsv)
    reduced_df = psm.create_PSM_dataframe(decoy_tag, db_type, options.level, taxon_graph, acc2tax_dict=acc_2_taxon_dict)
    fdr_pos, number_psms, decoys = psm.determine_FDR_position(reduced_df, options.fdr, is_decoy_column_set)
    print(fdr_pos, number_psms, decoys)
    # only identifications above fdr: x_tandem_result_tsv[0:psm.fdr_pos]
   # get_true_positive_and_true_negative(options.level, x_tandem_result_tsv[0:psm.fdr_pos])

if __name__ == '__main__':
    main()
