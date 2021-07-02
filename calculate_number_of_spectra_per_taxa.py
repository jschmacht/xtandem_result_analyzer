from collections import defaultdict

import pandas as pd
from pathlib import Path
import argparse
from create_reference_from_tsv_and_pepxml import ReferenceWriter
from create_PSM_df import PSM_FDR
from handling_acc_files import HelperMethod


def flatten_set(l):
    return {item for sublist in l for item in sublist}


def write_results(result_dict, number_all_identified_taxa, taxonIDs, taxids_of_level, result_output):
    with open(str(result_output), 'w') as output:
        output.write(f"tax \t taxa_level \t number_of_spectra \t % of all identified spectra\n")
        for taxon, taxon_level in zip(taxonIDs, taxids_of_level):
            output.write(f"{taxon}\t{taxon_level}\t{len(result_dict[taxon_level])}\t{len(result_dict[taxon_level])/number_all_identified_taxa*100}\n")
        output.write(f"{'CRAP'}\t{'CRAP'}\t{len(result_dict['CRAP'])}\t{len(result_dict['CRAP'])/number_all_identified_taxa*100}\n")
        output.write(f"{'DECOY'}\t{'DECOY'}\t{len(result_dict['DECOY'])}\t{len(result_dict['DECOY'])/number_all_identified_taxa*100}\n")


def main():
    parser = argparse.ArgumentParser(description='Read reduced .tsv')
    parser.add_argument('-r', '--result', dest='result', default=None, help='Path to reduced result tsv, Tax2Proteome X-Tandem results')
    parser.add_argument('-l', '--level', dest='level', choices=['species', 'genus', 'family', 'order'],
                        help='Level of database')
    parser.add_argument('-g', '--tax_graph', dest='tax_graph', help='Path to taxdump.tar.gz')
    parser.add_argument('-f', '--fdr', dest='fdr', type=float, default=0.05, help='FDR-rate, default  = 0.05')
    parser.add_argument('-s', '--taxonset', dest='taxonset', choices=['tanca', 'kleiner'], default=None,
                        help='Taxon dataset used for analysis.')
    parser.add_argument('-t', '--taxon', dest='taxon', type=int, nargs='+', action='append',
                        help='NCBI taxon ID/s for database extraction. Multiple taxonIDs seperated by space.')
    parser.add_argument('-o', '--output', dest='output', default=None, help='Path to result output')
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
    if options.taxon:
        taxonIDs.extend([taxID for taxonlist in options.taxon for taxID in taxonlist])
    taxon_graph = HelperMethod.load_taxa_graph(Path(options.tax_graph))
    path_to_result_file = Path(options.result)
    if options.output:
        path_to_output=Path(options.output)
    else:
        path_to_output = path_to_result_file.parent.joinpath(path_to_result_file.stem + '_identification_results.tsv')

    result_df = ReferenceWriter.read_csv_with_generic_function(path_to_result_file,
                                                               ['Protein', 'Hyperscore', 'decoy', 'taxID', f'taxID_{options.level}'])
    fdr_pos, number_psm, number_decoy, double_spectra  = PSM_FDR.determine_FDR_position(result_df, options.fdr, True)
    fdr_applied_result_df = result_df[0:fdr_pos]
    taxids_of_level = [taxon_graph.find_level_up(taxid, options.level) for taxid in taxonIDs]
    result_dict = defaultdict(set)
    for spectra, taxid_set in zip(fdr_applied_result_df['Title'].tolist(), fdr_applied_result_df[f'taxID_{options.level}'].tolist()):
        for taxid in taxid_set:
            result_dict[taxid].add(spectra)
    number_all_identified_spectra = len(flatten_set(list(result_dict.values())))
    write_results(result_dict, number_all_identified_spectra, taxonIDs, taxids_of_level, path_to_output)


if __name__ == '__main__':
    main()