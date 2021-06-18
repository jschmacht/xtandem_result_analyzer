from pathlib import Path
import argparse
from TaxonGraph import TaxonGraph
import pandas as pd
import pickle
from PSM_FDR import PSM_FDR
from AccessionSearcher import AccessionSearcherNCBI
from ReadAccTaxon import ReadAccTaxon
from SearchAccessions import SearchAccessions
from collections import defaultdict


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


def find_taxonID(accessions, acc_taxon_dict, crap_set,db):
    """
    :param accessions: for every xtandem spectrum list of accession -> [[]]
    :param acc_taxon_dict:  key=accession value=taxonID (from file generated from whole uniprot db)
    :param crap_set: all crap accessions
    """

    taxons_identified = []
    for acc_list in accessions:
        taxons = set()
        for acc in acc_list:
            if 'REVERSED' in acc:
                taxons.add(-1)
                continue
            try:
                if db=='ncbi':
                    taxons.add(int(acc_taxon_dict[acc]))
                elif db=='custom':
                    taxons.add(acc_taxon_dict[acc])
                else:
                    taxons.add(int(acc_taxon_dict[acc][0]))
            except KeyError:
                # acc in crap taxonID = 0
                if acc in crap_set:
                    taxons.add(0)
                else:
                    pass
                    # raise Exception ('Accession not known.')
                    #print('Accession %s not known.' % acc)
        taxons_identified.append(taxons)
    return taxons_identified

# result_TP contains for every protein_group one value (TP/FP/various)
def determine_species(taxons_identified, taxon_graph, taxonIDs, level):
    """
        :param taxons_identified: list of sets of taxon IDs, one set per identified spectrum
        :param taxon_graph object
        :param taxonIDs: list of taxon IDs from which database was created (z.B. Kleiner taxIDs)
        :param level: level for which the species are to be determined
        :return taxon_dict: key=taxon value=number appearance
        """
    print('Determine species of level ' + level)
    # taxonIDs specified to level
    if level != 'subspecies':
        taxonIDs_level = [taxon_graph.find_level_up(taxID, level) for taxID in taxonIDs]
    else:
        taxonIDs_level = taxonIDs

    # set of all descendant taxonIDs for every taxon ID of database creation
    taxonIDs_with_descendants = []
    for taxID in taxonIDs_level:
        taxonIDs_with_descendants.append(taxon_graph.find_taxIDs(taxID, childs_until_species=False))
    taxon_dict = {}
    for i, taxon_set in enumerate(taxonIDs_with_descendants):
        num = 0
        for identified_taxon_set in taxons_identified:
            # if one taxon of taxonIDs_with_descandants is in taxons_identified. if intersection of both sets is not empty
            if identified_taxon_set.intersection(taxon_set):
                num += 1
        taxon_dict[taxonIDs[i]] = num
    return taxon_dict


def write_result_file(file, number_psms, decoys, results_per_level, level_index):
    """
    :param file= path input tsv file
    :param number_psms = number identified psms = psm.number_psms
    :param decoys: psm.decoys
    :param results_per_level = list of taxon_dict: key = taxon value = number  appearance
    :param level_index = number, position level in level list
    """
    #percent or number?
    #sum_identified = sum(taxon_dict.values())
    #print(sum_identified)
    with open(file, 'w') as out:
        out.write('taxon' + '\t ' + 'number' + '\t ' + 'level_db' + '\t ' + 'level_analysis\n')
        for i, level_analysis in enumerate(levels[0:level_index + 1]):
            for key, value in results_per_level[i].items():
                # out.write(str(key) + ' \t  ' + str(value/sum_identified * 100) + '\t ' + level + '\t ' + level_analysis + '\n')
                out.write(str(key) + ' \t  ' + str(value) + '\t ' + levels[level_index] + '\t ' + level_analysis + '\n')
        out.write('\n\n')
        out.write('level_db\tPSM\tdecoy\n' + levels[level_index] + '\t' + str(number_psms) + '\t' + str(decoys) + '\n')


def write_sequence_file(file, spectra_sequence_dict, spectra_acc_dict):
    """
    :param file: path to output file
    :param spectra_sequence_dict: spectra ID to Peptide sequence dict
    :param spectra_acc_dict: spectra ID to [accs] dict
    :output file: spectrum ID, Peptide sequence, all accessions /t seperated,
    """
    with open(file, 'w') as out:
        out.write('#spectrum\t sequence\taccessions\n')
        for key, value in spectra_sequence_dict.items():
            if 'REVERSED' not in spectra_acc_dict[key]:
                out.write(str(key) + ' \t  ' + str(value) + '\t' + ' '.join(spectra_acc_dict[key]) + '\n')


def load_taxa_graph(tax_graph):
    """
    # Try load pre-builded taxonomy graph or built taxonomy graph now
    :param options: user input options
    :return: TaxonGraph object
    """

    if not (Path(tax_graph).parents[0] / 'taxon_graph_results').is_file():
        taxon_graph = TaxonGraph()
        print("Start building taxon graph.")
        taxon_graph.create_graph(str(Path(tax_graph)))
        print("Taxon graph successfully build.")
        # save TaxonGraph to harddrive:
        try:
            with open(str(Path(tax_graph).parents[0] / 'taxon_graph_results'), 'wb') as handle:
                pickle.dump(taxon_graph, handle, protocol=pickle.HIGHEST_PROTOCOL)
                print('Safe taxon graph to location: %s' % str(
                    Path(tax_graph).parents[0] / 'taxon_graph_results'))
        except FileNotFoundError:
            print('Error open tax_graph.')
            exit(1)
    # load Taxon Graph
    else:
        try:
            print('Load taxon graph from harddrive.')
            with open(str(Path(tax_graph).parents[0] / 'taxon_graph_results'), 'rb') as handle:
                taxon_graph = pickle.load(handle)
        except UnicodeDecodeError or EOFError:
            print(
                "Failed opening path to taxon graph / taxon_graph is corrupted. Delete %s file."
                % str(Path(tax_graph).parents[0] / 'taxon_graph'))
            exit(1)
    return taxon_graph


def create_all_accs_dict_and_write_seq(input, decoy, fdr, database):
    psm = PSM_FDR(input)
    psm.create_PSM_dataframe(decoy_tag=decoy)
    psm.determine_FDR_position(fdr)
    print('Create spectrum_acc dict.')
    # title_spectrum: [acc, acc] every spectrum: list annotated taxons
    psm.create_spectrum_acc_dict(database)
    write_sequence_file(input + '_sequence.txt', psm.spectra_sequence_dict, psm.spectra_acc_dict)
    # all accession dict
    all_acc_dict = psm.spectra_acc_dict
    print('length all_acc_dict: %d' %len(all_acc_dict))
    return all_acc_dict, psm

def get_group_accessions_ncbi(path_to_taxdump, accs, threads, spectra2accs_dict):
    """
    :param path_to_taxdump: Path to taxdump.tar.gz
    :param accs:
    :param threads:
    :param spectra2accs_dict: psm.spectra_acc_dict spectrum to [accs] dict
    :return:
    """
    # multispecies_acc: handmade ncbi multispecies acc file: multiacc /t all accs seperated by /t
    multiacc_search_ncbi = SearchAccessions(Path(path_to_taxdump).parents[0] / 'multispecies_acc', 'ncbi')
    multiacc_search_ncbi.read_database(accs, threads)
    multi_acc_2_taxa_dict = multiacc_search_ncbi.acc_dict
    for key, value in spectra2accs_dict.items():
        group_accessions = []
        for acc in accs:
            try:
                group_accessions = group_accessions + multi_acc_2_taxa_dict[acc]
            except KeyError:
                group_accessions.append(acc)
        spectra2accs_dict[key] = group_accessions
    return spectra2accs_dict


def find_multi_acc_in_list_of_accs(acc_list):
    multiacc_marks = ['WP_', 'XP_']
    multiacc = None
    for single_accs in acc_list:
        if single_accs.startswith(multiacc_marks[0]) or single_accs.startswith(multiacc_marks[1]):
            multiacc = single_accs
    return multiacc


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


def get_tsv_multspecies_acc_to_taxa_dict(multiacc2acc_dict, multi_acc_2_taxon_dict):
    tsv_multi_acc_to_taxon_dict = defaultdict(set)
    for acc, acc_list in multiacc2acc_dict.items():
        for multi_acc in acc_list:
            if multi_acc in multi_acc_2_taxon_dict.keys():
                tsv_multi_acc_to_taxon_dict[multi_acc].add(multi_acc_2_taxon_dict[multi_acc])
    return tsv_multi_acc_to_taxon_dict

def remove_accs_with_unsupported_taxa(multi_acc_2_taxon_dict, taxa):
    accs_with_not_matching_taxa = []
    for acc, taxon in multi_acc_2_taxon_dict.items():
        if int(taxon) not in taxa:
            accs_with_not_matching_taxa.append(acc)
    for acc in accs_with_not_matching_taxa:
        del multi_acc_2_taxon_dict[acc]
    return multi_acc_2_taxon_dict

def get_ncbi_acc2taxon_dict(x_tandem_result_tsv, db_path, db_type, taxa):
    """
    :param x_tandem_result_tsv:
    :param db_path:
    :param db_type:
    :param taxa: list of all relevant taxa, all taxa in taxon graph below specified level
    :return:
    """
    # Protein: generic|AC:2905764_REVERSED|VVN60222.1-REVERSED, generic|AC:4682148|PMZ83684.1 PMZ91489.1 WP_054593944.1 PMZ34596.1 ALI00421.1 PMZ37242.1...
    protein_accs_from_result_tsv = pd.read_csv(str(x_tandem_result_tsv), delimiter='\t')['Protein'].tolist()
    ncbi_accs_from_file = [acc.split('|')[2] for acc in protein_accs_from_result_tsv if 'REVERSED' not in acc]
    print(len(ncbi_accs_from_file))
    # for ncbi find all multispecies accessions with file ncbi acc
    path_to_multiaccs = '/home/jules/Documents/databases/databases_tax2proteome/multispecies_acc'
    multiaccs = set()
    complete_accs = set()
    for i, acc in enumerate(ncbi_accs_from_file):
        # if line ending with ...: not complete
        if len(acc.split()) > 1:
            multiaccs.add(acc.split()[0].strip())
        else:
            complete_accs.add(acc.strip())
    print('complete_accs: ', len(complete_accs))
    print('multiaccs: ', len(multiaccs))

    # key = first acc in result tsv
    multiacc2acc_dict = get_multispecies_accs(path_to_multiaccs, multiaccs)
    print('multiaccs ready', len(multiacc2acc_dict))
    multi_acc2tax_reader = ReadAccTaxon(db_path, db_type)
    # key = Multiacc or accs assigned to multiaccs
    multi_acc_2_taxon_dict = multi_acc2tax_reader.read_acc2tax(multiaccs)
    print('multi_acc_2_taxon_dict 1: ',len(multi_acc_2_taxon_dict), list(multi_acc_2_taxon_dict.values())[0:5])
    multi_acc_2_taxon_dict = remove_accs_with_unsupported_taxa(multi_acc_2_taxon_dict, taxa)
    print('multi_acc_2_taxon_dict 2: ', len(multi_acc_2_taxon_dict))
    tsv_multi_acc_to_taxa_dict = get_tsv_multspecies_acc_to_taxa_dict(multiacc2acc_dict, multi_acc_2_taxon_dict)
    # complete_accs.extend(flatten_list(multiacc2acc_dict.values()))
    print('tsv_multi_acc_to_taxon_dict: ', len(tsv_multi_acc_to_taxa_dict))
    acc2tax_reader=ReadAccTaxon(db_path, db_type)
    acc_2_taxon_dict = acc2tax_reader.read_acc2tax(complete_accs)
    # taxon to set(taxon), in tsv_multi_acc_to_taxa_dict multiple taxa in set possible
    acc_2_taxon_dict = {acc:{taxon} for (acc,taxon) in acc_2_taxon_dict.items()}
    print('acc2tax_dict ready1', len(acc_2_taxon_dict))
    acc_2_taxon_dict.update(tsv_multi_acc_to_taxa_dict)
    print('acc2tax_dict ready2', len(acc_2_taxon_dict))
    return acc_2_taxon_dict


def flatten_list(l):
    return [item for sublist in l for item in sublist]


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
                      536, 1407502]

    Tanca_taxIDs = [747, 5535, 655183, 1579, 1255, 4932, 1465, 1351, 562]
    if options.taxonset == 'kleiner':
        taxonIDs = Kleiner_taxIDs
        levels = ['subspecies', 'species', 'genus', 'family', 'order', 'superkingdom']
    elif options.taxonset == 'tanca':
        taxonIDs = Tanca_taxIDs
        levels = ['species', 'genus', 'family', 'order', 'superkingdom']

    db_path = Path(options.path)
    db_type = options.database
    x_tandem_result_tsv = options.input
    decoy_tag = options.decoy
    taxon_graph = load_taxa_graph(options.tax_graph)
    all_taxa_of_level_set = set(flatten_list(taxon_graph.get_all_taxids(taxonIDs, options.level)))
    # acc_2_taxon_dict for identification file identified accs
    if db_type == 'ncbi':
        acc_2_taxon_dict = get_ncbi_acc2taxon_dict(x_tandem_result_tsv, db_path, db_type, all_taxa_of_level_set)
    else:
        acc_2_taxon_dict = get_acc2taxon_dict(x_tandem_result_tsv, db_path, db_type)


    psm = PSM_FDR(x_tandem_result_tsv)
    reduced_df = psm.create_PSM_dataframe(decoy_tag, db_type, options.level, taxon_graph, acc2tax_dict=acc_2_taxon_dict)
    fdr_pos, number_psms, decoys = psm.determine_FDR_position(reduced_df, options.fdr)
    # only identifications above fdr: x_tandem_result_tsv[0:psm.fdr_pos]
   # get_true_positive_and_true_negative(options.level, x_tandem_result_tsv[0:psm.fdr_pos])

if __name__ == '__main__':
    main()
