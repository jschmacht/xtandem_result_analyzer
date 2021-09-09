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

def get_acc_taxon_dict(spectra2accs_dict, path, db_type, threads):
    accessions = set([item for sublist in spectra2accs_dict.values() for item in sublist])
    print('length accessions: %d' % len(accessions))
    acc2taxon = ReadAccTaxon(path, db_type)
    if db_type=='custom':
        acc_taxon_dict = acc2taxon.read_custom_acc2tax()
    else:
        acc_taxon_dict = acc2taxon.read_acc2tax(accessions, threads)
    print('length acc_taxon_dict: %d' % len(acc_taxon_dict))
    return acc_taxon_dict


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
            if fields[0] in accs:
                multiacc2accs_dict[fields[0]] = fields[1:]
        print('multi accs found')
    return multiacc2accs_dict


def flatten_list(l):
    return [item for sublist in l for item in sublist]


def find_psms_per_level(taxon_graph, taxonIDs, spectra2accs_dict, acc_2_taxon_dict, crap_set, db_type, level):
    taxons_identified = find_taxonID(spectra2accs_dict.values(), acc_2_taxon_dict, crap_set, db_type)
    results_per_level = []
    level_index = levels.index(level)
    # results_per_level = list of taxon_dict: key=taxon value=number appearance
    for level in levels[0:level_index + 1]:
        results_per_level.append(determine_species(taxons_identified, taxon_graph, taxonIDs, level))
    return results_per_level, taxons_identified

def load_ref_file(ref_file, level):
    level_to_column_nb_dict={'species': 4, 'genus': 5, 'family': 6, 'order': 7}
    spectraID_to_taxid_dict = defaultdict(list)
    with open(ref_file, 'r') as ref:
        ref.readline()
        for line in ref:
            fields = line.split()
            level_specific_taxid = fields[level_to_column_nb_dict[level]]
            spectraID = fields[0]
            spectraID_to_taxid_dict[spectraID].append(level_specific_taxid)
    return spectraID_to_taxid_dict


def get_all_spectra_IDs(ident_file):
    all_spec_IDs = set()
    with open(ident_file, 'r') as ident_file:
        for line in ident_file:
            if line.startswith('TITLE'):
                all_spec_IDs.add(line.split()[0].split('TITLE=')[1])
    return all_spec_IDs


def get_true_positive_and_true_negative(level, x_tandem_result_tsv):
    ref_file = '/home/jules/Documents/Tax2Proteome/benchmarking/spectra/Run1_U1_2000ng_spectra_ref.tsv'
    ident_file = '/home/jules/Documents/Tax2Proteome/benchmarking/spectra/Run1_U1_2000ng.mgf'
    all_spectra_IDs = get_all_spectra_IDs(ident_file)
    df = pd.read_csv(str(x_tandem_result_tsv +'.csv'))
    df['taxID_level'] = df['taxID_level'].apply(lambda taxid: taxid[1:-1].split(', '))
    df['Protein'] = df['Protein'].apply(lambda taxid: taxid[1:-1].split(', '))
    spectraID_to_taxid_dict = load_ref_file(ref_file, level)
    df['TP'] = df.apply(lambda row: 'FP' if row['Title'] not in spectraID_to_taxid_dict.keys() else (
        'TP' if (set(spectraID_to_taxid_dict[row['Title']]) & set(df['taxID_level'])) else ''))
    number_TN = 0
    number_FN = 0
    for spectra in all_spectra_IDs:
        if spectra not in df['Title'] and spectra not in spectraID_to_taxid_dict.keys():
            number_TN += 1
        if spectra not in df['Title'] and spectra in spectraID_to_taxid_dict.keys():
            number_FN += 1
    number_TP = df[df['TP'] == 'TP'].count()
    number_FP = df[df['TP'] == 'FP'].count()
    print(f"TP: {number_TP}, FP: {number_FP}, TN: {number_TN}, FN: {number_FN}")




def main():
    """
    :return:
    """
    parser = argparse.ArgumentParser(description='Read xtandem output .tsv')
    parser.add_argument('-i', '--input', dest='input', default=None, help='xtandem tsv with columns:'
                                '[spectra_file, Spectrum ID, Peptide, Protein Accession + Description, Hyperscore')
    parser.add_argument('-p', '--path', dest='path', default=None, help='Path to folder with acc2tax files. '
                                                                        '(Uniprot: acc2tax_uniprot, NCBI: pdb/prot.accession2taxid, '
                                                                        'Custom: acc2tax_custom')
    parser.add_argument('-c', '--crap', dest='crap', default=None, help='Crap-File')
    parser.add_argument('-g', '--tax_graph', dest='tax_graph', help='Path to taxdump.tar.gz')
    parser.add_argument('-l', '--level', dest='level', choices=['subspecies', 'species', 'genus', 'family', 'order', 'superkingdom'],
                        help='Level of database')
    parser.add_argument('-d', '--database', dest='database', choices=['ncbi', 'uniprot', 'custom'], default='uniprot',
                        help='Database format.')
    parser.add_argument('-r', '--non_redundant', dest='non_redundant', default=None, help='Path to tax2proteome generated '
                                                                                          'non redundant uniprot database.')
    parser.add_argument('-s', '--taxonset', dest='taxonset', choices=['tanca', 'kleiner'], default=None,
                        help='Taxon dataset used for analysis.')
    parser.add_argument('-t', '--taxon', dest='taxon', default=None, type=int, nargs='+', action='append',
                        help='NCBI taxon ID/s for database extraction. Multiple taxonIDs seperated by space.')
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
    global levels
    levels = ['subspecies', 'species', 'genus', 'family', 'superkingdom']
    if options.taxonset == 'kleiner':
        taxonIDs = Kleiner_taxIDs
        levels = ['subspecies', 'species', 'genus', 'family', 'order', 'superkingdom']
    elif options.taxonset == 'tanca':
        taxonIDs = Tanca_taxIDs
        levels = ['species', 'genus', 'family', 'order', 'superkingdom']
    elif options.taxon:
        taxonIDs = options.taxon[0]
        
    print('TaxonIDs of database: ' + ' '.join([str(x) for x in taxonIDs]))
    db_path = Path(options.path)
    db_type = options.database
    x_tandem_result_tsv = options.input
    decoy_tag = options.decoy
    path_to_uniprot_nr = options.non_redundant
    # TaxonGraph() object
    taxon_graph = load_taxa_graph(options.tax_graph)
    acc2tax_reader=ReadAccTaxon(db_path, db_type)
    if db_type == 'custom':
        acc_2_taxon_dict = acc2tax_reader.get_acc2taxonID_dict(db_path/'acc2tax_custom')
    elif db_type == 'uniprot':
        # taxids_level = taxon_graph.get_all_taxids(taxonIDs, options.level)
        accs = {acc.split()[0].split('|')[1] for acc in pd.read_csv(str(x_tandem_result_tsv), delimiter='\t')['Protein'].tolist()}
        acc_2_taxon_dict = acc2tax_reader.read_acc2tax(accs)
    elif db_type == 'ncbi':
        accs = set(flatten_list([[acc for acc in acc.split('|')[2].split() if 'REVERSED' not in acc]
                                 for acc in pd.read_csv(str(x_tandem_result_tsv), delimiter='\t')['Protein'].tolist()]))
        # for ncbi find all multispecies accessions with file ncbi acc
        # accs.union(set(flatten_list(get_multispecies_accs(db_path/'multispecies_acc', accs).values())))
        acc_2_taxon_dict = acc2tax_reader.read_acc2tax(accs)
    psm = PSM_FDR(x_tandem_result_tsv)
    psm.create_PSM_dataframe(decoy_tag, db_type, options.level, taxon_graph, acc2tax_dict=acc_2_taxon_dict)
    psm.determine_FDR_position(options.fdr)
    # only identifications above fdr: x_tandem_result_tsv[0:psm.fdr_pos]
    get_true_positive_and_true_negative(options.level, x_tandem_result_tsv[0:psm.fdr_pos])
       # if path_to_uniprot_nr and db_type == 'uniprot':
        #    multiacc_search_nr = SearchAccessions(path_to_uniprot_nr, 'uniprot')
         #   multiacc_search_nr.read_database(psm.spectra_acc_dict.values(), options.threads)

        # determine taxons for every set, taxon_list: list of sets of taxon IDs, one set per identified spectrum
     #   crap_set = read_crap(str(options.crap))

     #   results_per_level, taxons_identified = find_psms_per_level(taxon_graph, taxonIDs, spectra2accs_dict, acc_2_taxon_dict, crap_set, db_type, options.level)
     #   write_result_file(x_tandem_result_tsv+'_result.txt', psm.number_psms, psm.decoys, results_per_level, levels.index(options.level))


if __name__ == '__main__':
    main()
