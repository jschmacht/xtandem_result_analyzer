from pathlib import Path
import argparse
from TaxonGraph import TaxonGraph
import pandas as pd
import pickle
from PSM_FDR import PSM_FDR
from AccessionSearcher import AccessionSearcherNCBI
from ReadAccTaxon import ReadAccTaxon
from SearchAccessions import SearchAccessions


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
                else:
                    taxons.add(int(acc_taxon_dict[acc][0]))
            except KeyError:
                # acc in crap taxonID = 0
                if acc in crap_set:
                    taxons.add(0)
                else:
                    # raise Exception ('Accession not known.')
                    print('Accession %s not known.' % acc)
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


def writer(file, psm, results_per_level, level_index):
    """
            :param file= path input tsv file
            :param psm = number identified psms =psm.psm
            :param results_per_level = list of taxon_dict: key = taxon value = number  appearance
            :param level_index = number, position level in level list

            """
    #percent or number?
    #sum_identified = sum(taxon_dict.values())
    #print(sum_identified)
    with open(file+'_result.txt', 'w') as out:
        out.write('taxon' + '\t ' + 'number' + '\t ' + 'level_db' + '\t ' + 'level_analysis\n')
        print(levels[0:level_index + 1])
        for i, level_analysis in enumerate(levels[0:level_index + 1]):
            for key, value in results_per_level[i].items():
                # out.write(str(key) + ' \t  ' + str(value/sum_identified * 100) + '\t ' + level + '\t ' + level_analysis + '\n')
                out.write(str(key) + ' \t  ' + str(value) + '\t ' + levels[level_index] + '\t ' + level_analysis + '\n')
        out.write('\n\n')
        out.write('level_db\tPSM\tdecoy\n' + levels[level_index] + '\t' + str(psm.psm) + '\t' + str(psm.decoys) + '\n')
        print("Outfile written to " + file + '_result.txt')


def writer_sequence(file, psm):
    """
                :param file = path to folder
                :param psm = psm object
                                """
    with open(file + '_sequence.txt', 'w') as out:
        out.write('spectrum\t sequence\taccessions\n')
        for key, value in psm.spectra_sequence_dict.items():
            if 'REVERSED' not in psm.spectra_acc_dict[key]:
                out.write(str(key) + ' \t  ' + str(value) + '\t' + ' '.join(psm.spectra_acc_dict[key]) + '\n')
        print("Outfile written to " + file + '_sequence.txt')
    

def main():
    parser = argparse.ArgumentParser(description='Read xtandem output .tsv')
    parser.add_argument('-i', '--input', dest='input', default=None, help='xtandem tsv')
    parser.add_argument('-p', '--path', dest='path', default=None, help='Path to folder with acc2tax files.')
    parser.add_argument('-c', '--crap', dest='crap', default=None, help='Crap-File')
    parser.add_argument('-g', '--tax_graph', dest='tax_graph', help='Path to taxdump.tar.gz')
    parser.add_argument('-l', '--level', dest='level', choices=['subspecies', 'species', 'section', 'genus', 'family', 'order', 'superkingdom'],
                        help='Level of database')
    parser.add_argument('-d', '--database', dest='database', choices=['ncbi', 'uniprot'], default='uniprot',
                        help='Database format.')
    parser.add_argument('-r', '--non_redundant', dest='non_redundant', default=None, help='Path to non redundant uniprot database.')
    parser.add_argument('-s', '--taxonset', dest='taxonset', choices=['tanca', 'kleiner'], default=None,
                        help='Taxon dataset used for analysis.')
    parser.add_argument('-t', '--taxon', dest='taxon', default=None, type=int, nargs='+', action='append',
                        help='NCBI taxon ID/s for database extraction. Multiple taxonIDs seperated by space.')
    parser.add_argument('-f', '--fdr', dest='fdr', type=float, default=0.01, help='FDR-rate, default  = 0.01')
    parser.add_argument('-y', '--decoy', dest='decoy', default='REVERSED', help='Decoy_tag.')

    parser.add_argument('-x', '--threads', dest='threads', type=int, action="store", help='Number of threads.')

    options = parser.parse_args()
    Kleiner_taxIDs = [
        176299, 1004788, 224308, 266265, 3055, 536, 266264, 882, 511145, 228410, 44577, 926571, 323848, 1302247, 101570,
        1283336, 12022, 1977402, 10754, 1294143, 1114970, 1149133, 216596, 1041145, 99287, 46170, 93061, 1407502, 262724
    ]
    Tanca_taxIDs = [747, 5535, 655183, 1579, 1255, 4932, 1465, 1351, 562]
    global levels
    levels = ['subspecies', 'species', 'genus', 'family', 'superkingdom']
    if options.taxonset == 'kleiner':
        taxonIDs = Kleiner_taxIDs
        levels = ['subspecies', 'species','section', 'genus', 'family', 'order', 'superkingdom']
    elif options.taxonset == 'tanca':
        taxonIDs = Tanca_taxIDs
        levels = ['species', 'genus', 'family', 'order', 'superkingdom']
    elif options.taxon:
        taxonIDs = options.taxon[0]
        
    print('TaxonIDs of database: ' + ' '.join([str(x) for x in taxonIDs]))


    # Try load pre-builded taxonomy graph or built taxonomy graph now
    if not (Path(options.tax_graph).parents[0] / 'taxon_graph_results').is_file():
        taxon_graph = TaxonGraph()
        print("Start building taxon graph.")
        taxon_graph.create_graph(str(Path(options.tax_graph)))
        print("Taxon graph successfully build.")
        # save TaxonGraph to harddrive:
        try:
            with open(str(Path(options.tax_graph).parents[0] / 'taxon_graph_results'), 'wb') as handle:
                pickle.dump(taxon_graph, handle, protocol=pickle.HIGHEST_PROTOCOL)
                print('Safe taxon graph to location: %s' % str(
                    Path(options.tax_graph).parents[0] / 'taxon_graph_results'))
        except FileNotFoundError:
            print('Error open tax_graph.')
            exit(1)
    # load Taxon Graph
    else:
        try:
            print('Load taxon graph from harddrive.')
            with open(str(Path(options.tax_graph).parents[0] / 'taxon_graph_results'), 'rb') as handle:
                taxon_graph = pickle.load(handle)
        except UnicodeDecodeError or EOFError:
            print(
                "Failed opening path to taxon graph / taxon_graph is corrupted. Delete %s file."
                % str(Path(options.tax_graph).parents[0] / 'taxon_graph'))
            exit(1)

    psm = PSM_FDR(options.input)
    psm.create_PSM_dataframe()
    psm.determine_FDR_position(options.decoy, options.fdr)
    print('Create spectrum_acc dict.')
    # title_spectrum: [acc, acc] every spectrum: list annotated taxons
    psm.create_spectrum_acc_dict(options.database)
    writer_sequence(options.input, psm)

    # all accession dict
    all_acc_dict = psm.spectra_acc_dict
    print('length all_acc_dict: %d' %len(all_acc_dict))
    # for ncbi find all multispecies accessions with file ncbi acc : (uniprot only one acc in all_acc_dict)
    # and for uniprot non redundant
    if options.database == 'ncbi':
        #multiacc_search_ncbi = AccessionSearcherNCBI(Path(options.tax_graph).parents[0] / 'multispecies_acc')
        #multiacc_search_ncbi.search_ncbi_accession(all_acc_dict.values())
        #acc_multiacc_ncbi_dict = multiacc_search_ncbi.acc_dict
        
        #print(acc_multiacc_ncbi_dict['WP...'])
        multiacc_search_ncbi = SearchAccessions(Path(options.tax_graph).parents[0] / 'multispecies_acc', 'ncbi')
        multiacc_search_ncbi.read_database(psm.spectra_acc_dict.values(), options.threads)
        acc_multiacc_ncbi_dict = multiacc_search_ncbi.acc_dict
        print('length acc_multiacc_ncbi_dict: %d' % len(acc_multiacc_ncbi_dict))
        for key, value in all_acc_dict.items():
            group_accessions = []
            for acc in value:
                try:
                    group_accessions = group_accessions + acc_multiacc_ncbi_dict[acc]
                except KeyError:
                    group_accessions.append(acc)
            all_acc_dict[key] = group_accessions
    print('Spectrum: [Accessions] dict ready.')
    print('length all_acc_dict: %d' % len(all_acc_dict))
    #print(all_acc_dict['WP...'])
    if options.non_redundant and options.database == 'uniprot':
        multiacc_search_nr = SearchAccessions(options.non_redundant, 'uniprot')
        multiacc_search_nr.read_database(psm.spectra_acc_dict.values(), options.threads)

    accessions = set([item for sublist in all_acc_dict.values() for item in sublist])
    print('length accessions: %d' % len(accessions))
    acc2taxon = ReadAccTaxon(options.path)
    acc2taxon.read_acc2tax(options.database, accessions, options.threads)
    acc_taxon_dict = acc2taxon.acc_taxon_dict
    print('length acc_taxon_dict: %d' % len(acc_taxon_dict))
    #print(acc_taxon_dict['WP..'])
    print('Accession:taxon dict ready')

    # determine taxons for every set, taxon_list: list of sets of taxon IDs, one set per identified spectrum
    crap_set = read_crap(str(options.crap))
    taxons_identified = find_taxonID(all_acc_dict.values(), acc_taxon_dict, crap_set, options.database)

    results_per_level = []
    level_index = levels.index(options.level)
    # results_per_level = list of taxon_dict: key=taxon value=number appearance
    for level in levels[0:level_index + 1]:
        results_per_level.append(determine_species(taxons_identified, taxon_graph, taxonIDs, level))
    writer(options.input, psm, results_per_level, level_index)



if __name__ == '__main__':
    main()
