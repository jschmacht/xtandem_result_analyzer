from pyteomics import pepxml
from pathlib import Path
from collections import defaultdict
from TaxonGraph import TaxonGraph
import pandas as pd
import pickle

class ReferenceWriter():

    def __init__(self, path_to_file, path_to_tax, path_to_output, taxon_graph):
        """
        :param path_to_file: path to multiaccession file generated from NCBI-nr database
        one line one header with all accessions seperated by \t
        """
        self.path_to_xml = path_to_file
        self.output = Path(path_to_output)
        self.path_to_tax = Path(path_to_tax)
        self.taxon_graph = taxon_graph
        self.tax_level = ['species', 'genus', 'family', 'order']

    def read_pepXML(self):
        reader = pepxml.PepXML(source=str(self.path_to_xml), use_index=True, retrieve_refs=False, iterative=True)
        spectra_to_accs_dict = defaultdict(list)
        for spectrum_query in reader.iterfind("spectrum_query"):
            for search_hit in spectrum_query['search_hit']:
                percolator_q_value = search_hit['search_score']['Percolator q-Value']
                for protein in search_hit['proteins']:
                    protein_acc = protein['protein'].split()[0]
                    if '_WP_' in protein_acc:
                        protein_acc='WP_' + protein_acc.split('WP_')[1]
                    spectra_to_accs_dict[spectrum_query['spectrum']].append((protein_acc, percolator_q_value))
        return spectra_to_accs_dict

    def read_acc_to_taxid_file(self):
        acc_to_tax_dict = {}
        with open(self.path_to_tax, 'r') as tax:
            for line in tax.readlines():
                fields = line.split()
                taxid = fields[-1]
                acc = fields[1]
                if '_WP_' in acc:
                    acc = 'WP_' + acc.split('WP_')[1]
                acc_to_tax_dict[acc]=taxid
        return acc_to_tax_dict


    def load_taxa_graph(self):
        """
        # Try load pre-builded taxonomy graph or built taxonomy graph now
        :param options: user input options
        :return: TaxonGraph object
        """

        if not (self.path_to_taxdump.parents[0] / 'taxon_graph_results').is_file():
            taxon_graph = TaxonGraph()
            print("Start building taxon graph.")
            taxon_graph.create_graph(str(self.path_to_taxdump))
            print("Taxon graph successfully build.")
            # save TaxonGraph to harddrive:
            try:
                with open(str(self.path_to_taxdump.parents[0] / 'taxon_graph_results'), 'wb') as handle:
                    pickle.dump(taxon_graph, handle, protocol=pickle.HIGHEST_PROTOCOL)
                    print('Safe taxon graph to location: %s' % str(
                        self.path_to_taxdump.parents[0] / 'taxon_graph_results'))
            except FileNotFoundError:
                print('Error open tax_graph.')
                exit(1)
        # load Taxon Graph
        else:
            try:
                print('Load taxon graph from harddrive.')
                with open(str(self.path_to_taxdump.parents[0] / 'taxon_graph_results'), 'rb') as handle:
                    taxon_graph = pickle.load(handle)
            except UnicodeDecodeError or EOFError:
                print(
                    "Failed opening path to taxon graph / taxon_graph is corrupted. Delete %s file."
                    % str(self.path_to_taxdump.parents[0] / 'taxon_graph'))
                exit(1)
        return taxon_graph

    def determine_level_specific_taxIDs(self, taxid):
        level_specific_taxids = [taxid]
        for level in self.tax_level:
            if taxid == '0':
                level_specific_taxids.append(0)
            else:
                level_specific_taxids.append(self.taxon_graph.find_level_up(int(taxid), level))
        return level_specific_taxids


    def write_Kleiner_spectrum_reference_file(self, spectra_to_accs_dict):
        acc_to_tax_dict = self.read_acc_to_taxid_file()

        print('writing ...')
        with open(self.output, 'w') as output:
            output.write('SpectraID' + '\t' + 'Ref_ProteinAcc' + '\t' + 'Ref_Hyperscore' + '\t' + 'Ref_taxID_DB' + '\t'
                         + ('\t').join('Ref_taxid_' + level for level in self.tax_level) + '\n')
            for spectra, protein_list in spectra_to_accs_dict.items():
                level_specific_taxids = []
                for protein in protein_list:
                    taxid = acc_to_tax_dict[protein[0]]
                    if not level_specific_taxids:
                        level_specific_taxids = self.determine_level_specific_taxIDs(taxid)
                        level_specific_taxids = [{int(taxid)} for taxid in level_specific_taxids]
                    else:
                        l = self.determine_level_specific_taxIDs(taxid)
                        for i, taxid in enumerate(l):
                            level_specific_taxids[i].add(int(taxid))
                list_of_tax_str = [(', ' ).join([str(taxid) for taxid in taxid_set]) for taxid_set in level_specific_taxids]
                output.write('Run1_' + spectra + '\t' + protein[0] + '\t' + str(protein[1]) + '\t' +
                                  ('\t' ).join(list_of_tax_str) + '\n')

    @staticmethod
    def write_result_spectrum_reference_file(path_to_result, levels, taxon_graph, path_to_reference, fdr_pos):
        def get_taxa_of_level(taxon_graph, taxa_set, level):
            taxa_set_of_specified_level = set()
            for taxon in taxa_set:
                if taxon == 0 or taxon == "'DECOY'":
                    taxa_set_of_specified_level.add(0)
                else:
                    taxa_set_of_specified_level.add(taxon_graph.find_level_up(taxon, level))
            return taxa_set_of_specified_level

        def flatten_set(s):
            return {item for subset in s for item in subset}

        result_df = pd.read_csv(str(path_to_result), sep='\t')[0:fdr_pos]
        # remove all Decoys
        result_df = result_df[result_df.decoy != "'DECOY'"]
        reference_df = result_df[['Title', 'Peptide', 'Hyperscore', 'Protein', 'decoy', 'taxID']]
        #rename columns
        reference_df.columns = ['SpectraID', 'Ref_Peptide', 'Ref_Hyperscore', 'Ref_ProteinAcc', 'Ref_decoy', 'Ref_taxID_DB']
        reference_df['Ref_taxID_DB'] = reference_df['Ref_taxID_DB'].apply(lambda taxid_set: {int(taxon) for taxon in taxid_set[1:-1].split(', ') if taxon != "'DECOY'"})
        reference_df['Ref_decoy'] = reference_df['Ref_decoy'].apply(lambda decoy_set: {bool(decoy) for decoy in decoy_set[1:-1].split(', ')})
        for level in levels:
            reference_df[f'Ref_taxID_{level}'] = reference_df['Ref_taxID_DB'].apply(lambda taxid_set: get_taxa_of_level(taxon_graph, taxid_set, level))

        # merge all results for one spectrum
        reduced_df = reference_df.groupby(["SpectraID"], as_index=False).agg(
            {'Ref_Peptide': lambda acc: set(acc),
             'Ref_Hyperscore': lambda x: set(list(x)),
             'Ref_ProteinAcc': lambda acc_strings: flatten_set([set(acc_string[1:-1].split(', ')) for acc_string in acc_strings]),
             'Ref_decoy': lambda decoy_sets: flatten_set(decoy_sets),
             'Ref_taxID_DB': lambda taxid_sets: flatten_set(taxid_sets),
             'Ref_taxID_species': lambda taxid_sets: flatten_set(taxid_sets),
             'Ref_taxID_genus': lambda taxid_sets: flatten_set(taxid_sets),
             'Ref_taxID_family': lambda taxid_sets: flatten_set(taxid_sets),
             'Ref_taxID_order': lambda taxid_sets: flatten_set(taxid_sets),
             })
        print(f"writing refernece data frame to {path_to_reference}... ")
        reduced_df.to_csv(path_to_reference, sep='\t')

    def write_reduced_spectrum_reference_file(self, spectra_to_accs_dict):
        acc_to_tax_dict = self.read_acc_to_taxid_file()
        print('writing ...')
        with open(self.output, 'w') as output:
            output.write('#spectra_ID' + '\t' + 'protein_acc' + '\t' + 'score' + '\t' + 'taxid_DB' + '\t'
                         + ('\t').join('taxid_' + level for level in self.tax_level) + '\n')
            for spectra, protein_list in spectra_to_accs_dict.items():
                for protein in protein_list:
                    level_specific_taxids = self.determine_level_specific_taxIDs(acc_to_tax_dict[protein[0]])
                    output.write('Run1_' + spectra + '\t' + protein[0] + '\t' + str(protein[1]) + '\t' +
                                 ('\t' ).join(str(taxid) for taxid in level_specific_taxids) + '\n')


def main():
    reader = ReferenceWriter('/home/jules/Documents/Tax2Proteome/benchmarking/spectra/Run1_U1_2000ng.pep.xml',
                          '/home/jules/Documents/Tax2Proteome/benchmarking/Kleiner_ref_db/acc2tax_custom',
                          '/home/jules/Documents/Tax2Proteome/benchmarking/spectra/Run1_U1_2000ng_spectra_ref.tsv',
                          '/home/jules/Documents/databases/databases_tax2proteome/taxdump.tar.gz')
    spectra_to_accs_dict = reader.read_pepXML()
    reader.write_Kleiner_spectrum_reference_file(spectra_to_accs_dict)


if __name__ == '__main__':
    main()