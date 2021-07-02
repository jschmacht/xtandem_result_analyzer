from pyteomics import pepxml
from pathlib import Path
from collections import defaultdict
import pandas as pd
import ast
from create_PSM_df import PSM_FDR


class ReferenceWriter():

    def __init__(self, reference_df, path_to_reference, taxon_graph):
        """
        :param reference_df: reduced df of x-tandem search result against Kleiner DB
        :param path_to_reference: path to output spectrum reference file
        :param taxon_graph: loaded TaxonGraph
        :param fdr: choosen fdr value
        """
        self.reference_df = reference_df
        self.path_to_reference_output = Path(path_to_reference)
        self.taxon_graph = taxon_graph
        self.tax_level = ['species', 'genus', 'family', 'order']

    def read_pepXML(self, path_to_spectrum_identification_file):
        reader = pepxml.PepXML(source=str(path_to_spectrum_identification_file), use_index=True,
                               retrieve_refs=False, iterative=True)
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

    def read_acc_to_taxid_file(self, path_to_custom_tax):
        acc_to_tax_dict = {}
        with open(path_to_custom_tax, 'r') as tax:
            for line in tax.readlines():
                fields = line.split()
                taxid = fields[-1]
                acc = fields[1]
                if '_WP_' in acc:
                    acc = 'WP_' + acc.split('WP_')[1]
                acc_to_tax_dict[acc]=taxid
        return acc_to_tax_dict

    def determine_level_specific_taxIDs(self, taxid):
        level_specific_taxids = [taxid]
        for level in self.tax_level:
            if taxid == '0':
                level_specific_taxids.append(0)
            else:
                level_specific_taxids.append(self.taxon_graph.find_level_up(int(taxid), level))
        return level_specific_taxids

    def write_Kleiner_spectrum_reference_file(self, spectra_to_accs_dict, path_to_custom_tax):
        acc_to_tax_dict = self.read_acc_to_taxid_file(path_to_custom_tax)
        print('writing ...')
        with open(self.path_to_reference_output, 'w') as output:
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

    def get_taxa_of_level(self, taxa_set, level):
        taxa_set_of_specified_level = set()
        for taxon in taxa_set:
            if  taxon == 'DECOY' or taxon == 'DECOY/CRAP' or taxon == 'CRAP':
                taxa_set_of_specified_level.add(taxon)
            elif taxon == 0:
                print(taxa_set)
            else:
                taxa_set_of_specified_level.add(self.taxon_graph.find_level_up(taxon, level))
        return taxa_set_of_specified_level

    def flatten_set(self, s):
        return {item for subset in s for item in subset}

    def write_result_spectrum_reference_file(self):
        for level in self.tax_level:
            self.reference_df[f'taxID_{level}'] = self.reference_df['taxID'].apply(lambda taxid_set:
                                                                self.get_taxa_of_level(taxid_set, level))
        print(f"writing refernece data frame to {self.path_to_reference_output}... ")
        #reduced_df = reduced_df.sort_values(by=['Ref_Hyperscore', 'SpectraID'], ascending=False).reset_index(drop=True)
        self.reference_df.to_csv(str(self.path_to_reference_output), sep='\t')
        # merge all results for one spectrum
        # reduced_df = reference_df.groupby(["SpectraID"], as_index=False).agg(
        #     {'Ref_Peptide': lambda acc: set(acc),
        #      'Ref_Hyperscore': lambda scores: set(scores),
        #      'Ref_ProteinAcc': lambda acc_strings: self.flatten_set(acc_strings),
        #      'Ref_decoy': lambda decoy_sets: self.flatten_set(decoy_sets),
        #      'Ref_taxID_DB': lambda taxid_sets: self.flatten_set(taxid_sets),
        #      'Ref_taxID_species': lambda taxid_sets: self.flatten_set(taxid_sets),
        #      'Ref_taxID_genus': lambda taxid_sets: self.flatten_set(taxid_sets),
        #      'Ref_taxID_family': lambda taxid_sets: self.flatten_set(taxid_sets),
        #      'Ref_taxID_order': lambda taxid_sets: self.flatten_set(taxid_sets),
        #      })


    def create_reference_file_for_kleiner_pep_xml_identification(self, path_to_spectrum_identification_file):
        if path_to_spectrum_identification_file.suffixes == ['.pep', '.xml']:
            custom_acc2tax_file_based_on_Kleiner_DB = '/home/jules/Documents/Tax2Proteome/benchmarking/Kleiner_ref_db/acc2tax_custom'
            spectra_to_accs_dict = self.read_pepXML(path_to_spectrum_identification_file)
            self.write_Kleiner_spectrum_reference_file(spectra_to_accs_dict, custom_acc2tax_file_based_on_Kleiner_DB)

    @staticmethod
    def read_csv_with_generic_function(path_to_df, conv_columns):
        generic_read_csv_function = lambda x: ast.literal_eval(x)
        conv_dict = {}
        for column in conv_columns:
            conv_dict[column] = generic_read_csv_function
        df = pd.read_csv(str(path_to_df), sep='\t', converters=conv_dict)
        return df


