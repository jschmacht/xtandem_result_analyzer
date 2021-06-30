from pathlib import Path
from ReadAccTaxon import ReadAccTaxon
from TaxonGraph import TaxonGraph
import pickle
import pandas as pd

class Multiaccs():
    def __init__(self, path_to_multiacc_file, path_to_multiaccs_taxids_file):
        self.path_to_multiaccs_taxids_file = Path(path_to_multiaccs_taxids_file)
        if not self.path_to_multiaccs_taxids_file.exists():
            self.write_multaccs_to_taxids_file(path_to_multiacc_file)
        self.multiacc_marks = ['WP', 'XP']
        self.multiaccs_to_taxid_dict = self.read_multiaccs_file()

    def read_multiaccs_file(self):
        multiaccs_to_taxid_dict = {}
        with open(self.path_to_multiaccs_taxids_file) as multiacc_tax:
            for line in multiacc_tax:
                fields = line.split()
                multiaccs_to_taxid_dict[fields[0]] = fields[1:]
        return multiaccs_to_taxid_dict

    def is_multiacc(self, acc, multi_acc_dict):
        return True if acc in multi_acc_dict.keys() else False

    def get_multiacc(self, acc):
        accs = acc.split()
        acc_marks = [acc for acc in accs if acc.split('_')[0] in self.multiacc_marks]
        if len(acc_marks) == 1:
            return acc_marks[0]
        else:
            print(accs)

    def get_taxids_of_multiaccs(self, acc):
        taxid_list = self.multiaccs_to_taxid_dict(acc)
        return taxid_list

    def write_multaccs_to_taxids_file(self, path_to_multiacc_file):
        accessions = []
        all_marks = set()
        with open(path_to_multiacc_file, 'r') as multiacc:
            for line in multiacc:
                accessions = accessions.extend(line.split())
        accessions = set(accessions)
       # path_tax2proteome_db = '/home/jules/Documents/databases/databases_tax2proteome/'
        # acc2taxon = ReadAccTaxon(path_tax2proteome_db, 'ncbi')
        # acc_taxon_dict = acc2taxon.read_acc2tax(accessions)
        with open(path_to_multiacc_file, 'r') as multiacc:
            with open(self.path_to_multiaccs_taxids_file, 'r') as multiacc_tax:
                for line in multiacc:
                    fields = line.split()
                    multiacc = fields[0]
                    all_marks.add(multiacc.split('_')[0])
                   # taxa = {acc_taxon_dict[acc] for acc in fields}
                    accessions = accessions.union(set(line.split()))
                    multiacc_tax.write(multiacc + '\t' + ('\t').join(accessions))


class HelperMethod():

    @staticmethod
    def load_taxa_graph(path_to_taxdump):
        """
        # Try load pre-builded taxonomy graph or built taxonomy graph now
        :param options: user input options
        :return: TaxonGraph object
        """

        if not (path_to_taxdump.parents[0] / 'taxon_graph_results').is_file():
            taxon_graph = TaxonGraph()
            print("Start building taxon graph.")
            taxon_graph.create_graph(str(path_to_taxdump))
            print("Taxon graph successfully build.")
            # save TaxonGraph to harddrive:
            try:
                with open(str(path_to_taxdump.parents[0] / 'taxon_graph_results'), 'wb') as handle:
                    pickle.dump(taxon_graph, handle, protocol=pickle.HIGHEST_PROTOCOL)
                    print('Safe taxon graph to location: %s' % str(
                        path_to_taxdump.parents[0] / 'taxon_graph_results'))
            except FileNotFoundError:
                print('Error open tax_graph.')
                exit(1)
        # load Taxon Graph
        else:
            try:
                print('Load taxon graph from harddrive.')
                with open(str(path_to_taxdump.parents[0] / 'taxon_graph_results'), 'rb') as handle:
                    taxon_graph = pickle.load(handle)
            except UnicodeDecodeError or EOFError:
                print(
                    "Failed opening path to taxon graph / taxon_graph is corrupted. Delete %s file."
                    % str(path_to_taxdump.parents[0] / 'taxon_graph'))
                exit(1)
        return taxon_graph

    @staticmethod
    def get_taxid_specific_spectra(df, columnname, taxid, level1, level2, taxongraph):
        taxid_level1 = taxongraph.find_level_up(taxid, level1)
        taxid_level2 = taxongraph.find_level_up(taxid, level2)
        df = df[[columnname, f'taxID_{level1}', f'taxID_{level2}']]
        df = df[taxid_level1 in f'taxID_{level1}' | taxid_level2 in f'taxID_{level1}']
        return df[[df[columnname] == taxid_level]]

    @staticmethod
    def create_df_with_all_spectra_and_dfs_to_compare(all_spectra_list, df1, spectra_column_1, df2, spectra_column_2):
        df_all_spectra = pd.DataFrame(all_spectra_list, columns=['SpectraID'])
        df_with_all_spectra_and_df1 = pd.merge(df_all_spectra, df1, how="outer", left_on='SpectraID', right_on=spectra_column_1)
        df_with_all_spectra_and_df1_df2 = pd.merge(df_with_all_spectra_and_df1, df2, how="outer", left_on='SpectraID',
                                                   right_on=spectra_column_2)
        return df_with_all_spectra_and_df1_df2

    @staticmethod
    def create_df_with_all_reference_spectra(reference_df, reduced_result_df, spectra_column_1, spectra_column_2):
        df_with_all_reference_spectra_and_merged_results = pd.merge(reference_df, reduced_result_df, how="left",
                                                                    left_on=spectra_column_1, right_on=spectra_column_2)
        return df_with_all_reference_spectra_and_merged_results

    @staticmethod
    def get_rows_with_different_taxa(tax_column_1, tax_column_2, taxon_1, taxon_2):
        true_false_list=[]
        for tax_set_1, tax_set_2 in zip(tax_column_1, tax_column_2):
            if not pd.isna(tax_set_1) and not pd.isna(tax_set_2):
                true_false_list.append(taxon_1 in tax_set_1 and taxon_2 in tax_set_2)
            else:
                true_false_list.append(True)
        return true_false_list

    @staticmethod
    def get_difference_between_two_df_for_one_taxon(all_spectra_list, df1, spectra_column_1, column_of_interest_1, df2,
                                      spectra_column_2, column_of_interest_2, taxon, level1, level2, taxon_graph):
        df_with_all_spectra_and_df1_df2 = HelperMethod.create_df_with_all_spectra_and_dfs_to_compare(all_spectra_list, df1, spectra_column_1, df2, spectra_column_2)
        # remove empty lines
        df_with_all_spectra_and_df1_df2 = df_with_all_spectra_and_df1_df2[df_with_all_spectra_and_df1_df2.Hyperscore_x.notna() & df_with_all_spectra_and_df1_df2.Hyperscore_y.notna()]
        taxon_1 = taxon_graph.find_level_up(taxon, level1)
        taxon_2 = taxon_graph.find_level_up(taxon, level2)
        df_difference=df_with_all_spectra_and_df1_df2[HelperMethod.get_rows_with_different_taxa(df_with_all_spectra_and_df1_df2[column_of_interest_1].tolist(), df_with_all_spectra_and_df1_df2[column_of_interest_2].tolist(), taxon_1, taxon_2)]
        return df_difference
