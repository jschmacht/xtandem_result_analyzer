import pandas as pd
from pathlib import Path
import argparse
from handling_acc_files import HelperMethod
from create_reference_from_tsv_and_pepxml import ReferenceWriter
from collections import defaultdict
from create_PSM_df import PSM_FDR

class DeterminatorSpecificitySensitivityFirstTry():

    def __init__(self, level, fdr_applied_df, reference_df, spectra_file):
        """
        :param fdr_applied_df:
        :param refernce_df: column names =
        :param spectra_file: 'Run1_U1_2000ng.mgf'
        """
        self.tax_level = ['species', 'genus', 'family', 'order']
        self.level = level
        self.result_df = fdr_applied_df[['Title', 'Peptide', 'Hyperscore', 'Protein', 'decoy', 'taxID', f'taxID_{level}']]
        self.reference_df = reference_df[['SpectraID', 'Ref_Peptide', 'Ref_Hyperscore', 'Ref_ProteinAcc', 'Ref_decoy', 'Ref_taxID_DB', f'Ref_taxID_{level}']]
        self.all_spectra_set = self.get_all_spectra_IDs(spectra_file)

    def create_df_with_all_spectra_reference_and_result_taxa(self, path_to_out):
        """
        :param path_to_out:
        :return: df_with_all_spectra_and_reference_and_results: in result_reduced.tsv file multiple entries for spectra possible
        (e.g. different Peptides, different Hyperscores (only one Hyperscore per cell for sorting) in reference_file only one row per spectra
        in df_with_all_spectra_and_reference_and_results in Ref columns same entries, if spectra occures multiple times
        """
        df_with_all_spectra_and_reference_and_results =  HelperMethod.create_df_with_all_spectra_and_dfs_to_compare(
            self.all_spectra_set, self.result_df, 'Title', self.reference_df, 'SpectraID')
        print(f"write df_with_all_spectra_and_reference_and_results {path_to_out}... ")
        df_with_all_spectra_and_reference_and_results.to_csv(str(path_to_out), sep='\t')
        return df_with_all_spectra_and_reference_and_results

    def create_df_with_all_reference_spectra_and_matching_tax2proteome_results(self):
        df_with_all_reference_spectra = HelperMethod.create_df_with_all_reference_spectra
        return df_with_all_reference_spectra

    def calculate_sensitivity(self, TP, FN):
        return TP/(TP+FN) * 100

    def calculate_specificity(self, FP, TN):
        return TN/(TN + FP) * 100

    def calculate_sensitivity_and_specificity(self, path_to_out):
        df_with_all_spectra_and_reference_and_results = self.create_df_with_all_spectra_reference_and_result_taxa(path_to_out)
        print('calculate TP, FP, TN, FN')
        TP, FP, TN, FN = self.get_true_positive_and_true_negative(df_with_all_spectra_and_reference_and_results)
        sensitivity = self.calculate_sensitivity(TP, FN)
        specificity = self.calculate_specificity(FP, TN)
        print(f'sensitivity: {sensitivity}%, specificity: {specificity}%')

    def load_ref_file(self, ref_file, level):
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

    def get_all_spectra_IDs(self, ident_file):
        all_spec_IDs = set()
        with open(ident_file, 'r') as ident_file:
            for line in ident_file:
                if line.startswith('TITLE'):
                    all_spec_IDs.add(line.split()[0].split('TITLE=')[1])
        return all_spec_IDs

    def check_for_TP(self, taxid_set, taxid_ref_set):
        # ignore Decoy crap entries from result_reduced (not contained in reference)
        decoy_set = {'DECOY', 'DECOY/CRAP', 0}
        taxid_set = taxid_set.difference(decoy_set)
        if len(taxid_set) == 0 and pd.isna(taxid_ref_set): #only decoy entries
            return False
        return taxid_set.issubset(taxid_ref_set)

    def check_for_FP(self, taxid_set, taxid_ref_set):
        for taxid in taxid_set:
            # ignore Decoy crap entries from result_reduced (not contained in reference)
            if taxid == 'DECOY/CRAP' or taxid == 'DECOY':
                continue
            else:
                if taxid not in taxid_ref_set:
                    return True
        return False

    def compare_tax_sets(self, taxid_set, taxid_ref_set, is_FP):
        if not pd.isna(taxid_set) and not pd.isna(taxid_ref_set):
            if is_FP:
                return self.check_for_FP(taxid_set, taxid_ref_set)
            else:
                return self.check_for_TP(taxid_set, taxid_ref_set)
        return False

    def check_taxid_in_reference(self, taxid_level_column, taxid_level_ref_column, is_FP):
        true_false_list = []
        for taxid_set, taxid_ref_set in zip(taxid_level_column, taxid_level_ref_column):
            true_false_list.append(self.compare_tax_sets(taxid_set, taxid_ref_set, is_FP))
        return true_false_list

    def get_true_positive_and_true_negative(self, df_with_all_spectra_and_reference_and_results):
        pd.set_option("display.max_rows", None, "display.max_columns", None)
        df_taxid = df_with_all_spectra_and_reference_and_results[['SpectraID','taxID','Ref_taxID_DB']]
        df_taxid_level = df_with_all_spectra_and_reference_and_results[['SpectraID', f'taxID_{self.level}',f'Ref_taxID_{self.level}']]

        df_TN = df_taxid[df_taxid.taxID != {'DECOY/CRAP'} & df_taxid.Ref_taxID_DB.isna()]
        df_TN = df_TN[df_TN.taxID.isna() & df_TN.Ref_taxID_DB.isna()]
        TN=len(set(df_TN.SpectraID.tolist()))
        df_TP = df_taxid_level[self.check_taxid_in_reference(df_taxid_level[f'taxID_{self.level}'].tolist(), df_taxid_level[f'Ref_taxID_{self.level}'].tolist(), is_FP=False)]
        TP = len(set(df_TP.SpectraID.tolist()))

        df_FN = df_taxid[(df_taxid.taxID.notna() ) & df_taxid.Ref_taxID_DB.isna()]
        df_FN = df_FN[df_FN.taxID != {'DECOY/CRAP'}]
        FN=len(set(df_FN.SpectraID.tolist()))


        df_FP = df_taxid_level[self.check_taxid_in_reference(df_taxid_level[f'taxID_{self.level}'].tolist(), df_taxid_level[f'Ref_taxID_{self.level}'].tolist(), is_FP=True)]
        FP = len(set(df_FP.SpectraID.tolist()))
        df_s = df_with_all_spectra_and_reference_and_results[df_with_all_spectra_and_reference_and_results.taxID != {'DECOY/CRAP'}]
        print(f"TP: {TP}, FP: {FP}, TN: {TN}, FN: {FN}, number of all spectra without decoy/crap spectra: "
              f"{len(set(df_s.SpectraID))}, number of TP+TN+FP+FN: {TP+FN+FP+TN}")
        return TP, FP, TN, FN


def main():
    parser = argparse.ArgumentParser(description='Read xtandem output .tsv')
    parser.add_argument('-i', '--input', dest='input', default=None, help='xtandem reduced.tsv with columns:'
                                                                          '[spectra_file, Spectrum ID, Peptide, Protein Accession + Description, Hyperscore, Evalue,'
                                                                          'or pep.xml file for create_ref')
    parser.add_argument('-r', '--reference', dest='reference', default=None, help='Path to reference tsv')
    parser.add_argument('-g', '--tax_graph', dest='tax_graph', help='Path to taxdump.tar.gz')
    options = parser.parse_args()

    path_to_taxdump = Path('/home/jules/Documents/databases/databases_tax2proteome/taxdump.tar.gz')
    taxon_graph = HelperMethod.load_taxa_graph(path_to_taxdump)
    path_to_input = Path('/home/jules/Documents/Tax2Proteome/benchmarking/results_searchgui_kleiner_db/Run1_U1_2000ng.t.xml_reduced.tsv')
    path_to_reference = Path('/home/jules/Documents/Tax2Proteome/benchmarking/results_searchgui_kleiner_db/Run1_U1_2000ng.t.xml_reduced_reference.tsv')
    # create reference file
    reference_df = ReferenceWriter.read_csv_with_generic_function(path_to_input,
                                                                  ['Protein', 'Hyperscore', 'decoy', 'taxID', f'taxID_species'])
    reader = ReferenceWriter(reference_df, path_to_reference, taxon_graph)
    reader.write_result_spectrum_reference_file()


    #  TP analysis first try
    #  path_to_all_info_tsv = path_to_identification_file.parent.joinpath(path_to_identification_file.stem + '_' + path_to_reference.stem + '.tsv')
    #  print(path_to_all_info_tsv)
    #  result_df = ReferenceWriter.read_csv_with_generic_function(options.input,
    #                                         ['Protein', 'decoy', 'taxID', f'taxID_{options.level}'])
    #  reference_df = ReferenceWriter.read_csv_with_generic_function(options.reference,
    #                                          ['Ref_Peptide', 'Ref_ProteinAcc', 'Ref_Hyperscore', 'Ref_decoy','Ref_taxID_DB', f'Ref_taxID_{options.level}'])
    #  determinator = DeterminatorSpecificitySensitivity(options.level, fdr_applied_df, fdr_applied_reference_df, options.spectra_file)
    #  determinator.calculate_sensitivity_and_specificity(path_to_all_info_tsv)
        #determinator.get_true_positive_and_true_negative(options.level, fdr_applied_df, reference_df, options.spectra_file)

if __name__ == '__main__':
    main()
