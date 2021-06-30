import pandas as pd
from pathlib import Path
import argparse
from handling_acc_files import HelperMethod
from create_reference_from_tsv_and_pepxml import ReferenceWriter
from collections import defaultdict
from create_PSM_df import PSM_FDR

class SensitivityAndSpecificity():

    def __init__(self, df_reference, df_result, spectra_file, level, fdr):
        """
        :param df_reference: X-tandem search result with Kleiner database on Run1_2000.mgf (the reduced df)
        :param df_result: X-tandem search result with Tax2Proteome database (the reduced df)
        :param spectra_file: complete Run1_2000.mgf file
        :param level: level of tax2proteome generated database
        """
        self.level = level
        self.result_df = df_result
        self.reference_df = self.rename_reference_df_columns(df_reference)
        self.all_spectra_set = self.get_all_spectra_IDs(spectra_file)
        self.fdr_pos_result, self.number_psm_result, self.number_decoy_result, self.double_spectra_result  = PSM_FDR.determine_FDR_position(self.result_df, fdr, True)
        self.fdr_pos_reference, self.number_psm_reference, self.number_decoy_reference, self.double_spectra_reference = PSM_FDR.determine_FDR_position(self.reference_df, fdr, True, 'Ref_decoy')
        self.decoy_list = [{'DECOY'}, {'DECOY/CRAP'}, {0}]

    def get_all_spectra_IDs(self, ident_file):
        all_spec_IDs = set()
        with open(ident_file, 'r') as ident_file:
            for line in ident_file:
                if line.startswith('TITLE'):
                    all_spec_IDs.add(line.split()[0].split('TITLE=')[1])
        return all_spec_IDs

    def rename_reference_df_columns(self, reference_df):
        reference_df = reference_df[['Title', 'Peptide', 'Hyperscore', 'Protein', 'decoy', 'taxID', f'taxID_{self.level}']]
        reference_df.columns = ['Title', 'Ref_Peptide', 'Ref_Hyperscore', 'Ref_ProteinAcc', 'Ref_decoy', 'Ref_taxID_DB', f'Ref_taxID_{self.level}']
        return reference_df

    def get_rows_with_no_decoy_or_nan_in_reference_and_result(self, Ref_taxID_level_column, taxID_level_column):
        true_false_list = []
        for ref_taxID_set, taxID_set in zip(Ref_taxID_level_column, taxID_level_column):
            if ref_taxID_set in self.decoy_list and taxID_set in self.decoy_list:
                true_false_list.append(False)
            # type(taxID_set) != set = is nan
            elif (ref_taxID_set in self.decoy_list and type(taxID_set) != set) or (type(ref_taxID_set) != set and taxID_set in self.decoy_list):
                true_false_list.append(False)
            else:
                true_false_list.append(True)
        return true_false_list

    def get_rows_with_nan_or_decoy_in_reference_and_result(self, ref_taxID_column, taxID_column):
        true_false_list = []
        for ref_taxID_set, taxID_set in zip(ref_taxID_column, taxID_column):
            if type(ref_taxID_set) != set and type(taxID_set) != set:
                true_false_list.append(True)
            elif (ref_taxID_set == {'DECOY'} and type(taxID_set) != set) or (type(ref_taxID_set) != set and taxID_set == {'DECOY'}):
                true_false_list.append(True)
            else:
                true_false_list.append(False)
        return true_false_list

    def get_row_with_reference_nan_or_decoy_and_result_identified(self, ref_taxID_column, taxID_column):
        true_false_list = []
        for ref_taxID_set, taxID_set in zip(ref_taxID_column, taxID_column):
            if type(ref_taxID_set) != set and type(taxID_set) == set:
                if taxID_set not in self.decoy_list:
                    true_false_list.append(True)
                else:
                    true_false_list.append(False)
            elif ref_taxID_set in self.decoy_list and taxID_set not in self.decoy_list:
                true_false_list.append(True)
            else:
                true_false_list.append(False)
        return true_false_list

    def check_for_TP(self, taxid_set, taxid_ref_set):
        # TP: if all taxa from ref set contained in result set
        # ignore Decoy crap entries from result_reduced (not contained in reference)
        decoy_set = {'DECOY', 'DECOY/CRAP', 0}
        taxid_set = taxid_set.difference(decoy_set)
        if len(taxid_set) == 0 and pd.isna(taxid_ref_set): #only decoy entries
            return False
        return taxid_ref_set.issubset(taxid_ref_set)

    def check_for_FP(self, taxid_set, taxid_ref_set):
        # FP: if all taxa from result_df not in reference_set
        if len(taxid_set.intersection(taxid_ref_set)) == 0:
            return True
        else:
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

    def remove_all_decoy_and_nan_rows(self, df_with_all_reference_spectra_and_merged_results_in_fdr):
        return df_with_all_reference_spectra_and_merged_results_in_fdr[
            self.get_rows_with_no_decoy_or_nan_in_reference_and_result(df_with_all_reference_spectra_and_merged_results_in_fdr[f'Ref_taxID_{self.level}'],
                                                                       df_with_all_reference_spectra_and_merged_results_in_fdr[f'taxID_{self.level}'])
        ]

    def get_df_with_all_unidentified_spectra_in_reference_and_result(self):
        df_with_all_spectra_and_reference_and_results =  HelperMethod.create_df_with_all_spectra_and_dfs_to_compare(
            self.all_spectra_set, self.reference_df, 'Title', self.result_df, 'Title')
        df_with_all_unidentified_spectra = df_with_all_spectra_and_reference_and_results[
            self.get_rows_with_nan_or_decoy_in_reference_and_result(df_with_all_spectra_and_reference_and_results.Ref_taxID_DB,
                                                                    df_with_all_spectra_and_reference_and_results.taxID)]
        return df_with_all_unidentified_spectra

    def get_true_positive_and_true_negative(self, df_with_all_unidentified_spectra, df_with_all_reference_spectra_and_merged_results):
        pd.set_option("display.max_rows", None, "display.max_columns", None)
        df_taxid = df_with_all_reference_spectra_and_merged_results[['Title','taxID','Ref_taxID_DB']]
        df_taxid_level = df_with_all_reference_spectra_and_merged_results[['Title', f'taxID_{self.level}',f'Ref_taxID_{self.level}']]
        # both nan or one nan one DECOY
        df_TN = df_with_all_unidentified_spectra
        # print(df_TN.head())
        TN=len(set(df_TN.SpectraID.tolist()))

        df_TP = df_taxid_level[self.check_taxid_in_reference(df_taxid_level[f'taxID_{self.level}'].tolist(),
                                                             df_taxid_level[f'Ref_taxID_{self.level}'].tolist(), is_FP=False)]
        TP = len(set(df_TP.Title.tolist()))
        # FN: not identified in result, but identified in referernce
        df_FN = df_taxid[(df_taxid.taxID.isna() ) & df_taxid.Ref_taxID_DB.notna()]
        df_FN = df_FN[df_FN.Ref_taxID_DB != {'DECOY/CRAP'}]
        FN=len(set(df_FN.Title.tolist()))


        df_FP = df_taxid_level[self.check_taxid_in_reference(df_taxid_level[f'taxID_{self.level}'].tolist(),
                                                              df_taxid_level[f'Ref_taxID_{self.level}'].tolist(), is_FP=True)]
        FP = len(set(df_FP.Title.tolist()))
        df_s = df_with_all_reference_spectra_and_merged_results[df_with_all_reference_spectra_and_merged_results.taxID != {'DECOY/CRAP'}]
        print(f"TP: {TP}, FP: {FP}, TN: {TN}, FN: {FN}, number of all spectra without decoy/crap spectra: "
              f"{len(set(df_s.Title))}, number of TP+TN+FP+FN: {TP+FN+FP+TN}")
        return TP, FP, TN, FN

    @staticmethod
    def calculate_sensitivity(TP, FN):
        return TP/(TP+FN) * 100

    @staticmethod
    def calculate_specificity(FP, TN):
        return TN/(TN + FP) * 100

    def get_df_with_identified_spectra_in_result_df_but_not_in_reference_df(self):
        df_with_all_result_spectra_and_merged_reference_in_fdr = pd.merge(self.result_df[0:self.fdr_pos_result],
                                                                          self.reference_df[0:self.fdr_pos_reference],
                                                                          how="left", on='Title')
        df_with_only_in_result_identified_spectra = df_with_all_result_spectra_and_merged_reference_in_fdr[
            self.get_row_with_reference_nan_or_decoy_and_result_identified(df_with_all_result_spectra_and_merged_reference_in_fdr.Ref_taxID_DB,
                                                                           df_with_all_result_spectra_and_merged_reference_in_fdr.taxID)]
        return df_with_only_in_result_identified_spectra