import numpy as np
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
        self.fdr_pos_result, self.number_psm_result, self.number_decoy_result, self.double_spectra_result,\
        self.score_last_item_result  = PSM_FDR.determine_FDR_position(self.result_df, fdr, True)
        self.fdr_pos_reference, self.number_psm_reference, self.number_decoy_reference, self.double_spectra_reference, \
        self.score_last_item_reference = PSM_FDR.determine_FDR_position(self.reference_df, fdr, True, 'Ref_decoy')
        self.decoy_list = [{'DECOY'}, {0}]

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
            if (ref_taxID_set in self.decoy_list or type(ref_taxID_set) != set) and (taxID_set in self.decoy_list or type(taxID_set) != set):
                true_false_list.append(False)
            else:
                true_false_list.append(True)
        return true_false_list

    def get_rows_with_nan_or_decoy_in_reference_and_result(self, ref_taxID_column, taxID_column):
        true_false_list = []
        for ref_taxID_set, taxID_set in zip(ref_taxID_column, taxID_column):
            # NaN Rows
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


    def check_for_TP(self, taxid_level_column, taxid_level_ref_column):
        decoy_set = {'DECOY', 0}
        true_false_list = []
        for taxid_set, taxid_ref_set in zip(taxid_level_column, taxid_level_ref_column):
            # taxid_set.difference(decoy_set) remove decoy entries
            # both tax sets empty or with decoy
            if len(taxid_set.difference(decoy_set)) == 0 and len(taxid_ref_set.difference(decoy_set)) == 0:
                true_false_list.append(False)
            else:
                # at least on taxon of result contained in refernce
                true_false_list.append(len(taxid_set.difference(decoy_set).intersection(taxid_ref_set)) > 0)
        return true_false_list

    def check_for_FP(self, taxid_level_column, taxid_level_ref_column):
        decoy_set = {'DECOY', 0}
        true_false_list = []
        for taxid_set, taxid_ref_set in zip(taxid_level_column, taxid_level_ref_column):
            # remove decoy entries
            taxid_set = taxid_set.difference(decoy_set)
            if len(taxid_set) > 0 and (len(taxid_ref_set) == 0 or taxid_ref_set in self.decoy_list) :
                is_FP = True
            elif len(taxid_set) == 0:
                is_FP = False
            else:
                # no taxon of result taxa set contained in refernce
                is_FP = len(taxid_set.intersection(taxid_ref_set)) == 0
            true_false_list.append(is_FP)
        return true_false_list

    def check_for_FN(self, taxID_column, Ref_taxID_DB_column):
        decoy_set = {'DECOY', 0}
        true_false_list = []
        for taxid_set, taxid_ref_set in zip(taxID_column, Ref_taxID_DB_column):
            if len(taxid_set.difference(decoy_set)) == 0 and len(taxid_ref_set.difference(decoy_set)) > 0 :
                 true_false_list.append(True)
            else:
                true_false_list.append(False)
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

    def replace_nan_by_empty_set(self, df_all):
        df_all.loc[df_all['taxID'].isnull(),['taxID']] = df_all.loc[df_all['taxID'].isnull(), 'taxID'].apply(lambda x: [])
        df_all.loc[df_all['Ref_taxID_DB'].isnull(),['Ref_taxID_DB']] = df_all.loc[df_all['Ref_taxID_DB'].isnull(), 'Ref_taxID_DB'].apply(lambda x: [])
        df_all.loc[df_all[f'taxID_{self.level}'].isnull(),[f'taxID_{self.level}']] = df_all.loc[df_all[f'taxID_{self.level}'].isnull(), f'taxID_{self.level}'].apply(lambda x: [])
        df_all.loc[df_all[f'Ref_taxID_{self.level}'].isnull(),[f'Ref_taxID_{self.level}']] = df_all.loc[df_all[f'Ref_taxID_{self.level}'].isnull(), f'Ref_taxID_{self.level}'].apply(lambda x: [])
        return df_all

    def get_true_positive_and_true_negative(self, df_with_all_unidentified_spectra, df_with_all_reference_spectra_and_merged_results):
        df_with_all_reference_spectra_and_merged_results.to_csv('/home/jules/Documents/Tax2Proteome/benchmarking/all_df.tsv', sep='\t')

        df_with_all_reference_spectra_and_merged_results = self.replace_nan_by_empty_set(df_with_all_reference_spectra_and_merged_results)

        df_taxid = df_with_all_reference_spectra_and_merged_results[['Title','taxID','Ref_taxID_DB']]
        df_taxid = df_taxid.groupby(['Title'], as_index=False).agg({f'taxID': lambda x: self.flatten_set(x),
                                                                    f'Ref_taxID_DB': lambda x: self.flatten_set(x)})

        df_taxid_level = df_with_all_reference_spectra_and_merged_results[['Title', f'taxID_{self.level}',f'Ref_taxID_{self.level}']]
        df_taxid_level = df_taxid_level.groupby(['Title'], as_index=False).agg({
            f'taxID_{self.level}': lambda x: self.flatten_set(x), f'Ref_taxID_{self.level}': lambda x: self.flatten_set(x)})
        # both nan or one nan one DECOY
        df_TN = df_with_all_unidentified_spectra
        df_TN.to_csv('/home/jules/Documents/Tax2Proteome/benchmarking/TN.tsv', sep='\t')
        # print(df_TN.head())
        TN=len(set(df_TN.SpectraID.tolist()))

        df_TP = df_taxid_level[self.check_for_TP(df_taxid_level[f'taxID_{self.level}'].tolist(),
                                                             df_taxid_level[f'Ref_taxID_{self.level}'].tolist())]
        df_TP.to_csv('/home/jules/Documents/Tax2Proteome/benchmarking/TP.tsv', sep='\t')
        TP = len(set(df_TP.Title.tolist()))
        # FN: not identified in result, but identified in referernce
        df_FN = df_taxid[self.check_for_FN(df_taxid.taxID, df_taxid.Ref_taxID_DB)]
        df_FN.to_csv('/home/jules/Documents/Tax2Proteome/benchmarking/FN.tsv', sep='\t')
        FN=len(set(df_FN.Title.tolist()))


        df_FP = df_taxid_level[self.check_for_FP(df_taxid_level[f'taxID_{self.level}'].tolist(),
                                                              df_taxid_level[f'Ref_taxID_{self.level}'].tolist())]
        df_FP.to_csv('/home/jules/Documents/Tax2Proteome/benchmarking/FP.tsv', sep='\t')
        FP = len(set(df_FP.Title.tolist()))
        df_s = df_with_all_reference_spectra_and_merged_results[df_with_all_reference_spectra_and_merged_results.taxID != {'DECOY'}]
        print(f"TP: {TP}, FP: {FP}, TN: {TN}, FN: {FN}, number of all spectra without decoy/crap spectra: "
              f"{len(set(df_s.Title))}, number of TP+TN+FP+FN: {TP+FN+FP+TN}")
        return TP, FP, TN, FN

    @staticmethod
    def calculate_sensitivity(TP, FN):
        return TP/(TP+FN) * 100

    @staticmethod
    def calculate_specificity(FP, TN):
        return TN/(TN + FP) * 100

    def flatten_set(self, s):
        return {item for sublist in s for item in sublist}

    def get_df_with_identified_spectra_in_result_df_but_not_in_reference_df(self):
        result_df = self.result_df[0:self.fdr_pos_result].groupby(["#SpecFile", 'Title'], as_index=False).agg(
            {'Peptide': lambda pep: set(pep), 'Hyperscore': lambda score: set(score), 'Protein': lambda acc: self.flatten_set(acc),
              'decoy': lambda x: self.flatten_set(x),
             'taxID': lambda taxid: self.flatten_set(taxid), f'taxID_{self.level}': lambda x: self.flatten_set(x)})
        reference_df = self.reference_df[0:self.fdr_pos_reference].groupby(['Title'], as_index=False).agg(
            {'Ref_Peptide': lambda pep: set(pep), 'Ref_Hyperscore': lambda score: set(score), 'Ref_ProteinAcc': lambda acc: self.flatten_set(acc),
             'Ref_decoy': lambda x: self.flatten_set(x),
             'Ref_taxID_DB': lambda taxid: self.flatten_set(taxid), f'Ref_taxID_{self.level}': lambda x: self.flatten_set(x)})
        df_with_all_result_spectra_and_merged_reference_in_fdr = pd.merge(result_df,
                                                                          reference_df,
                                                                          how="left", on='Title')
        df_with_only_in_result_identified_spectra = df_with_all_result_spectra_and_merged_reference_in_fdr[
            self.get_row_with_reference_nan_or_decoy_and_result_identified(df_with_all_result_spectra_and_merged_reference_in_fdr.Ref_taxID_DB,
                                                                           df_with_all_result_spectra_and_merged_reference_in_fdr.taxID)]
        return df_with_only_in_result_identified_spectra