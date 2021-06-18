from pathlib import Path
import pandas as pd
import argparse
import pickle
from collections import defaultdict
from handling_acc_files import Multiaccs

class PSM_FDR:
    def __init__(self, path_to_file):
        """
        :param path_to_file: path to xtandem output .xml converted to tsv
        one line one header with all accessions seperated by \t
        """
        self.path_to_file = path_to_file
        self.sorted_xtandem_df = None
        self.spectra_acc_dict = defaultdict(list)
        self.spectra_sequence_dict = {}
        self.fdr_pos = 0
        self.number_psms = 0
        self.decoys = 0

    def flatten_set(self, s):
        flatten_set = {item for sublist in s for item in sublist}
        print(s)
        print(flatten_set)
        return flatten_set

    def flatten_sets_to_one_set(self, *sets):
        result_set = set()
        for s in sets:
            result_set.union(s)
        return result_set


    def add_level_specific_taxid_column_and_group(self, taxon_graph, level):
        self.sorted_xtandem_df[f'taxID_{level}'] = self.sorted_xtandem_df.apply(lambda row:
                                                taxon_graph.find_level_up(int(row['taxID']), level)
                                                if row['taxID'] != 'DECOY' and row['taxID'] != 0 else 'DECOY',
                                                                                axis=1)
        reduced_df = self.sorted_xtandem_df.groupby(["#SpecFile", 'Title', 'Peptide', 'Hyperscore'], as_index=False).agg(
            {'Protein': lambda acc: set(acc), 'EValue': lambda x: set(list(x)), 'decoy': lambda x: set(x),
             'taxID': lambda taxid: set(taxid), f'taxID_{level}': lambda x: set(x)})
        return reduced_df


    def create_PSM_dataframe_for_uniprot_accs(self, acc2tax_dict, taxon_graph, level):
        self.sorted_xtandem_df['Protein'] = self.sorted_xtandem_df.apply(lambda row: row['Protein'].split()[0], axis=1)
        self.sorted_xtandem_df['taxID'] = self.sorted_xtandem_df.apply(lambda row:
                                                                       acc2tax_dict[row['Protein'].split('|')[1]]
                                                                       if row['Protein'].split('|')[1] in acc2tax_dict else 'DECOY',
                                                                       axis=1)
        reduced_df = self.add_level_specific_taxid_column_and_group(taxon_graph, level)
        return reduced_df

    def create_PSM_dataframe_for_ncbi_accs(self, acc2tax_dict, taxon_graph, level):
        pd.set_option("display.max_rows", None, "display.max_columns", None)
        # Multiaccs in Protein column: 'PNW76085.1 BAB64417.1 BAB64413.1 XP_001693987.1'
        self.sorted_xtandem_df['Protein'] = self.sorted_xtandem_df['Protein'].apply(lambda acc: acc.split('|')[2]
                                                                                    .split()[0] if acc.startswith('generic') else acc)
        self.sorted_xtandem_df['taxID'] = self.sorted_xtandem_df['Protein'].apply(lambda acc:
                                                                                  acc2tax_dict[acc] if acc in acc2tax_dict.keys() else {'DECOY'} )
        # 'Protein': generic|AC:4260012|KKY45073.1 WP_046834534.1, generic|AC:6172351|WP_086660554.1 ...
        # NCBI: 'taxID': {83334}, {562}
        print(self.sorted_xtandem_df['taxID'][0:5])
        self.sorted_xtandem_df[f'taxID_{level}'] = self.sorted_xtandem_df['taxID'].apply(lambda taxID_set:
                                                                                         {taxon_graph.find_level_up(int(taxID), level)
                                                                                          if taxID != 'DECOY' and taxID != 0 else 'DECOY'
                                                                                          for taxID in taxID_set})
        print(self.sorted_xtandem_df[f'taxID_{level}'][0:5])
        print(self.sorted_xtandem_df.head())
        reduced_df = self.sorted_xtandem_df.groupby(["#SpecFile", 'Title', 'Peptide', 'Hyperscore'], as_index=False).agg(
            {'Protein': lambda acc: set(acc), 'EValue': lambda x: set(list(x)), 'decoy': lambda decoy: set(decoy),
             'taxID': lambda taxid_sets: self.flatten_set(taxid_sets),
             f'taxID_{level}': lambda taxid_sets: self.flatten_set(taxid_sets)})
        print(reduced_df.head())
        return reduced_df

    def create_PSM_dataframe_for_custom_accs(self, acc2tax_dict, decoy_tag, taxon_graph, level):
        self.sorted_xtandem_df['taxID'] = self.sorted_xtandem_df.apply(lambda row:
                                                                       acc2tax_dict[row['Protein'].strip()]
                                                                       if not decoy_tag in row['Protein'] else 'DECOY', axis=1)
        reduced_df = self.add_level_specific_taxid_column_and_group(taxon_graph, level)
        return reduced_df

    def create_PSM_dataframe(self, decoy_tag, db_type, level=None, taxon_graph=None, acc2tax_dict=None, multiaccs_dict=None):
        """

        :return: Hyperscore sorted dataframe, Protein with highest Hyperscore first, Title = spectra_file_spectraID
        """
        xtandem_df = pd.read_csv(str(self.path_to_file), delimiter='\t')
        xtandem_df['Protein'] = xtandem_df['Protein'].apply(lambda acc: acc.strip())
        print('Sorting panda dataframe by hyperscore. Find taxa, and level taxa')
        self.sorted_xtandem_df = xtandem_df.sort_values(by=['Hyperscore', 'Title'], ascending=False).reset_index(drop=True)
        self.sorted_xtandem_df['decoy'] = self.sorted_xtandem_df.apply(lambda row: True if decoy_tag in row['Protein'] else False, axis=1)

        if db_type == 'uniprot':
            reduced_df = self.create_PSM_dataframe_for_uniprot_accs(acc2tax_dict, taxon_graph, level)
        elif db_type == 'ncbi':
            reduced_df = self.create_PSM_dataframe_for_ncbi_accs(acc2tax_dict, taxon_graph, level)
        elif db_type == 'custom':
            reduced_df = self.create_PSM_dataframe_for_custom_accs(acc2tax_dict, decoy_tag, taxon_graph, level)

        # change spectra Title
        reduced_df['Title'] = reduced_df['Title'].apply(lambda row: row.split(' File')[0])
        # sort by Hyperscore
        reduced_df = reduced_df.sort_values(by=['Hyperscore', 'Title'], ascending=False).reset_index(drop=True)
        print(f"writing data frame to {self.path_to_file+'.tsv'}... ")
        reduced_df.to_csv(self.path_to_file+'.tsv', sep='\t')
        print(f'entries in df: {len(self.sorted_xtandem_df)} '
              f'decoys in df: {len(self.sorted_xtandem_df[self.sorted_xtandem_df.decoy==True])} '
              f'hits in df: {len(self.sorted_xtandem_df[self.sorted_xtandem_df.decoy==False])}')
        return reduced_df

    @staticmethod
    def determine_FDR_position(sorted_xtandem_df, fdr):
        """
        :param decoy_tag: tag of decoy entries, for example REVERSED
        :param fdr: false discovery rate, for example 0.01
        """
        repeatedly_identified_spectra = set()
        number_multiple_identified_spectra = 0
        hits = decoy = 0
        FDR_position_not_set = True
        title_set = set()
        spectra_header = [(x, y) for x, y in zip(sorted_xtandem_df['Title'], sorted_xtandem_df['decoy'])]
        for elem in spectra_header:
            if elem[0] in title_set:
                repeatedly_identified_spectra.add(elem[0])
                number_multiple_identified_spectra += 1
                continue
            else:
                title_set.add(elem[0])
            if elem[1]:
                decoy += 1
            else:
                hits += 1
            if decoy / (hits + decoy) > fdr and FDR_position_not_set:
                fdr_pos = hits + decoy - 1 + number_multiple_identified_spectra
                number_psms = hits
                decoys = decoy - 1
                FDR_position_not_set = False
                break
            if decoy / (hits + decoy) <= fdr:
                continue
                #and not FDR_position_not_set:
                # self.fdr_pos = hits + decoy + number_multiple_identified_spectra
                # self.number_psms = hits
                # self.decoys = decoy - 1
                # FDR_position_not_set = True
        print('Number of PSMs: %d' % number_psms)
        print('Number of decoys: %d' % decoys)
        print(f"double identified spetra {number_multiple_identified_spectra}")
        print('Position FDR border/Number of PSMs: %d' % fdr_pos)
        return fdr_pos, number_psms, decoys

    # for every spectrum (key) all by xtandem found identifications = accessions (value)
    def create_spectrum_acc_dict(self, db):
        """
        :param db: db-format, 'uniprot' or 'ncbi' or 'custom'
        """
        # 'Title' = spectrum ID, 'Protein' = protein acc
        spectra_header_psm = [(x, y, z) for x, y, z in zip(self.sorted_xtandem_df['Title'], self.sorted_xtandem_df['Protein'], self.sorted_xtandem_df['decoy'])][0:self.fdr_pos]

        for elem in spectra_header_psm:
            if elem[2]:
                self.spectra_acc_dict[elem[0]].append('REVERSED')
            else:
                if db == 'uniprot':
                    self.spectra_acc_dict[elem[0]].append(elem[1].split('|')[1])
                elif db == 'ncbi':
                    self.spectra_acc_dict[elem[0]].append('|'.join(elem[1].split('|')[2:]).split(' ')[0])
                elif db == 'custom':
                    self.spectra_acc_dict[elem[0]].append(elem[1])

        spectra_sequence_psm = [(x, y) for x, y in
                              zip(self.sorted_xtandem_df['Title'], self.sorted_xtandem_df['Peptide'])][0:self.fdr_pos]

        for elem in spectra_sequence_psm:
            if elem[0] in self.spectra_sequence_dict:
                continue
            else:
                if 'REVERSED' in elem[1]:
                    self.spectra_sequence_dict[elem[0]] = ['REVERSED']
                else:
                    self.spectra_sequence_dict[elem[0]] = elem[1]


