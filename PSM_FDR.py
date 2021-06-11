from pathlib import Path
import pandas as pd
import argparse
import pickle
from collections import defaultdict

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

    def create_PSM_dataframe(self, decoy_tag, db_type, level=None, taxon_graph=None, acc2tax_dict=None):
        """

        :return: Hyperscore sorted dataframe, Protein with highest Hyperscore first, Title = spectra_file_spectraID
        """
        xtandem_df = pd.read_csv(str(self.path_to_file), delimiter='\t')
        print('Sorting panda dataframe by hyperscore. Find taxa, and level taxa')
        self.sorted_xtandem_df = xtandem_df.sort_values(by=['Hyperscore', 'Title'], ascending=False).reset_index(drop=True)
        self.sorted_xtandem_df['decoy'] = self.sorted_xtandem_df.apply(lambda row: True if decoy_tag in row['Protein'] else False, axis=1)
        print('df with decoy column')
        if db_type == 'uniprot':
            self.sorted_xtandem_df['Protein'] = self.sorted_xtandem_df.apply(lambda row: row['Protein'].split()[0], axis=1)
            self.sorted_xtandem_df['taxID'] = self.sorted_xtandem_df.apply(lambda row:
                                                acc2tax_dict[row['Protein'].split('|')[1]]
                                                if row['Protein'].split('|')[1] in acc2tax_dict else 'DECOY',
                                                                           axis=1)
        elif db_type == 'ncbi':
            self.sorted_xtandem_df['Protein'] = self.sorted_xtandem_df.apply(lambda row: row['Protein'].split('|')[2].split(), axis=1)
            self.sorted_xtandem_df['taxID'] = self.sorted_xtandem_df.apply(lambda row:
                                                                           acc2tax_dict[row['Protein']]
                                                                           if row['Protein'] in acc2tax_dict else 'DECOY',
                                                                           axis=1)
        elif db_type == 'custom':
            self.sorted_xtandem_df['taxID'] = self.sorted_xtandem_df.apply(lambda row:
                                                acc2tax_dict[row['Protein'].strip()]
                                                if not decoy_tag in row['Protein'] else 'DECOY',
                                                                           axis=1)
        print('df with taxID column')
        self.sorted_xtandem_df['taxID_level'] = self.sorted_xtandem_df.apply(lambda row:
                                                    taxon_graph.find_level_up(int(row['taxID']), level)
                                                    if row['taxID'] != 'DECOY' and row['taxID'] != 0 else 'DECOY',
                                                                             axis=1)
        print('df with taxID_level column')
        reduced_df = self.sorted_xtandem_df.groupby(["#SpecFile", 'Title', 'Peptide', 'Hyperscore'], as_index=False).agg(
            {'Protein': lambda x: set(x), 'EValue': lambda x: set(list(x)), 'decoy': lambda x: set(x),
            'taxID': lambda x: set(x), 'taxID_level': lambda x: set(x)})
        print('reduced df')
        reduced_df = reduced_df.sort_values(by=['Hyperscore', 'Title'], ascending=False).reset_index(drop=True)
        print(f"writing data frame to {self.path_to_file+'.csv'}... ")
        reduced_df.to_csv(self.path_to_file+'.csv', sep='\t')
        print(f'entries in df: {len(self.sorted_xtandem_df)} '
              f'decoys in df: {len(self.sorted_xtandem_df[self.sorted_xtandem_df.decoy==True])} '
              f'hits in df: {len(self.sorted_xtandem_df[self.sorted_xtandem_df.decoy==False])}')


    def determine_FDR_position(self, fdr):
        """
        :param decoy_tag: tag of decoy entries, for example REVERSED
        :param fdr: false discovery rate, for example 0.01
        """
        repeatedly_identified_spectra = set()
        number_multiple_identified_spectra = 0
        hits = decoy = 0
        FDR_position_not_set = True
        title_set = set()
        spectra_header = [(x, y) for x, y in zip(self.sorted_xtandem_df['Title'], self.sorted_xtandem_df['decoy'])]
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
                self.fdr_pos = hits + decoy - 1 + number_multiple_identified_spectra
                self.number_psms = hits
                self.decoys = decoy - 1
                FDR_position_not_set = False
                break
            if decoy / (hits + decoy) <= fdr:
                continue
                #and not FDR_position_not_set:
                # self.fdr_pos = hits + decoy + number_multiple_identified_spectra
                # self.number_psms = hits
                # self.decoys = decoy - 1
                # FDR_position_not_set = True
        print('Number of PSMs: %d' % self.number_psms)
        print('Number of decoys: %d' % self.decoys)
        print(f"double identified spetra {number_multiple_identified_spectra}")
        print('Position FDR border/Number of PSMs: %d' % self.fdr_pos)

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


