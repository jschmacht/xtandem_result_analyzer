from pathlib import Path
import pandas as pd
import argparse
import pickle


class PSM_FDR:
    def __init__(self, path_to_file):
        """
        :param path_to_file: path to xtandem output .xml converted to tsv
        one line one header with all accessions seperated by \t
        """
        self.path_to_file = path_to_file
        self.sorted_xtandem_df = None
        self.spectra_acc_dict = {}
        self.spectra_sequence_dict = {}
        self.fdr_pos = 0
        self.psm = 0
        self.decoys = 0

    def create_PSM_dataframe(self):
        print('Reading xtandem tsv result as panda dataframe.')
        xtandem_df = pd.read_csv(str(self.path_to_file), delimiter='\t')
        print('Sorting panda dataframe by hyperscore.')
        self.sorted_xtandem_df = xtandem_df.sort_values(by=['Hyperscore', 'Title'], ascending=False)
        # print(sorted_xtandem_df['Hyperscore'][0:3])

    def determine_FDR_position(self, decoy_tag, fdr):
        """
        :param decoy_tag: tag of decoy entries, for example REVERSED
        :param fdr: false discovery rate, for example 0.01
        """
        doubles = pos = decoy = 0
        FDR_position_not_set = True
        title_set = set()
        spectra_header = [(x, y) for x, y in zip(self.sorted_xtandem_df['Title'], self.sorted_xtandem_df['Protein'])]
        for elem in spectra_header:
            if elem[0] in title_set:
                doubles += 1
                continue
            else:
                title_set.add(elem[0])
            if decoy_tag in (elem[1]):
                decoy += 1
            else:
                pos += 1
            if decoy / (pos + decoy) > fdr and FDR_position_not_set:
                self.fdr_pos = pos + decoy - 1 + doubles
                self.psm = pos
                self.decoys = decoy - 1
                FDR_position_not_set = False
            if decoy / (pos + decoy) <= fdr and not FDR_position_not_set:
                self.fdr_pos = pos + decoy + doubles
                self.psm = pos
                self.decoys = decoy - 1
                FDR_position_not_set = True
        print('Number of PSMs: %d' % self.psm)
        print('Number of decoys: %d' % self.decoys)
        print('Position FDR border/Number of PSMs: %d' % self.fdr_pos)

    # for every spectrum (key) all by xtandem found identifications = accessions (value)
    def create_spectrum_acc_dict(self, db):
        """
        :param db: db-format, 'uniprot' or 'ncbi'
        """
        spectra_header_psm = [(x, y) for x, y in zip(self.sorted_xtandem_df['Title'], self.sorted_xtandem_df['Protein'])][0:self.fdr_pos]

        for elem in spectra_header_psm:
            if elem[0] in self.spectra_acc_dict:
                if 'REVERSED' in elem[1]:
                    self.spectra_acc_dict[elem[0]].append('REVERSED')
                else:
                    if db == 'uniprot':
                        self.spectra_acc_dict[elem[0]].append(elem[1].split('|')[1])
                    else:
                        self.spectra_acc_dict[elem[0]].append('|'.join(elem[1].split('|')[2:]).split(' ')[0])
            else:
                if 'REVERSED' in elem[1]:
                    self.spectra_acc_dict[elem[0]] = ['REVERSED']
                else:
                    if db == 'uniprot':
                        self.spectra_acc_dict[elem[0]] = [elem[1].split('|')[1]]
                    else:
                        self.spectra_acc_dict[elem[0]] = ['|'.join(elem[1].split('|')[2:]).split(' ')[0]]

        print('Length spectra_acc_dict: %d ' % len(self.spectra_acc_dict))
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


