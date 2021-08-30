from collections import defaultdict
from pathlib import Path
import pandas as pd
from PSM_FDR import PSM_FDR
from create_reference_from_tsv_and_pepxml import ReferenceWriter
from handling_acc_files import HelperMethod


class PsmNumberPerTaxIDs():

    def __init__(self, taxa_set, path_to_reduced_df, level, fdr=0.1, add_taxa=None, custom_taxa=None):
        self.kleiner_taxIDs = [536, 882, 44577, 228410, 323848, 46170, 93061, 224308, 99287, 511145, 176299, 216596, 1041145,
                               262724, 266264, 266265, 1004788, 1114970, 1149133, 1294143, 1407502, 3055, 1302247, 926571,
                               # virus
                               10754, 101570, 1985310, 329852, 1283336, 12022, 1977402]
        self.kleiner_taxIDs_virus = [10754, 101570, 1985310, 329852, 1283336, 12022, 1977402]
        self.tanca_taxIDs = [747, 5535, 655183, 1579, 1255, 4932, 1465, 1351, 562]
        self.taxa_set = taxa_set
        if taxa_set == 'kleiner':
            self.taxIDs = self.kleiner_taxIDs
        elif taxa_set == 'tanca':
            self.taxIDs = self.tanca_taxIDs
        if add_taxa:
            self.taxIDs.extend(add_taxa)
        if custom_taxa:
            self.taxIDs = custom_taxa
        self.level = level
        self.reduced_df_in_fdr = self.get_df_in_fdr(path_to_reduced_df, fdr)
        self.nb_all_identified_spectra_for_uniprot = len(set(self.reduced_df_in_fdr["Title"]))
        self.taxon_graph = HelperMethod.load_taxa_graph(Path("/home/jules/Documents/databases/databases_tax2proteome/taxdump.tar.gz"))

    def get_df_in_fdr(self, path_to_reduced_df, fdr):
        reduced_df = ReferenceWriter.read_csv_with_generic_function(path_to_reduced_df,['Protein', 'Hyperscore', 'decoy', 'taxID'])
        fdr_pos_result, number_psm_result, number_decoy_result, double_spectra_result, score_last_item_result =PSM_FDR.determine_FDR_position(reduced_df, fdr, True)
        return number_psm_result, reduced_df[0:fdr_pos_result]

    def count_spectra_per_taxon(self, reduced_df_in_fdr, level, taxon=None, taxa_list=None):
        spectra = set()
        for spectrum_ID, taxa_set in zip(list(reduced_df_in_fdr['Title']), list(reduced_df_in_fdr[f'taxID_{level}'])):
            if taxon:
                if taxon in taxa_set:
                    spectra.add(spectrum_ID)
            elif taxa_list:
                for taxon in taxa_list:
                    if taxon in taxa_set:
                        spectra.add(spectrum_ID)
        return (spectra)

    def count_all_taxa(self, level):
        result_spectra_dict = {}
        for taxon in self.taxIDs:
            taxon = self.taxon_graph.find_level_up(taxon, level)
            result_spectra_dict[str(taxon)]=self.count_spectra_per_taxon(self.reduced_df_in_fdr, taxon)
        if self.taxa_set == 'kleiner':
            result_spectra_dict['virus']=set()
            for virus_taxon in self.kleiner_taxIDs_virus:
                result_spectra_dict['virus'].union(result_spectra_dict[str(virus_taxon)])
        return result_spectra_dict

    def flatten_list(self, l):
        return [item for sublist in l for item in sublist]

    def get_percentage_psm_per_taxon(self, level, taxon_str):
        all_counts = self.count_all_taxa(level)
        final_spectra_set= set()
        for taxon in taxon_str.split(', '):
            final_spectra_set=final_spectra_set.union(all_counts[taxon])
        psm_percentage = len(final_spectra_set)/self.nb_all_identified_spectra_for_uniprot*100
        return psm_percentage

    def count_row_by_row(self):
        taxID_to_spectra_dict = defaultdict(set)
        for spectrum, taxID_set in zip(self.reduced_df_in_fdr['Title'], self.reduced_df_in_fdr[f'taxID_{self.level}']):
            for taxID in taxID_set:
                taxID_to_spectra_dict[taxID].add(spectrum)
        print(taxID_to_spectra_dict.keys())

