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
        self.column_of_interest = f'taxID_{level}' if level != 'subspecies' else 'taxID_species'
        self.reduced_df_in_fdr, self.psm_count = self.get_df_in_fdr_and_psm_count(path_to_reduced_df, level, fdr)
        self.nb_all_identified_spectra_for_uniprot = len(set(self.reduced_df_in_fdr["Title"]))
        self.taxon_graph = HelperMethod.load_taxa_graph(Path("/home/jules/Documents/databases/databases_tax2proteome/taxdump.tar.gz"))

    def get_df_in_fdr_and_psm_count(self, path_to_reduced_df, level, fdr):
        reduced_df = ReferenceWriter.read_csv_with_generic_function(path_to_reduced_df,['Protein', 'Hyperscore', 'decoy', 'taxID', self.column_of_interest])
        fdr_pos, number_psms, decoys =PSM_FDR.determine_FDR_position(reduced_df, fdr, True)
        return reduced_df[0:fdr_pos], number_psms

    def count_spectra_per_taxon(self, reduced_df_in_fdr, level, taxon=None, taxa_list=None):
        spectra = set()
        for spectrum_ID, taxa_set in zip(list(reduced_df_in_fdr['Title']), list(reduced_df_in_fdr[self.column_of_interest])):
            if taxon:
                if taxon in taxa_set:
                    spectra.add(spectrum_ID)
            elif taxa_list:
                for taxon in taxa_list:
                    if taxon in taxa_set:
                        spectra.add(spectrum_ID)
        return spectra

    def count_all_taxa(self, level):
        result_spectra_dict = {}
        for taxon in self.taxIDs:
            taxon = self.taxon_graph.find_level_up(taxon, level)
            result_spectra_dict[str(taxon)]=self.count_spectra_per_taxon(self.reduced_df_in_fdr, taxon)
        if self.taxa_set == 'kleiner':
            result_spectra_dict['viruses']=set()
            for virus_taxon in self.kleiner_taxIDs_virus:
                result_spectra_dict['viruses'].union(result_spectra_dict[str(virus_taxon)])
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
        for spectrum, taxID_set in zip(self.reduced_df_in_fdr['Title'], self.reduced_df_in_fdr[self.column_of_interest]):
            for taxID in taxID_set:
                taxID_to_spectra_dict[taxID].add(spectrum)
        return taxID_to_spectra_dict

    def get_percentage(self, taxID_to_spectra_dict):
        taxID_to_percentage_dict = {}
        for taxID, spectra_set in taxID_to_spectra_dict.items():
            taxID_to_percentage_dict[taxID] = len(spectra_set)/self.psm_count*100
        return taxID_to_percentage_dict

    def get_virus_spectra(self, taxID_to_spectra_dict):
        for taxID in self.kleiner_taxIDs_virus:
            try:
                taxID_to_spectra_dict['viruses'] =  taxID_to_spectra_dict['viruses'].union(taxID_to_spectra_dict[taxID])
            except KeyError:
                print(taxID)
                continue
        return taxID_to_spectra_dict

def main():
    uniprot_species_reduced = "/home/jules/Documents/Tax2Proteome/benchmarking/results_reanalysis_uniprot/Run1_U1_2000ng_uniprot_species_nr.t.xml_reduced.tsv"
    obj = PsmNumberPerTaxIDs('kleiner', uniprot_species_reduced, 'species')

    taxID_to_spectra_dict = obj.count_row_by_row()
    print(taxID_to_spectra_dict)
    taxID_to_percentage_dict = obj.get_percentage(taxID_to_spectra_dict)
    print(taxID_to_percentage_dict)

if __name__ == '__main__':
    main()