import pandas as pd
from pathlib import Path
import argparse
from TaxonGraph import TaxonGraph
import pickle
from PSM_FDR import PSM_FDR
from handling_acc_files import HelperMethod
from create_reference_from_tsv_and_pepxml import ReferenceWriter
from AccessionSearcher import AccessionSearcherNCBI
from ReadAccTaxon import ReadAccTaxon
from SearchAccessions import SearchAccessions
from collections import defaultdict


class DeterminatorSpecificitySensitivity():

    def __init__(self, level, fdr_applied_df, reference_df, spectra_file):
        """
        :param fdr_applied_df:
        :param refernce_df: column names =
        :param spectra_file: 'Run1_U1_2000ng.mgf'
        """
        self.tax_level = ['species', 'genus', 'family', 'order']
        self.result_df = fdr_applied_df[["Title", "Peptide", "Hyperscore", 'Protein', "decoy", f"taxID"]]
        self.result_df.columns = ['SpectraID', 'Peptide', 'Hyperscore', 'ProteinAcc', 'decoy', 'taxID_DB']
        self.reference_df = reference_df[['SpectraID', 'Ref_Peptide', 'Ref_Hyperscore', 'Ref_ProteinAcc', 'Ref_decoy', 'Ref_taxID_DB', f"Ref_taxID_{level}"]]
        self.all_spectra_list = self.get_all_spectra_IDs(spectra_file)

    def create_df_with_all_spectra_reference_and_result_taxa(self, level, taxon_graph, path_to_out):
        pd.set_option("display.max_rows", None, "display.max_columns", None)
        df_all_spectra = pd.DataFrame(self.all_spectra_list, columns=['SpectraID'])
        self.result_df[f'taxID_{level}'] = self.result_df['taxID_DB'].apply(lambda taxid: taxon_graph.find_level_up(taxid, level))
        print(df_all_spectra.head())
        df_with_all_spectra_and_reference_and_results = pd.merge(df_all_spectra, self.result_df, how="outer", on=['SpectraID'])
        print(df_with_all_spectra_and_reference_and_results.columns)
        print(df_with_all_spectra_and_reference_and_results.head())
        df_with_all_spectra_and_reference_and_results = pd.merge(df_with_all_spectra_and_reference_and_results,
                                                                 self.reference_df, how="outer", on=['SpectraID'])
        print(df_with_all_spectra_and_reference_and_results.columns)
        print(df_with_all_spectra_and_reference_and_results.head())
        print(f"df_with_all_spectra_and_reference_and_results {path_to_out}... ")
        df_with_all_spectra_and_reference_and_results.to_csv(str(path_to_out), sep='\t')


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

    def get_true_positive_and_true_negative(self, level, fdr_applied_df, refernce_df, spectra_file):
        all_spectra_IDs = self.get_all_spectra_IDs(spectra_file)

        print(fdr_applied_df.columns)
        fdr_applied_df['taxID_level'] = fdr_applied_df[f'taxID_{level}'].apply(lambda taxid: taxid[1:-1].split(', '))
        fdr_applied_df['Protein'] = fdr_applied_df['Protein'].apply(lambda taxid: taxid[1:-1].split(', '))
        spectraID_to_taxid_dict = self.load_ref_file(refernce_df, level)
        fdr_applied_df['TP'] = fdr_applied_df.apply(lambda row: 'FP' if row['Title'] not in spectraID_to_taxid_dict.keys() else (
            'TP' if (set(spectraID_to_taxid_dict[row['Title']]) & set(fdr_applied_df['taxID_level'])) else ''), axis=1)
        number_TN = 0
        number_FN = 0
        for spectra in all_spectra_IDs:
            if spectra not in fdr_applied_df['Title'] and spectra not in spectraID_to_taxid_dict.keys():
                number_TN += 1
            if spectra not in fdr_applied_df['Title'] and spectra in spectraID_to_taxid_dict.keys():
                number_FN += 1
        number_TP = fdr_applied_df[fdr_applied_df['TP'] == 'TP'].count()
        number_FP = fdr_applied_df[fdr_applied_df['TP'] == 'FP'].count()
        print(f"TP: {number_TP}, FP: {number_FP}, TN: {number_TN}, FN: {number_FN}")


def get_decoy_set(decoy_str):
    if 'False' in decoy_str and 'True' in decoy_str:
        return {True, False}
    elif 'True' in decoy_str:
        return {True}
    elif 'False' in decoy_str:
        return {False}


def main():
    parser = argparse.ArgumentParser(description='Read xtandem output .tsv')
    parser.add_argument('-i', '--input', dest='input', default=None, help='xtandem tsv.tsv with columns:'
                                                                          '[spectra_file, Spectrum ID, Peptide, Protein Accession + Description, Hyperscore, Evalue,'
                                                                          'or pep.xml file for create_ref')
    parser.add_argument('-r', '--reference', dest='reference', default=None, help='Path to reference tsv')
    parser.add_argument('-g', '--tax_graph', dest='tax_graph', help='Path to taxdump.tar.gz')
    parser.add_argument('-l', '--level', dest='level', choices=['subspecies', 'species', 'genus', 'family', 'order', 'superkingdom'],
                        help='Level of database')
    parser.add_argument('-d', '--database', dest='database', choices=['ncbi', 'uniprot', 'custom'], default='uniprot',
                        help='Database format.')
    parser.add_argument('-s', '--spectra_file', dest='spectra_file', default=None,
                        help='path to Run1_U1_2000ng.mgf spectra file')
    parser.add_argument('-c', '--create_reference', dest='create_reference', action='store_true', default=False,
                        help='create reference from pep.xml')
    parser.add_argument('-f', '--fdr', dest='fdr', type=float, default=0.01, help='FDR-rate, default  = 0.01')
    parser.add_argument('-y', '--decoy', dest='decoy', default='REVERSED', help='Decoy_tag.')
    options = parser.parse_args()

    taxon_graph = HelperMethod.load_taxa_graph(Path(options.tax_graph))
    path_to_result = Path(options.input)
    path_to_reference = Path(options.reference)
    if options.create_reference and path_to_result.suffixes == ['.pep', '.xml']:
        custom_acc2tax_file_based_on_Kleiner_DB = '/home/jules/Documents/Tax2Proteome/benchmarking/Kleiner_ref_db/acc2tax_custom'
        reader = ReferenceWriter(path_to_result, custom_acc2tax_file_based_on_Kleiner_DB, path_to_reference, taxon_graph)
        spectra_to_accs_dict = reader.read_pepXML()
        reader.write_Kleiner_spectrum_reference_file(spectra_to_accs_dict)
    elif options.create_reference and path_to_result.suffix == '.tsv':
        reduced_df = pd.read_csv(str(path_to_result), sep='\t')
        reduced_df['decoy']= reduced_df['decoy'].apply(lambda decoy_str: get_decoy_set(decoy_str) )
        psm = PSM_FDR('')
        fdr_pos, number_psms, decoys = psm.determine_FDR_position(reduced_df, options.fdr, True)
        ReferenceWriter.write_result_spectrum_reference_file(path_to_result, ['species', 'genus', 'family', 'order'], taxon_graph, path_to_reference, fdr_pos)

    else:
        path_to_all_info_tsv = path_to_result.parent.joinpath(path_to_result.stem + '_' + path_to_reference.stem + '.tsv')
        print(path_to_all_info_tsv)
        result_df = pd.read_csv(str(options.input), sep='\t')
        reference_df = pd.read_csv(str(options.reference), sep='\t')
        psm = PSM_FDR(options.input)
        fdr_pos, number_psms, decoys = psm.determine_FDR_position(result_df, options.fdr)
        fdr_applied_df = result_df[0:fdr_pos]

        determinator = DeterminatorSpecificitySensitivity(options.level, fdr_applied_df, reference_df, options.spectra_file)
        determinator.create_df_with_all_spectra_reference_and_result_taxa(options.level, taxon_graph, path_to_all_info_tsv)
        #determinator.get_true_positive_and_true_negative(options.level, fdr_applied_df, reference_df, options.spectra_file)

if __name__ == '__main__':
    main()
