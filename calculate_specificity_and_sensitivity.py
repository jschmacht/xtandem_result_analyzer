import pandas as pd
from pathlib import Path
import argparse
from create_reference_from_tsv_and_pepxml import ReferenceWriter
from sensitivity_calculator import SensitivityAndSpecificity


def write_output(path_to_result_file, fdr, specificity, sensitivity, TP, FP, TN, FN, s_u_s, df_only_in_result_identified_spectra, ignore_unidentified=False):
    if not ignore_unidentified:
        output_file=f"{str(path_to_result_file)}_{fdr}_sensitivity_new"
    else:
        output_file=f"{str(path_to_result_file)}_ignore_unclassified_{fdr}_sensitivity_new"
    with open(output_file, 'w') as output:
        output.write('FDR:' + '\t' + str(fdr) + '\n')
        output.write('specificity:'+ '\t'  + str(specificity) + '\n')
        output.write('sensitivity:'+ '\t'  + str(sensitivity) + '\n')
        output.write('TP:' + '\t'  + str(TP) + '\n')
        output.write('FP:' + '\t'  + str(FP) + '\n')
        output.write('TN:' + '\t'  + str(TN) + '\n')
        output.write('FN:' + '\t'  + str(FN) + '\n')
        output.write('Number of spectra identified in result but not in reference:' + '\t'  + str(len(set(df_only_in_result_identified_spectra.Title.tolist()))) + '\n')
        output.write('Number of spectra identified in result:' + '\t' + str(len(set(s_u_s.result_df[0:s_u_s.fdr_pos_result].Title.tolist()))) + '\n')
        output.write('Number of spectra identified in reference:' + '\t' + str(len(set(s_u_s.reference_df[0:s_u_s.fdr_pos_reference].Title.tolist())))+ '\n')
        output.write('FDR Position reference:' + '\t'  + str(s_u_s.fdr_pos_result) + '\n')
        output.write('Number of PSMs in result:' + '\t'  + str(s_u_s.number_psm_result)+ '\n')
        output.write('Number of Decoys in result:' + '\t'  + str(s_u_s.number_decoy_result)+ '\n')
        output.write('Number of double identified spectra in result:' + '\t'  + str(s_u_s.double_spectra_result) + '\n')
        output.write('Hyperscore of last item in FDR boundaries:'+ '\t' + str(s_u_s.score_last_item_result)  + '\n')

def main():
    parser = argparse.ArgumentParser(description='Read reduced .tsv')
    parser.add_argument('-r', '--result', dest='result', default=None, help='Path to reduced result tsv, Tax2Proteome X-Tandem results')
    parser.add_argument('-l', '--level', dest='level', choices=['species', 'genus', 'family', 'order'],
                        help='Level of database')
    parser.add_argument('-f', '--fdr', dest='fdr', type=float, default=0.05, help='FDR-rate, default = 0.05')
    parser.add_argument('-i', '--ignore_unidentified_spectra', action='store_true', default=False,
                        help='all spectra unidentified in reference and result file are ignored for '
                                            'specificity and sensitivity analysis, default = False')

    options = parser.parse_args()

    path_to_result_file = Path(options.result)
    path_to_reference = Path('/home/jules/Documents/Tax2Proteome/benchmarking/results_searchgui_kleiner_db/Run1_U1_2000ng.t.xml_reduced_reference.tsv')
    path_to_spectra = Path('/home/jules/Documents/Tax2Proteome/benchmarking/spectra/Run1_U1_2000ng.mgf')
    result_df = ReferenceWriter.read_csv_with_generic_function(path_to_result_file,
                                                               ['Protein', 'Hyperscore', 'decoy', 'taxID', f'taxID_{options.level}'])
    reference_df = ReferenceWriter.read_csv_with_generic_function(path_to_reference,
                                                                  ['Protein', 'Hyperscore', 'decoy', 'taxID', f'taxID_{options.level}'])
    print('Dataframes loaded.')
    s_u_s = SensitivityAndSpecificity(reference_df, result_df, path_to_spectra, options.level, options.fdr)
    print('Initialized.')
    # only Spectra identified in Reference
    df_with_all_reference_spectra_and_merged_results_in_fdr = pd.merge(s_u_s.reference_df[0:s_u_s.fdr_pos_reference],
                                                                s_u_s.result_df[0:s_u_s.fdr_pos_result], how="left",
                                                                left_on='Title', right_on='Title')
    print('Merged dataframe')
    if not options.ignore_unidentified_spectra:
        df_with_all_unidentified_spectra =  s_u_s.get_df_with_all_unidentified_spectra_in_reference_and_result()
        TP, FP, TN, FN = s_u_s.get_true_positive_and_true_negative(df_with_all_reference_spectra_and_merged_results_in_fdr,
                                                                   df_with_all_unidentified_spectra)
    else:
        TP, FP, TN, FN = s_u_s.get_true_positive_and_true_negative(df_with_all_reference_spectra_and_merged_results_in_fdr)
    sensitivity = SensitivityAndSpecificity.calculate_sensitivity(TP, FN)
    specificity = SensitivityAndSpecificity.calculate_specificity(FP, TN)
    df_only_in_result_identified_spectra = s_u_s.get_df_with_identified_spectra_in_result_df_but_not_in_reference_df()

    write_output(path_to_result_file, options.fdr, specificity, sensitivity, TP, FP, TN, FN, s_u_s, df_only_in_result_identified_spectra, options.ignore_unidentified_spectra)



if __name__ == '__main__':
    main()