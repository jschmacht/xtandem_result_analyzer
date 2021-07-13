import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path
import argparse
from create_reference_from_tsv_and_pepxml import ReferenceWriter
from sensitivity_calculator import SensitivityAndSpecificity
from create_PSM_df import PSM_FDR

def read_all_results_into_dict(file):
    result_dict= {}
    with open(file, 'r') as inp:
        firstline=inp.readline()
        db_names = [name.strip() for name in firstline.strip().split('\t')[1:]]
        for db in db_names:
            result_dict[db] = {}
        for line in inp.readlines():
            result_name = line.split('\t')[0]
            result_fields = [name.strip() for name in line.split('\t')[1:]]
            for i, db in enumerate(db_names):
                result_dict[db][result_name]=result_fields[i]
    return result_dict

def read_analysis_results_into_dict(file):
    result_dict= {}
    with open(file, 'r') as inp:
        for line in inp.readlines():
            result_dict[line.split('\t')[0]] = line.split('\t')[1].strip()
    #print(result_dict)
    return result_dict

def read_TP(file):
    TP = int(read_analysis_results_into_dict(file)['TP:'])
    return TP

def read_TN(file):
    TN = int(read_analysis_results_into_dict(file)['TN:'])
    return TN

def read_FP(file):
    TP = int(read_analysis_results_into_dict(file)['FP:'])
    return TP

def read_FN(file):
    TN = int(read_analysis_results_into_dict(file)['FN:'])
    return TN

def read_sensitivity(file):
    sen = float(read_analysis_results_into_dict(file)['sensitivity:'])
    return sen

def read_specificity(file):
    spe = float(read_analysis_results_into_dict(file)['specificity:'])
    return spe

def read_last_score(file):
    last_score = read_analysis_results_into_dict(file)['Hyperscore of last item in FDR boundaries:']
    return last_score

def get_score_results(filtered_dict):
    db_size = []
    score_border = []
    result_dict = {}
    for file_name, db_size_file_tuple in filtered_dict.items():
        result_dict[file_name] = {'db_size': db_size_file_tuple[0]/100000, 'score_border': read_last_score(db_size_file_tuple[1])}

    for k, v in result_dict.items():
        db_size.append(v['db_size'])
        score_border.append(v['score_border'])
    return result_dict, db_size, score_border

def get_specificity_results(filtered_dict):
    db_size = []
    score_border = []
    result_dict = {}
    for file_name, db_size_file_tuple in filtered_dict.items():
        result_dict[file_name] = {'db_size': db_size_file_tuple[0]/100000, 'specificity': read_specificity(db_size_file_tuple[1])}
    for k, v in result_dict.items():
        db_size.append(v['db_size'])
        score_border.append(v['specificity'])
    return result_dict, db_size, score_border

def get_decoy_rows(decoy_column):
    return [True if d=={True} else False for d in decoy_column]

def get_tax_rows(tax_column, taxid):
    return [True if tax_set == {taxid} else False for tax_set in tax_column]

def get_decoy_rows(protein_column):
    return [True if 'REVERSED' in p else False for p in protein_column]


def get_hit_rows(protein_column):
    return [False if 'REVERSED' in p else True for p in protein_column]

def get_decoy_rows2(decoy_column):
    return [True if d_set in [{True, False}, {True}] else False for d_set in decoy_column]

def get_hit_rows2(decoy_column):
    return [True if d_set in  [{True, False}, {False}] else False for d_set in decoy_column]

def create_df_for_lineplot_from_rows_df2(column_names, rows):
    df_spe_sens = pd.DataFrame(rows,columns=column_names)
    return df_spe_sens

def plot_lineplot_sens_spec2(df_spe_sens_wide):
    plt.figure(figsize=(10,6))
    sns.set_theme(style="whitegrid")
    sns.lineplot(x='level', y='value', data=df_spe_sens_wide, hue='database', style='sens/spe')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
    plt.ylabel("specificity/sensitivity")
    plt.savefig('/home/jules/Documents/Tax2Proteome/benchmarking/plots/specificity_sensitivity.png')
    return plt

def main():
    # Files

    species_reduced_df=["/home/jules/Documents/Tax2Proteome/benchmarking/results_searchgui_xtandem_analyzer_bachelor_thesis/uniprot/x_tandem_tsv/Run1_U1_2000ng_uniprot_species.t.xml_reduced.tsv",
                        "/home/jules/Documents/Tax2Proteome/benchmarking/results_searchgui_xtandem_analyzer_bachelor_thesis/uniprot/x_tandem_tsv/Run1_U1_2000ng_uniprot_species_species.t.xml_reduced.tsv",
                        "/home/jules/Documents/Tax2Proteome/benchmarking/results_searchgui_xtandem_analyzer_bachelor_thesis/uniprot/x_tandem_tsv/Run1_U1_2000ng_uniprot_species_nr.t.xml_reduced.tsv",
                        "/home/jules/Documents/Tax2Proteome/benchmarking/results_searchgui_xtandem_analyzer_bachelor_thesis/swissprot/x_tandem_tsv/Run1_U1_2000ng_swissprot_species.t.xml_reduced.tsv",
                        "/home/jules/Documents/Tax2Proteome/benchmarking/results_searchgui_xtandem_analyzer_bachelor_thesis/ncbi_kleiner/x_tandem_tsv/Run1_U1_2000ng_ncbi_kleiner_species_species.t.xml_reduced.tsv",
                        "/home/jules/Documents/Tax2Proteome/benchmarking/results_searchgui_xtandem_analyzer_bachelor_thesis/ncbi_kleiner/x_tandem_tsv/Run1_U1_2000ng_ncbi_kleiner_species.t.xml_reduced.tsv",
                        "/home/jules/Documents/Tax2Proteome/benchmarking/results_searchgui_kleiner_reference_aradiopsis/Run1_U1_2000ng_kleiner_aradiopsis.t.xml_reduced.tsv",
                        "/home/jules/Documents/Tax2Proteome/benchmarking/results_searchgui_kleiner_db/Run1_U1_2000ng.t.xml_reduced.tsv"]

    genus_reduced_df=["/home/jules/Documents/Tax2Proteome/benchmarking/results_searchgui_xtandem_analyzer_bachelor_thesis/uniprot/x_tandem_tsv/Run1_U1_2000ng_uniprot_genus.t.xml_reduced.tsv",
                      "/home/jules/Documents/Tax2Proteome/benchmarking/results_searchgui_xtandem_analyzer_bachelor_thesis/uniprot/x_tandem_tsv/Run1_U1_2000ng_uniprot_genus_species.t.xml_reduced.tsv",
                      "/home/jules/Documents/Tax2Proteome/benchmarking/results_searchgui_xtandem_analyzer_bachelor_thesis/uniprot/x_tandem_tsv/Run1_U1_2000ng_uniprot_genus_nr.t.xml_reduced.tsv",
                      "/home/jules/Documents/Tax2Proteome/benchmarking/results_searchgui_xtandem_analyzer_bachelor_thesis/swissprot/x_tandem_tsv/Run1_U1_2000ng_swissprot_genus.t.xml_reduced.tsv",
                      "/home/jules/Documents/Tax2Proteome/benchmarking/results_searchgui_xtandem_analyzer_bachelor_thesis/ncbi_kleiner/x_tandem_tsv/Run1_U1_2000ng_ncbi_kleiner_genus.t.xml_reduced.tsv",
                      "/home/jules/Documents/Tax2Proteome/benchmarking/results_searchgui_xtandem_analyzer_bachelor_thesis/ncbi_kleiner/x_tandem_tsv/Run1_U1_2000ng_ncbi_kleiner_genus_species.t.xml_reduced.tsv"]

    family_reduced_df=["/home/jules/Documents/Tax2Proteome/benchmarking/results_searchgui_xtandem_analyzer_bachelor_thesis/uniprot/x_tandem_tsv/Run1_U1_2000ng_uniprot_family_nr.t.xml_reduced.tsv",
                       "/home/jules/Documents/Tax2Proteome/benchmarking/results_searchgui_xtandem_analyzer_bachelor_thesis/swissprot/x_tandem_tsv/Run1_U1_2000ng_swissprot_family.t.xml_reduced.tsv"]

    order_reduced_df=["/home/jules/Documents/Tax2Proteome/benchmarking/results_searchgui_xtandem_analyzer_bachelor_thesis/swissprot/x_tandem_tsv/Run1_U1_2000ng_swissprot_order.t.xml_reduced.tsv"]
    fdr = 0.05
    dict_databases_to_size_and_result_file={'reference': (123088, f"{species_reduced_df[7]}_{fdr}_sensitivity"),
                                            "reference_Arabidopsis_thaliana": (138980, f"{species_reduced_df[6]}_{fdr}_sensitivity"),
                                            'ncbi_species': (8702135, f"{species_reduced_df[5]}_{fdr}_sensitivity"),
                                            "ncbi_genus":(27804893, f"{genus_reduced_df[4]}_{fdr}_sensitivity"),

                                            "uniprot_species": (4683371, f"{species_reduced_df[0]}_{fdr}_sensitivity"),
                                            "uniprot_species_species": (2093157, f"{species_reduced_df[1]}_{fdr}_sensitivity"),
                                            "uniprot_species_nr": (2991727, f"{species_reduced_df[2]}_{fdr}_sensitivity"),
                                            "uniprot_genus": (18352148,  f"{genus_reduced_df[0]}_{fdr}_sensitivity"),
                                            "uniprot_genus_species": (13068285,  f"{genus_reduced_df[1]}_{fdr}_sensitivity"),
                                            "uniprot_genus_nr": (13210287,  f"{genus_reduced_df[2]}_{fdr}_sensitivity"),

                                            "uniprot_family_nr": (22509624, f"{family_reduced_df[0]}_{fdr}_sensitivity"),
                                            "swissprot_species": (58505, f"{species_reduced_df[3]}_{fdr}_sensitivity"),
                                            "swissprot_genus": (88164, f"{genus_reduced_df[3]}_{fdr}_sensitivity" ),
                                            "swissprot_family": (124044, f"{family_reduced_df[1]}_{fdr}_sensitivity"),
                                            "swissprot_order": (181725, f"{order_reduced_df[0]}_{fdr}_sensitivity")}

    uniprot_dict =dict(filter(lambda item: 'uniprot' in item[0], dict_databases_to_size_and_result_file.items()))
    uniprot_nr_dict =dict(filter(lambda item: '_nr' in item[0], uniprot_dict.items()))
    uniprot_species_dict=dict(filter(lambda item: item[0] in ['uniprot_species_species', 'uniprot_genus_species'], uniprot_dict.items()))
    ncbi_dict =dict(filter(lambda item: 'ncbi' in item[0], dict_databases_to_size_and_result_file.items()))
    swissprot_dict=dict(filter(lambda item: 'swiss' in item[0], dict_databases_to_size_and_result_file.items()))

    column_names2 = ['database', 'level', 'value', 'sens/spe']
    rows2 =[['uniprot_nr', 'species', 98.323234, 'specificity'], ['uniprot_nr', 'genus', 98.432606, 'specificity'], ['uniprot_nr', 'family', 98.305966, 'specificity'],
            ['ncbi_nr', 'species', 98.562885, 'specificity'], ['ncbi_nr', 'genus', 98.320025, 'specificity'],
            ['swissprot', 'species', 97.838305, 'specificity'], ['swissprot', 'genus', 98.001373, 'specificity'], ['swissprot', 'family', 97.848361, 'specificity'],
            ['uniprot_nr', 'species',  85.891404, 'sensitivity' ], ['uniprot_nr', 'genus',  78.391253, 'sensitivity'], ['uniprot_nr', 'family',  74.921506, 'sensitivity' ],
            ['ncbi_nr', 'species',  85.726948 , 'sensitivity'], ['ncbi_nr', 'genus', 76.974026, 'sensitivity'],
            ['swissprot', 'species',  50.595859, 'sensitivity' ], ['swissprot', 'genus', 51.592057, 'sensitivity'], ['swissprot', 'family', 50.350655, 'sensitivity']]
    df_all=create_df_for_lineplot_from_rows_df2(column_names2, rows2)
    plt = plot_lineplot_sens_spec2(df_all)


if __name__ == '__main__':
    main()