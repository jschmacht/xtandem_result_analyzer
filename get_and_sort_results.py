
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

rows_sensitivity_analysis = ["specificity", "sensitivity", "TP", "FP", "TN", "FN", "Number of spectra identified in result but not in reference",
           "Number of spectra identified in result", "Number of spectra identified in reference", "FDR Position reference",
           "Number of PSMs in result", "Number of Decoys in result", "Number of double identified spectra in result",
           "Hyperscore of last item in FDR boundaries"]
fdr = 0.05
# size crap = 116 *2 =232
dbs = ['reference',"reference_Arabidopsis_thaliana", 'ncbi_species', "ncbi_genus", "uniprot_species", "uniprot_species_species",
      "uniprot_species_nr", "uniprot_genus",  "uniprot_genus_species",  "uniprot_genus_nr",  "uniprot_family_nr", "swissprot_species",
      "swissprot_genus", "swissprot_family", "swissprot_order"]
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

dict_uncomplete = {"uniprot_family": (29773074, ''),
                   'ncbi_species_species': (0, f"{species_reduced_df[4]}_{fdr}_sensitivity"),
                   "ncbi_genus_species":(0, f"{genus_reduced_df[5]}_{fdr}_sensitivity"),
                   }


def read_sensitivity_result_files(input):
    with open(input, 'r') as input:
        # discard fdr line
        result = []
        input.readline()
        for line in input.readlines():
            result.append(line.split('\t')[1])

    result_dict = {}
    for i, item in enumerate(result):
        result_dict[rows_sensitivity_analysis[i]] = item

    return result_dict

def write_all_results_into_one_file(output_file, dbs, results, dict_databases_to_size_and_result_file=None):
    with open(output_file, 'w') as out:
        out.write('db_name\t')
        for db in dbs:
            out.write(db + '\t')
        out.write('\n')

        for row_name in rows_sensitivity_analysis:
            out.write(row_name + '\t')
            for db in dbs:
                out.write(results[db][row_name].strip()+ '\t')
            out.write('\n')

        if dict_databases_to_size_and_result_file:
            out.write('db_size\t')
            for db in dbs:
                out.write(str(dict_databases_to_size_and_result_file[db][0]) + '\t')
            out.write('\n')


output_file = "/home/jules/Documents/Tax2Proteome/benchmarking/sensitivity_analysis.tsv"
results = {}
for db in dbs:
    results[db] = read_sensitivity_result_files(dict_databases_to_size_and_result_file[db][1])

write_all_results_into_one_file(output_file, dbs, results, dict_databases_to_size_and_result_file)

file = 'uniprot_species_nr'
output_file_fdr = f"{species_reduced_df[2]}_all_fdr_sensitivity.tsv"
fdrs = [0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
result_fdr = {}
filename_fdr = "uniprot_species_nr"
for fdr in fdrs:
    result_fdr[f"{filename_fdr}_{fdr}"] = read_sensitivity_result_files(f"{species_reduced_df[2]}_{fdr}_sensitivity")
dbs_FDR = [f"{filename_fdr}_{fdr}" for fdr in fdrs]
write_all_results_into_one_file(output_file_fdr, dbs_FDR, result_fdr)
