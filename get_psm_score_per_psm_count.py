from create_reference_from_tsv_and_pepxml import ReferenceWriter
from create_PSM_df import PSM_FDR

# psm score vs psm count
def get_psm_score_count(level_list, file_list, in_fdr=True, fdr=0.05):
    result = {}
    for level, file in zip(level_list, file_list):
        result_df = ReferenceWriter.read_csv_with_generic_function(file, ['Hyperscore', 'decoy'])
        result_df = result_df[['Title',  'decoy']]
        if in_fdr:
            fdr_pos, number_psm, number_decoy, double_spectra, score_last_item = PSM_FDR.determine_FDR_position(result_df, fdr, True)
            result_df = result_df[0:fdr_pos]
        reversed_df = result_df.sort_values(by=['Hyperscore', 'Title'], ascending=True).reset_index(drop=True)
        # print(reversed_df.head(10))
        print('Dataframe loaded.')
        #psm count = len set(spectra)
        spectra_set=set()
        decoy_set=set()
        psm_score = []
        psm_count = []
        decoy_count = []
        decoy_score = []
        last_psm_score = 0
        last_decoy_score = 0
        for (spectra, score, decoy) in zip(reversed_df.Title.tolist(), reversed_df.Hyperscore.tolist(), reversed_df.decoy.tolist()):
            if decoy == {True}:
                decoy_set.add(spectra)
                if score == last_decoy_score:
                    decoy_count[-1]=len(decoy_set)
                else:
                    decoy_score.append(score)
                    decoy_count.append(len(decoy_set))
                last_decoy_score = score
            else:
                spectra_set.add(spectra)
                if score == last_psm_score:
                    psm_count[-1]=len(spectra_set)
                else:
                    psm_score.append(score)
                    psm_count.append(len(spectra_set))
                last_psm_score = score
        result[level] = {'psm_count': psm_count, 'psm_score': psm_score, 'decoy_count': decoy_count, 'decoy_score':decoy_score}
    return result


def main():

    species_reduced_df= ["/home/jules/Documents/Tax2Proteome/benchmarking/results_searchgui_xtandem_analyzer_bachelor_thesis/uniprot/x_tandem_tsv/Run1_U1_2000ng_uniprot_species.t.xml_reduced.tsv",
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

    family_reduced_df=["/home/jules/Documents/Tax2Proteome/benchmarking/results_searchgui_xtandem_analyzer_bachelor_thesis/uniprot/x_tandem_tsv/Run1_U1_2000ng_uniprot_family_nr.t.xml_reduced.tsv"
    "/home/jules/Documents/Tax2Proteome/benchmarking/results_searchgui_xtandem_analyzer_bachelor_thesis/swissprot/x_tandem_tsv/Run1_U1_2000ng_swissprot_family.t.xml_reduced.tsv"]


    order_reduced_df=["/home/jules/Documents/Tax2Proteome/benchmarking/results_searchgui_xtandem_analyzer_bachelor_thesis/swissprot/x_tandem_tsv/Run1_U1_2000ng_swissprot_order.t.xml_reduced.tsv"]

    uniprot_list = ["/home/jules/Documents/Tax2Proteome/benchmarking/results_searchgui_xtandem_analyzer_bachelor_thesis/uniprot/x_tandem_tsv/Run1_U1_2000ng_uniprot_species_nr.t.xml_reduced.tsv",
                    "/home/jules/Documents/Tax2Proteome/benchmarking/results_searchgui_xtandem_analyzer_bachelor_thesis/uniprot/x_tandem_tsv/Run1_U1_2000ng_uniprot_genus_nr.t.xml_reduced.tsv",
                    "/home/jules/Documents/Tax2Proteome/benchmarking/results_searchgui_xtandem_analyzer_bachelor_thesis/uniprot/x_tandem_tsv/Run1_U1_2000ng_uniprot_family_nr.t.xml_reduced.tsv"]
    level_list = ['species', 'genus', 'family']
    result_uniprot_all = get_psm_score_count(level_list, uniprot_list, in_fdr=False)
    result_uniprot_in_fdr = get_psm_score_count(level_list, uniprot_list, in_fdr=True)


    level_list = ['species', 'genus']
    ncbi_list = ["/home/jules/Documents/Tax2Proteome/benchmarking/results_searchgui_xtandem_analyzer_bachelor_thesis/ncbi_kleiner/x_tandem_tsv/Run1_U1_2000ng_ncbi_kleiner_species.t.xml_reduced.tsv",
                 "/home/jules/Documents/Tax2Proteome/benchmarking/results_searchgui_xtandem_analyzer_bachelor_thesis/ncbi_kleiner/x_tandem_tsv/Run1_U1_2000ng_ncbi_kleiner_genus.t.xml_reduced.tsv"]
    result_ncbi_all = get_psm_score_count(level_list, ncbi_list, in_fdr=False)
    result_ncbi_in_fdr = get_psm_score_count(level_list, ncbi_list)



if __name__ == '__main__':
    main()