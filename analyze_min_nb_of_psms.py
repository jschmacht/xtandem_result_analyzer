import pandas as pd
import numpy as np
from pathlib import Path
from create_reference_from_tsv_and_pepxml import ReferenceWriter
from create_PSM_df import PSM_FDR
from collections import defaultdict
from ReadAccTaxon import ReadAccTaxon
from TaxonGraph import TaxonGraph
from DataframeMethods import remove_empty_rows, remove_acc_row, get_taxon_rows, get_taxa_rows, \
    get_all_exclusive_rows, get_hit_rows2


class defined_numer_of_psms_per_protein:
    def __init__(self, path_to_taxdump, level_to_tsv_path_dict, species_taxa_list, path_to_result,
                 level_to_taxa_list_dict=None, name_to_taxa_group_dict=None, fdr=0.05):
        self.level_to_tsv_path_dict = level_to_tsv_path_dict
        self.species_taxa_list = species_taxa_list
        self.name_to_taxa_group_dict = name_to_taxa_group_dict
        self.level_to_taxa_list_dict = level_to_taxa_list_dict
        self.path_to_result = path_to_result
        self.fdr = fdr
        self.taxon_graph = TaxonGraph()
        self.taxon_graph.create_graph(path_to_taxdump)

    # functions for loading xtandem results into df
    def get_psm_and_df_in_fdr(self, file, fdr, remove_one_charged_spectra=False, columns=None):
        cs = ['Protein', 'Hyperscore', 'decoy', 'taxID']
        if columns:
            cs = cs + columns
        reduced_df = ReferenceWriter.read_csv_with_generic_function(file, cs, remove_one_charged_spectra)
        fdr_pos_result, number_psm_result, number_decoy_result, double_spectra_result, score_last_item_result = PSM_FDR.determine_FDR_position(reduced_df, fdr)
        return number_psm_result, reduced_df[0:fdr_pos_result]

    def get_df_in_fdr_without_decoy_rows(self, file, fdr, remove_one_charged_spectra=True, columns=None):
        df = self.get_psm_and_df_in_fdr(file, fdr, remove_one_charged_spectra, columns)[1]
        df = df[get_hit_rows2(df.decoy)]
        return df

    def load_dfs(self):
        result_dfs= {}
        for level in self.level_to_tsv_path_dict.keys():
            tsv = self.level_to_tsv_path_dict[level]
            df_in_fdr = self.get_df_in_fdr_without_decoy_rows(tsv, self.fdr, columns=[f'taxID_{level}'])
            df_in_fdr.Protein = df_in_fdr.Protein.apply(lambda acc_set: {acc.split('|')[1] for acc in acc_set})
            result_dfs[level] = df_in_fdr
        return result_dfs

    def get_acc2taxid_dict(self, all_accs):
        final_accs = set()
        for acc in all_accs:
            try:
                final_accs.add(acc.split('|')[1])
            except:
                final_accs.add(acc)
        acc2tax_reader=ReadAccTaxon("/home/jules/Documents/Metaproteomics/databases/databases_tax2proteome/", "uniprot")
        acc_to_taxid_dict = acc2tax_reader.read_acc2tax(final_accs)
        acc_to_taxid_dict = {key: int(taxid) for key, taxid in acc_to_taxid_dict.items()}
        return acc_to_taxid_dict


    #um PSMs mit weniger als threshold nb PSMs zu entfernen muss acca zu taxid passende listen erstellen
    def get_taxid_of_acc(self, acc, level, taxon_graph, acc_to_taxid_dict):
        if 'REVERSED' in acc:
            return "DECOY"
        try:
            taxid = taxon_graph.find_level_up(acc_to_taxid_dict[acc], level)
        except KeyError:
            taxid = "CRAP"
        return taxid

    def get_taxids_of_accs(self, acc_list, level, taxon_graph, acc_to_taxid_dict):
        taxid_list = []
        for acc in acc_list:
            taxid_list.append(self.get_taxid_of_acc(acc, level, taxon_graph, acc_to_taxid_dict))
        return taxid_list

    def sort_taxid_and_acc_in_df(self, df_dict, acc_to_taxid_dict):
        for level, df_unsorted in df_dict.items():
            df_sorted = df_unsorted.copy(deep=True)
            df_sorted.Protein = df_sorted.Protein.apply(lambda acc_set: sorted(list(acc_set)))
            df_sorted[f"taxID_{level}"] = df_sorted.Protein.apply(lambda acc_list: self.get_taxids_of_accs(acc_list, level, self.taxon_graph,
                                                                                             acc_to_taxid_dict))
            df_dict[level]=df_sorted
        return df_dict

    def filter_all_protein_accs_from_dfs(self, df_dict):
        accs_set = set()
        for level, df in df_dict.items():
            accs_set= accs_set.union({item for sublist in df.Protein for item in sublist})
        return accs_set

    # functions for reducing df to PSMs wit at least x numbers of matching PSMs
    #nb_peptides 1 - 4
    def add_tax_information(self, acc_list, acc2tax_dict, level):
        tax_list=[]
        for acc in acc_list:
            if 'REVERSED' in acc:
                tax_list.append('DECOY')
            else:
                try:
                    tax_list.append(self.taxon_graph.find_level_up(acc2tax_dict[acc], level))
                except KeyError:
                    tax_list.append('CRAP')
        return tax_list


    def count_all_identified_spectra_per_taxa(self, df, level, taxa_list):
        taxon_to_spectra_count_dict={}
        for taxon in taxa_list:
            taxon_to_spectra_count_dict[taxon] = len(set(df[get_taxon_rows(df[f"taxID_{level}"], taxon)].Title))
        return taxon_to_spectra_count_dict

    def remove_all_accs_with_less_then_x_peptides(self, df, level, nb_peptides, acc2tax_dict):
        acc_to_peptide_dict=defaultdict(set)
        accs_to_remove= set()
        accs_to_keep = set()
        for index, row in df.iterrows():
            for acc in row['Protein']:
                if acc!="DECOY" and acc!="CRAP":
                    acc_to_peptide_dict[acc].add(row["Peptide"])
        for acc, pep_set in acc_to_peptide_dict.items():
            if len(pep_set) < nb_peptides:
                accs_to_remove.add(acc)
            else:
                accs_to_keep.add(acc)
        df_no_x_hits=df[remove_acc_row(df.Protein, accs_to_remove)]
        df_no_x_hits.Protein = df_no_x_hits.Protein.apply(lambda acc_list: set(acc_list).intersection(accs_to_keep))
        df_no_x_hits.Protein = df_no_x_hits.Protein.apply(lambda acc_set: sorted(list(acc_set)))
        df_no_x_hits = df_no_x_hits[remove_empty_rows(df_no_x_hits.Protein)]
        df_no_x_hits[f"taxID_{level}"] = df_no_x_hits.Protein.apply(lambda acc_list: self.add_tax_information(acc_list,
                                                                                                    acc2tax_dict, level))
        return df_no_x_hits

    def write_header(self, levels):
        fields = []
        for level in levels:
            fields.append(f"name {level}")
        for level in levels:
            fields.append(f"spectra count {level}")
            if level == "species":
                for group_level in self.level_to_taxa_list_dict.keys():
                    fields.append(f"spectra count species - {group_level}")
        for level in levels:
            fields.append(f"unique spectra count {level}")

        with open(self.path_to_result, "w") as out:
            out.write("\t".join(fields))
            out.write("\n")

    def write_results(self, all_spectra_dict, unique_spectra_dict, levels):
        """
        param all_spectra_dict: {min_nb_psms:{taxon_str:{level:spectra_set}}}, spectra set includes all spectra identified for taxon
        param unique_spectra_dict: {min_nb_psms:{taxon_str:{level:spectra_set}}}, spectra set includes only spectra unique identified for taxon
            """
        with open(self.path_to_result, "a") as out:
            for min_nb_of_psms in all_spectra_dict.keys():
                out.write(f"Minimum number of PSMS per Protein to be counted: {min_nb_of_psms}\n")
                # row: name species, name genus, name family, spectra count species, ..., unique spectra count species, ...
                for species_taxon in self.species_taxa_list:
                    row = [self.taxon_graph.get_scientific_name(species_taxon)]
                    # names --> corect
                    for level in levels[1:]:
                        level_taxon = self.taxon_graph.find_level_up(species_taxon, level)
                        row.append(self.taxon_graph.get_scientific_name(level_taxon))
                    #spectral count
                    for level in levels:
                        row.append(str(len(all_spectra_dict[min_nb_of_psms][species_taxon][level])))
                        #add spectra count multiple related species for interesting level together
                        if level == "species":
                            for lev in levels[1:]:
                                try:
                                    if any(species_taxon in x for x in self.level_to_taxa_list_dict[lev]):
                                        taxa_list = [tax_list for tax_list in self.level_to_taxa_list_dict[lev] if
                                                     species_taxon in tax_list][0]
                                        spectra_set = set()
                                        for tax in taxa_list:
                                            spectra_set = spectra_set.union(all_spectra_dict[min_nb_of_psms][tax]['species'])
                                        row.append(str(len(spectra_set)))
                                    else:
                                        row.append("")
                                # more levels then level in level_taxa-groups
                                except KeyError:
                                    continue
                    #unique spectral count
                    for level in levels:
                        row.append(str(len(unique_spectra_dict[min_nb_of_psms][species_taxon][level])))
                    out.write("\t".join(row))
                    out.write("\n")
                # taxa groups
                if self.name_to_taxa_group_dict:
                    for name, taxa_list in self.name_to_taxa_group_dict.items():
                        row = [name]*len(levels)
                        for level in levels:
                            row.append(str(len(all_spectra_dict[min_nb_of_psms][name][level])))
                            #add spectra count multiple related species for interesting level together
                            if level == "species":
                                for lev in levels[1:]:
                                    row.append("")
                        for level in levels:
                            row.append(str(len(unique_spectra_dict[min_nb_of_psms][name][level])))
                        out.write("\t".join(row))
                        out.write("\n")
                #last line
                out.write(f"all spectra\t")
                for i in range(len(levels)-1):
                    out.write("\t")
                for level in levels:
                    out.write(f"{all_spectra_dict[min_nb_of_psms]['all'][level]}\t")
                    if level == 'species':
                        for i in range(len(levels)-1):
                            out.write("\t")
                out.write("\n\n")

    def extract_information_from_df(self, level_to_df_dict, acc_to_taxid_dict, min_nb_psms=4):
        """
        return all_spectra_dict: {min_nb_psms:{taxon_str:{level:spectra_set}}}, spectra set includes all spectra identified for taxon
        e.g. {2: {'562': {'species': {'Run1_U1_2000ng.41650.41650.2','Run1_U1_2000ng.76581.76581.3',...)
        return unique_spectra_dict: {min_nb_psms:{taxon_str:{level:spectra_set}}}, spectra set includes only spectra unique identified for taxon
        """
        # min_nb_of_psms : {level: {taxon: spectra_count}}
        # level_to_taxa_tuple_list_dict : {level: (name, taxa list species, taxa list level, % protein)}
        all_spectra_dict = {}
        unique_spectra_dict = {}
        for min_nb_of_psms in range(1,min_nb_psms+1):
            print(f"Processing results for minimum {min_nb_of_psms} PSMS per Protein to be counted.")
            all_spectra_dict[min_nb_of_psms] = defaultdict(dict)
            unique_spectra_dict[min_nb_of_psms] = defaultdict(dict)
            #create df with all accs removed with less than min_nb_of_psms spectra hits per protein
            for level in level_to_df_dict.keys():
                df_no_x_hits=self.remove_all_accs_with_less_then_x_peptides(level_to_df_dict[level].copy(deep=True),
                                                                       level, min_nb_of_psms, acc_to_taxid_dict)
                df_no_x_hits_unique=df_no_x_hits[get_all_exclusive_rows(df_no_x_hits[f"taxID_{level}"])]
                print(f"dfs cleaned from unique hits for level {level}")
                for taxon in self.species_taxa_list:
                    # all species taxa, determine genus ans family taxa
                    taxon_level = self.taxon_graph.find_level_up(taxon, level)
                    all_spectra_dict[min_nb_of_psms][taxon][level] = \
                        set(df_no_x_hits[get_taxa_rows(df_no_x_hits[f"taxID_{level}"], taxon_level)].Title)
                    unique_spectra_dict[min_nb_of_psms][taxon][level] = \
                        set(df_no_x_hits_unique[get_taxa_rows(df_no_x_hits_unique[f"taxID_{level}"], taxon_level)].Title)

                    # virus all taxa counted together
                if self.name_to_taxa_group_dict:
                    for group_name, taxa_group in self.name_to_taxa_group_dict.items():
                        taxa_level = [self.taxon_graph.find_level_up(taxon, level) for taxon in taxa_group]
                        all_spectra_dict[min_nb_of_psms][group_name][level] = \
                            set(df_no_x_hits[get_taxa_rows(df_no_x_hits[f"taxID_{level}"], taxa_level)].Title)
                        unique_spectra_dict[min_nb_of_psms][group_name][level] = \
                            set(df_no_x_hits_unique[get_taxa_rows(df_no_x_hits_unique[f"taxID_{level}"], taxa_level)].Title)
                # all spectra
                all_spectra_dict[min_nb_of_psms]['all'][level] = len(set(df_no_x_hits.Title))
        return all_spectra_dict, unique_spectra_dict

    def analyse_identification_files(self):
        level_to_df_dict = self.load_dfs()
        print("Dataframes loaded")
        all_accs_set = self.filter_all_protein_accs_from_dfs(level_to_df_dict)
        acc_to_taxid_dict = self.get_acc2taxid_dict(all_accs_set)
        print("Protein accessions and taxon IDs loaded.")
        level_to_df_dict = self.sort_taxid_and_acc_in_df(level_to_df_dict, acc_to_taxid_dict)
        print("Dataframes sorted.")
        all_spectra_dict, unique_spectra_dict = self.extract_information_from_df(level_to_df_dict, acc_to_taxid_dict)
        self.write_header(level_to_df_dict.keys())
        self.write_results(all_spectra_dict, unique_spectra_dict, list(level_to_df_dict.keys()))
        print(f"Results stored in {self.path_to_result}")


def main():
    path_to_kleiner_results_bachelor = "/home/jules/Documents/Metaproteomics/Tax2Proteome/benchmarking/results_searchgui_xtandem_analyzer_bachelor_thesis/"
    path_to_tanca_results_bachelor = "/home/jules/Documents/Metaproteomics/Tax2Proteome/benchmarking/results_searchgui_xtandem_analyzer_bachelor_thesis/9MM_FASP/x_tandem_tsv"
    uniprot_nr_reduced_tsv = {
        'subspecies': path_to_kleiner_results_bachelor + "/uniprot_kleiner/x_tandem_tsv/Run1_U1_2000ng_uniprot_subspecies.t.xml_new_reduced.tsv",
        'species': path_to_kleiner_results_bachelor + "/uniprot_kleiner/x_tandem_tsv/Run1_U1_2000ng_uniprot_species_nr.t.xml_new_reduced.tsv",
        'genus': path_to_kleiner_results_bachelor + "/uniprot_kleiner/x_tandem_tsv/Run1_U1_2000ng_uniprot_genus_nr.t.xml_new_reduced.tsv",
        'family': path_to_kleiner_results_bachelor + "/uniprot_kleiner/x_tandem_tsv/Run1_U1_2000ng_uniprot_family_nr.t.xml_new_reduced.tsv",
    }

    tanca_uniprot_tsv = {
        'species': path_to_tanca_results_bachelor + "/9MM_FASP_uniprot_Tanca_species_nr.t.xml_new_reduced.tsv",
        'genus': path_to_tanca_results_bachelor + "/9MM_FASP_uniprot_Tanca_genus.t.xml_new_reduced.tsv",
        'family': path_to_tanca_results_bachelor + "/9MM_FASP_uniprot_Tanca_family_nr.t.xml_new_reduced.tsv",
        'order': path_to_tanca_results_bachelor + "/9MM_FASP_uniprot_Tanca_order_nr.t.xml_new_reduced.tsv"
    }
    Kleiner_species_taxa = [
        ("Escherichia coli", [562], [562], 5.788),
        ("Salmonella enterica", [28901], [28901], 33.773),
        ("Bacillus subtilis", [1423], [1423], 0.788),
        ("Staphylococcus aureus", [1280], [1280], 2.605),
        ("Desulfovibrio vulgaris", [881], [881], 0.946),
        ("Thermus thermophilus", [274], [274], 1.68),
        ("Chlamydomonas reinhardtii", [3055], [3055], 3.996),
        ("Paracoccus denitrificans", [266], [266], 0.922),
        ("Nitrososphaera viennensis", [1034015], [1034015], 0.819),
        ("Stenotrophomonas maltophilia", [40324], [40324], 8.021),
        ("Altermonas macleodii", [28108], [28108], 0.954),
        ("Chromobacterium violaceum", [536], [536], 1.259),
        ("Paraburkholderia xenovorans", [36873], [36873], 0.433),
        ("Cupriavidus metallidurans", [119219], [119219], 15.519),
        ("Nitrosomonas ureae", [44577], [44577], 0.543),
        ("Nitrosomonas europaea", [915], [915], 0.082),
        ("Nitrosospira multiformis", [1231], [1231], 0.209),
        ("Agrobacterium fabrum", [1176649], [1176649], 5.647),
        ("Rhizobium leguminosarum", [384], [384], 3.172),
        ("Pseudomonas fluorescens", [294], [294], 6.696),
        ("Pseudomonas furukawaii", [1149133], [1149133], 1.165),
        ("Pseudomonas sp. ATCC 13867", [1294143], [1294143], 2.871)
    ]
    Kleiner_groups = {'virus': [10754, 101570, 1985310, 329852, 1977402]}

    Kleiner_level_groups = {'genus': [[44577, 915], [294, 1149133, 1294143]],
                                 'family': [[28901, 562], [36873, 119219],  [44577, 915, 1231],  [1176649, 384],
                                            [294, 1149133, 1294143]]
                                 }
    Tanca_species_taxa = [
        ("Pasteurella multocida", [747], 11.1), ("Brevibacillus laterosporus", [1465], 11.1),
        ("Lactobacillus acidophilus", [1579], 11.1), ("Lactobacillus casei", [1582], 11.1),
        ("Pediococcus pentosaceus", [1255], 11.1), ("Enterococcus faecalis", [1351], 11.1),
        ("Rhodotorula glutinis", [5535], 11.1), ("Saccharomyces cerevisiae", [4932], 11.1),
        ("Escherichia coli", [562], 11.1)
    ]
    Tanca_level_groups = {'genus': [[1579,1582]], 'family': [[1579,1582,1255]]}

    path_to_taxdump = "/home/jules/Documents/Metaproteomics/databases/databases_tax2proteome/taxdump.tar.gz"
    path_to_kleiner_result = f"/home/jules/Documents/Metaproteomics/Tax2Proteome/benchmarking/results_analysis_min_nb_psms/Kleiner_all_spectra_more_than_x_hits_2.tsv"
    Kleiner_species_taxa_list = [tax_tuple[1][0] for tax_tuple in Kleiner_species_taxa]
    level_to_tsv_dict= {k: df for k,df in uniprot_nr_reduced_tsv.items() if k!='subspecies'}
    MinNbAnalyserKleiner =  defined_numer_of_psms_per_protein(path_to_taxdump, level_to_tsv_dict, Kleiner_species_taxa_list,
                                                              path_to_kleiner_result, Kleiner_level_groups, Kleiner_groups)
    MinNbAnalyserKleiner.analyse_identification_files()

    path_to_tanca_result = f"/home/jules/Documents/Metaproteomics/Tax2Proteome/benchmarking/results_analysis_min_nb_psms/Tanca_all_spectra_more_than_x_hits_2.tsv"
    species_taxa_list_tanca = [tax_tuple[1][0] for tax_tuple in Tanca_species_taxa]

    MinNbAnalyserTanca =  defined_numer_of_psms_per_protein(path_to_taxdump, tanca_uniprot_tsv, species_taxa_list_tanca,
                                                            path_to_tanca_result, Tanca_level_groups)
    MinNbAnalyserTanca.analyse_identification_files()


if __name__ == "__main__":
    main()
