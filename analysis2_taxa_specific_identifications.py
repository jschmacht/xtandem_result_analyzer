import pandas as pd
import numpy as np
from pathlib import Path
from create_reference_from_tsv_and_pepxml import ReferenceWriter
from create_PSM_df import PSM_FDR
from collections import defaultdict
from ReadAccTaxon import ReadAccTaxon
# load taxon graph
from TaxonGraph import TaxonGraph
taxon_graph = TaxonGraph()
taxon_graph.create_graph("/home/jules/Documents/Metaproteomics/databases/databases_tax2proteome/taxdump.tar.gz")

path_to_bachelor_results = "/home/jules/Documents/Metaproteomics/Tax2Proteome/benchmarking/results_searchgui_xtandem_analyzer_bachelor_thesis/"

uniprot_nr_reduced_tsv = {
    'subspecies': path_to_bachelor_results + "/uniprot_kleiner/x_tandem_tsv/Run1_U1_2000ng_uniprot_subspecies.t.xml_new_reduced.tsv",
    'species': path_to_bachelor_results + "/uniprot_kleiner/x_tandem_tsv/Run1_U1_2000ng_uniprot_species_nr.t.xml_new_reduced.tsv",
    'genus': path_to_bachelor_results + "/uniprot_kleiner/x_tandem_tsv/Run1_U1_2000ng_uniprot_genus_nr.t.xml_new_reduced.tsv",
    'family': path_to_bachelor_results + "/uniprot_kleiner/x_tandem_tsv/Run1_U1_2000ng_uniprot_family_nr.t.xml_new_reduced.tsv",
}

Kleiner_2 = [("Desulfovibrio vulgaris", [881]), ("Nitrosomonas ureae", [44577]),("Thermus thermophilus", [274]),
             ("Stenotrophomonas maltophilia", [40324]), ("Chlamydomonas reinhardtii", [3055]),  ("Nitrososphaera viennensis", [1034015])]

Kleiner_3 = [("Enterobacteriaceae", [28901, 562],  [543]), ("Nitrosomonadaceae", [384, 1176649],  [82115]),
             ("Rhizobium", [44577, 1231, 915],  [206379]), ("Pseudomonas", [294, 1149133, 1294143], [135621]),
             ("viruses", [10754, 101570, 1283336, 12022, 1977402], [10744, 10699, 10662, 11989]),
             ("paraburkholderia", [36873, 119219], [119060]), ("altermonas macleodii", [28108], [72275]),
             ("chromobacterium violaceum", [536], [1499392]), ("Paracoccus denitrificans", [266], [31989]),
             ("staphylococcus aureus", [1280], [90964]), ("bacillus subtilis", [1423], [186817])]

def get_hit_rows2(decoy_column):
    return [True if d_set in  [{True, False}, {False}] else False for d_set in decoy_column]

def get_psm_and_df_in_fdr(file, fdr, remove_one_charged_spectra=False, columns=None):
    cs = ['Protein', 'Hyperscore', 'decoy', 'taxID']
    if columns:
        cs = cs + columns
    reduced_df = ReferenceWriter.read_csv_with_generic_function(file, cs, remove_one_charged_spectra)
    fdr_pos_result, number_psm_result, number_decoy_result, double_spectra_result, score_last_item_result = PSM_FDR.determine_FDR_position(reduced_df, fdr)
    return number_psm_result, reduced_df[0:fdr_pos_result]

def get_df_in_fdr_without_decoy(file, fdr, remove_one_charged_spectra=True, columns=None):
    df = get_psm_and_df_in_fdr(file, fdr, remove_one_charged_spectra, columns)[1]
    df = df[get_hit_rows2(df.decoy)]
    return df

def get_taxa_rows(column, taxID):
    if type(taxID)==int:
        return [True if taxID in t_set else False for t_set in column]
    elif type(taxID)==list:
        return [True if len(set(taxID).intersection(set(t_set)))>0 else False for t_set in column]

def get_exclusive_rows(column, taxid):
    if type(taxid)==int:
        return [True if {taxid} == t_set else False for t_set in column]
    elif type(taxid)==list:
        return [True if len(set(t_set).difference(taxid))==0 else False for t_set in column]

def remove_acc_row(column, acc_to_remove_set):
    return [False if len(set(accs).difference(acc_to_remove_set))==0 else True for accs in column]

def remove_empty_rows(column):
    return [False if len(accs)==0 else True for accs in column]

def remove_all_accs_with_only_one_peptide(df):
    acc_to_peptide_dict=defaultdict(set)
    accs_to_remove= set()
    accs_to_keep = set()
    for index, row in df.iterrows():
        for acc in row['Protein']:
            acc_to_peptide_dict[acc].add(row["Peptide"])
    for acc, pep_set in acc_to_peptide_dict.items():
        if len(pep_set)<2:
            accs_to_remove.add(acc)
        else:
            accs_to_keep.add(acc)
    df_no_one_hits=df[remove_acc_row(df.Protein, accs_to_remove)]
    df_no_one_hits.Protein = df_no_one_hits.Protein.apply(lambda acc_list: set(acc_list).intersection(accs_to_keep))
    df_no_one_hits = df_no_one_hits[remove_empty_rows(df_no_one_hits.Protein)]
    return df_no_one_hits

def remove_unrelated_accs_and_taxa(df, taxids, level, acc_2_taxid_dict):
    df.Protein = df.Protein.apply(lambda acc_list: [acc for acc in acc_list if \
                                                    get_taxid_of_acc(acc, level, taxon_graph, acc_2_taxid_dict) in taxids])

    df[f"taxID_{level}"] = df.Protein.apply(lambda acc_list: [get_taxid_of_acc(acc, level, taxon_graph, acc_2_taxid_dict) for acc in acc_list])
    return df


def taxon_exclusive_spectra(df, df_no_one_hits, column_name, taxon, level, acc_2_taxid_dict):
    df_all_taxa_spec_identified_spectra = df[get_taxa_rows(df[column_name], taxon)]
    if df_all_taxa_spec_identified_spectra.empty:
        df_all_taxa_spec_identified_spectra = pd.DataFrame(columns=["Title", "Protein", f"taxID_species"])
    all_identified_spectra = set(df_all_taxa_spec_identified_spectra.Title)
    df_exclusive_spectra = df_all_taxa_spec_identified_spectra[get_exclusive_rows(df_all_taxa_spec_identified_spectra[column_name], taxon)]
    if df_exclusive_spectra.empty:
        df_exclusive_spectra = pd.DataFrame(columns=["Title", "Protein", f"taxID_species"])
    exclusive_spectra = set(df_exclusive_spectra.Title)

    df_all_taxa_spec_identified_spectra_no_one_hits = df_no_one_hits[get_taxa_rows(df_no_one_hits[column_name], taxon)]
    df_all_taxa_spec_identified_spectra_no_one_hits = remove_unrelated_accs_and_taxa(df, taxon, level, acc_2_taxid_dict)
    if df_all_taxa_spec_identified_spectra_no_one_hits.empty:
        df_all_taxa_spec_identified_spectra_no_one_hits = pd.DataFrame(columns=["Title", "Protein", f"taxID_species"])
    all_identified_spectra_no_one_hits = set(df_all_taxa_spec_identified_spectra_no_one_hits.Title)
    df_exclusive_spectra_no_one_hits = df_all_taxa_spec_identified_spectra_no_one_hits[get_exclusive_rows(df_all_taxa_spec_identified_spectra_no_one_hits[column_name], taxon)]
    if df_exclusive_spectra_no_one_hits.empty:
        df_exclusive_spectra_no_one_hits = pd.DataFrame(columns=["Title", "Protein", f"taxID_species"])
    exclusive_spectra_no_one_hits = set(df_exclusive_spectra_no_one_hits.Title)
    return all_identified_spectra, all_identified_spectra_no_one_hits, exclusive_spectra, exclusive_spectra_no_one_hits


def get_acc2taxid_dict(df_in_fdr_uniprot_species, df_in_fdr_uniprot_family):
    all_accs={item for sublist in df_in_fdr_uniprot_species.Protein for item in sublist}
    all_accs=all_accs.union({item for sublist in df_in_fdr_uniprot_family.Protein for item in sublist})
    final_accs = set()
    for acc in all_accs:
        try:
            final_accs.add(acc.split('|')[1])
        except:
            final_accs.add(acc)
    print(len(final_accs))
    acc2tax_reader=ReadAccTaxon("/home/jules/Documents/Metaproteomics/databases/databases_tax2proteome/", "uniprot")
    acc_to_taxid_dict = acc2tax_reader.read_acc2tax(final_accs)
    acc_to_taxid_dict = {key: int(taxid) for key, taxid in acc_to_taxid_dict.items()}
    return acc_to_taxid_dict

def get_taxids_of_accs(acc_list, level, taxon_graph, acc_to_taxid_dict):
    taxid_list = []
    for acc in acc_list:
        taxid_list.append(get_taxid_of_acc(acc, level, taxon_graph, acc_to_taxid_dict))
    return taxid_list

def get_taxid_of_acc(acc, level, taxon_graph, acc_to_taxid_dict):
    try:
        taxid = taxon_graph.find_level_up(acc_to_taxid_dict[acc], level)
    except KeyError:
        taxid = "CRAP"
    return taxid

def sort_taxid_and_acc_in_df(df, level, acc_to_taxid_dict):
    df.Protein = df.Protein.apply(lambda acc_set: sorted(list(acc_set)))
    df[f"taxID_{level}"] = df.Protein.apply(lambda acc_list: get_taxids_of_accs(acc_list, level, taxon_graph, acc_to_taxid_dict))
    return df

def count_taxa(df_column):
    taxa_count_dict = {}
    for tax_set in df_column:
        for taxon in tax_set:
            if taxon in taxa_count_dict:
                taxa_count_dict[taxon]+= 1
            else:
                taxa_count_dict[taxon]=1
    return taxa_count_dict

def count_accs(acc_list):
    acc_list.sort()
    count_dict={1:0, 2:0, 3:0, '>3':0}
    try:
        acc_before = acc_list[0]
        count=0
        for acc in acc_list:
            if acc == acc_before:
                count+=1
            else:
                if count > 3:
                    count_dict['>3'] = count_dict['>3']+1
                else:
                    count_dict[count] = count_dict[count]+1
                count=1
            acc_before = acc
        return count_dict
    except IndexError:
        return count_dict

def get_taxon_specific_df(df, level, taxids, acc_2_taxid_dict):
    df = df[get_taxa_rows(df[f"taxID_{level}"], taxids)]
    df.Protein = df.Protein.apply(lambda acc_list: [acc for acc in acc_list if \
                                                   get_taxid_of_acc(acc, level, taxon_graph, acc_2_taxid_dict) in taxids])

    df[f"taxID_{level}"] = df.Protein.apply(lambda acc_list: [get_taxid_of_acc(acc, level, taxon_graph, acc_2_taxid_dict) for acc in acc_list])
    #  df.Protein = df.Protein.apply(lambda acc_list:[x for i, x in enumerate(acc_list) if i in get_index_list(df[f"taxID_{level}"], taxids)])
    # df[f"taxID_{level}"]= df[f"taxID_{level}"].apply(lambda tax_list:[taxid for taxid in tax_list if taxid in taxids])
    #remove rows with no accs
    df = df[remove_empty_rows(df.Protein)]
    if df.empty:
        df = pd.DataFrame(columns=["Title", "Protein", f"taxID_{level}"])
    return df


def describe_difference_between_species_and_family_identification(name, taxids_species, taxid_family,
                                                                  df_in_fdr_species, df_in_fdr_family,
                                                                  df_in_fdr_species_no_one_hits,
                                                                  df_in_fdr_family_no_one_hits, acc_to_taxid_dict):
    print(name)
    for df_spec, df_fam, name in [(df_in_fdr_species, df_in_fdr_family, 'FULL'), (df_in_fdr_species_no_one_hits, df_in_fdr_family_no_one_hits, "REDUCED")]:
        if name=='FULL':
            print("RESULTS FOR FULL DF")
        else:
            print()
            print("RESULS FOR DF without ONE HITS")
        # keep only specific accs
        df_in_fdr_species_tax_acc = df_spec.copy(deep=True)
        df_in_fdr_family_tax_acc = df_fam.copy(deep=True)
        # keep acc with correct assigned tax
        df_in_fdr_species_tax_acc = get_taxon_specific_df(df_in_fdr_species_tax_acc, 'species', taxids_species, acc_to_taxid_dict)
        df_in_fdr_family_tax_acc = get_taxon_specific_df(df_in_fdr_family_tax_acc, 'family', taxid_family, acc_to_taxid_dict)

        #print(df_in_fdr_family_186817_acc.head(10))
        # all accs
        try:
            accs_spec = sorted([item for sublist in df_in_fdr_species_tax_acc.Protein for item in sublist])
            accs_fam = sorted([item for sublist in df_in_fdr_family_tax_acc.Protein for item in sublist])
        except:
            print('error: empty df? ', df_in_fdr_species_tax_acc)

        print(f"identified spectra species for taxon {taxids_species}: {len(set(df_in_fdr_species_tax_acc.Title))}\nidentified spectra family for taxon {taxid_family}: {len(set(df_in_fdr_family_tax_acc.Title))}")
        print(f"identified accs species for taxon {taxids_species}: {len(set(accs_spec))}\nidentified accs family for taxon {taxid_family}: {len(set(accs_fam))}")

        #count
        count_dict_spec = count_accs(accs_spec)
        print(f"Anzahl Protein-Identifikationen auf Species Level für taxon {taxids_species}\n\
        mit 1 identifikation:{count_dict_spec[1]}\n\
        mit 2 identifikationen: {count_dict_spec[2]},\n\
        mit 3 identifikationen: {count_dict_spec[3]}, \n\
        mit mehr als 3 Identifikationen: {count_dict_spec['>3']}")

        count_dict_fam = count_accs(accs_fam)
        print(f"Anzahl Protein-Identifikationen auf Family Level für taxon {taxid_family}\n\
        mit 1 identifikation:{count_dict_fam[1]}\n\
        mit 2 identifikationen: {count_dict_fam[2]},\n\
        mit 3 identifikationen: {count_dict_fam[3]}, \n\
        mit mehr als 3 Identifikationen: {count_dict_fam['>3']}")

        print()

        print("Eigenschaften der neu identifizierten Spektra:")
        df_family_tax_fam_not_in_species_df=df_in_fdr_family_tax_acc[~df_in_fdr_family_tax_acc.Title.isin(df_in_fdr_species_tax_acc.Title)]

        print(f"all new identified spectra family: {len(set(df_family_tax_fam_not_in_species_df.Title))}")
        df_only_one_taxa_tax_fam = df_family_tax_fam_not_in_species_df[get_taxa_rows(df_family_tax_fam_not_in_species_df.taxID_family, taxid_family)]
        df_only_one_taxa_tax_fam = df_family_tax_fam_not_in_species_df[get_taxa_rows(df_family_tax_fam_not_in_species_df.taxID_family, taxid_family)]
        taxa_count_dict_tax_fam = count_taxa(df_family_tax_fam_not_in_species_df.taxID_family)
        if len(df_only_one_taxa_tax_fam) != 0:
            print(f"identified exclusive for taxon {taxid_family}: {len(set(df_only_one_taxa_tax_fam.Title))}\nnot exclusive identified for taxon {taxid_family}, taxid to PSM number dict: \n{taxa_count_dict_tax_fam}")
        if len(df_family_tax_fam_not_in_species_df) != 0:
            df_family_tax_fam_in_species_df = df_family_tax_fam_not_in_species_df[df_family_tax_fam_not_in_species_df.Title.isin(df_in_fdr_uniprot_species.Title)]
            if len(df_family_tax_fam_in_species_df) != 0:
                print(f"Anzahl Spekten bereits vorher für andere Spezien identifiziert: {len(set(df_family_tax_fam_in_species_df.Title))}")

def load_spec_and_fam_df():
    fdr=0.05
    df_in_fdr_uniprot_species = get_df_in_fdr_without_decoy(uniprot_nr_reduced_tsv['species'], fdr, columns=['taxID_species'])
    df_in_fdr_uniprot_family =  get_df_in_fdr_without_decoy(uniprot_nr_reduced_tsv['family'], fdr, columns=['taxID_family'])
    df_in_fdr_uniprot_species.Protein = df_in_fdr_uniprot_species.Protein.apply(lambda acc_set: {acc.split('|')[1] \
                                                                                                 for acc in acc_set})
    df_in_fdr_uniprot_family.Protein = df_in_fdr_uniprot_family.Protein.apply(lambda acc_set: {acc.split('|')[1] \
                                                                                               for acc in acc_set})

    return df_in_fdr_uniprot_species, df_in_fdr_uniprot_family

def write_results(out, all_identified_spectra_spe, all_identified_spectra_no_one_hits_spe, exclusive_spectra_spe,
                  exclusive_spectra_no_one_hits_spe, all_identified_spectra_fam, all_identified_spectra_no_one_hits_fam,
                  exclusive_spectra_fam, exclusive_spectra_no_one_hits_fam, new_identified_spectra, spectra_in_species_df):
    out.write(f"all identified_spectra for taxa level species\t{len(all_identified_spectra_spe)}\n")
    out.write(f"all identified_spectra for taxa level family\t{len(all_identified_spectra_fam)}\n")
    out.write(f"exclusive spectra for taxa level species\t{len(exclusive_spectra_spe)}\n")
    out.write(f"exclusive spectra for taxa level family\t{len(exclusive_spectra_fam)}\n")
    out.write(f"intersection species and family exclusive spectra\t{len(exclusive_spectra_spe.intersection(exclusive_spectra_fam))}")
    out.write("\n")
    out.write(f"all identified_spectra for taxa level species no one hits\t{len(all_identified_spectra_no_one_hits_spe)}\n")
    out.write(f"all identified_spectra for taxa level family no one hits\t{len(all_identified_spectra_no_one_hits_fam)}\n")
    out.write(f"exclusive spectra for taxa level species no one hits\t{len(exclusive_spectra_no_one_hits_spe)}\n")
    out.write(f"exclusive spectra for taxa level family no one hits\t{len(exclusive_spectra_no_one_hits_fam)}\n")
    out.write(f"intersection species and family exclusive spectra no one hits\t{len(exclusive_spectra_no_one_hits_spe.intersection(exclusive_spectra_no_one_hits_fam))}\n")

    out.write(f"new_identified_spectra no one hits\t{len(new_identified_spectra)}\n")
    out.write(f"new_identified_spectra bereits in species_df no one hit: {len(spectra_in_species_df)}\n")
    out.write(f"spectra unknown in species df no one hits\t{len(new_identified_spectra) - len(spectra_in_species_df)}\n")
    out.write("\n")


def main():
    df_in_fdr_uniprot_species, df_in_fdr_uniprot_family = load_spec_and_fam_df()
    acc_to_taxid_dict = get_acc2taxid_dict(df_in_fdr_uniprot_species, df_in_fdr_uniprot_family)
    df_in_fdr_uniprot_species = sort_taxid_and_acc_in_df(df_in_fdr_uniprot_species, 'species', acc_to_taxid_dict)
    df_in_fdr_uniprot_family = sort_taxid_and_acc_in_df(df_in_fdr_uniprot_family, 'family', acc_to_taxid_dict)

    df_in_fdr_uniprot_species_no_one_hits = remove_all_accs_with_only_one_peptide(df_in_fdr_uniprot_species.copy(deep=True))
    df_in_fdr_uniprot_family_no_one_hits = remove_all_accs_with_only_one_peptide(df_in_fdr_uniprot_family.copy(deep=True))

    with open("/home/jules/Documents/Metaproteomics/Tax2Proteome/benchmarking/results_analysis/results_per_taxa.tsv", "w") as out:
        out.write("Analysis of identifed spectra\n")

    with open("/home/jules/Documents/Metaproteomics/Tax2Proteome/benchmarking/results_analysis/results_per_taxa.tsv", "a") as out:
        for name, taxon_species in Kleiner_2:
            taxon_family = taxon_graph.find_level_up(taxon_species[0], "family")
            out.write(f"{name}\t{taxon_species[0]}\t{taxon_family}\n")
            all_identified_spectra_spe, all_identified_spectra_no_one_hits_spe, exclusive_spectra_spe, \
            exclusive_spectra_no_one_hits_spe = taxon_exclusive_spectra(df_in_fdr_uniprot_species, df_in_fdr_uniprot_species_no_one_hits, "taxID_species", taxon_species, level, acc_2_taxid_dict)
            all_identified_spectra_fam, all_identified_spectra_no_one_hits_fam, exclusive_spectra_fam, exclusive_spectra_no_one_hits_fam = \
                taxon_exclusive_spectra(df_in_fdr_uniprot_family, df_in_fdr_uniprot_family_no_one_hits, "taxID_family", taxon_family, level, acc_2_taxid_dict)
            new_identified_spectra, spectra_in_species_df = spectra_identified_in_df(df_in_fdr_uniprot_species_no_one_hits,
                                                                                     df_in_fdr_uniprot_family_no_one_hits, taxon_species, taxon_family)
            write_results(out, all_identified_spectra_spe, all_identified_spectra_no_one_hits_spe, exclusive_spectra_spe,
                          exclusive_spectra_no_one_hits_spe, all_identified_spectra_fam, all_identified_spectra_no_one_hits_fam,
                          exclusive_spectra_fam, exclusive_spectra_no_one_hits_fam, new_identified_spectra, spectra_in_species_df)

        for name, taxon_species, taxon_family in Kleiner_3:
            out.write(f"{name}\t{taxon_species[0]}\t{taxon_family}\n")
            all_identified_spectra_spe, all_identified_spectra_no_one_hits_spe, exclusive_spectra_spe, exclusive_spectra_no_one_hits_spe = \
                taxon_exclusive_spectra(df_in_fdr_uniprot_species, df_in_fdr_uniprot_species_no_one_hits, "taxID_species", taxon_species)
            all_identified_spectra_fam, all_identified_spectra_no_one_hits_fam, exclusive_spectra_fam, exclusive_spectra_no_one_hits_fam = \
                taxon_exclusive_spectra(df_in_fdr_uniprot_family, df_in_fdr_uniprot_family_no_one_hits, "taxID_family", taxon_family)
            new_identified_spectra, spectra_in_species_df = spectra_identified_in_df(df_in_fdr_uniprot_species_no_one_hits,
                                                                                     df_in_fdr_uniprot_family_no_one_hits, taxon_species, taxon_family)
            write_results(out, all_identified_spectra_spe, all_identified_spectra_no_one_hits_spe, exclusive_spectra_spe,
                          exclusive_spectra_no_one_hits_spe, all_identified_spectra_fam, all_identified_spectra_no_one_hits_fam,
                          exclusive_spectra_fam, exclusive_spectra_no_one_hits_fam, new_identified_spectra, spectra_in_species_df)
