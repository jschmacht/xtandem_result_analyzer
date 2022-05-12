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

Kleiner_2 = [("Desulfovibrio vulgaris", [881]), ("Thermus thermophilus", [274]),
             ("Stenotrophomonas maltophilia", [40324]), ("Chlamydomonas reinhardtii", [3055]),
             ("Nitrososphaera viennensis", [1034015]), ("altermonas macleodii", [28108]),
             ("chromobacterium violaceum", [536]), ("Paracoccus denitrificans", [266]),
             ("staphylococcus aureus", [1280]), ("bacillus subtilis", [1423])]

Kleiner_3 = [("Enterobacteriaceae", [28901, 562],  [543]), ("Nitrosomonadaceae", [44577, 915, 1231],  [206379]),
             ("Rhizobium", [1176649, 384],  [82115]), ("Pseudomonas", [294, 1149133, 1294143], [135621]),
             ("viruses", [10754, 101570, 1985310, 329852, 1977402], [10744, 10699, 10662, 11989, 10860]),
             ("paraburkholderia", [36873, 119219], [119060])]

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

def get_acc_rows(column, acc):
    return [True if acc in acc_list else False for acc_list in column ]

def get_taxa_rows(column, taxID):
    if type(taxID)==int:
        return [True if taxID in t_set else False for t_set in column]
    elif type(taxID)==list:
        return [True if len(set(taxID).intersection(set(taxa_set)))>0 else False for taxa_set in column]

def get_exclusive_rows(column, taxid):
    if type(taxid)==int:
        return [True if {taxid} == t_set else False for t_set in column]
    elif type(taxid)==list:
        return [True if len(set(t_set).difference(taxid))==0 else False for t_set in column]

def remove_acc_row(column, acc_to_remove_set):
    return [False if len(set(accs).difference(acc_to_remove_set))==0 else True for accs in column]

def remove_empty_rows(column):
    return [False if len(accs)==0 else True for accs in column]

def add_tax_information(acc_list, acc2tax_dict, level):
    tax_list=[]
    for acc in acc_list:
        try:
            tax_list.append(taxon_graph.find_level_up(acc2tax_dict[acc], level))
        except KeyError:
            tax_list.append('CRAP')
    return tax_list

def remove_all_accs_with_less_then_x_peptides(df, level, nb_peptides, acc2tax_dict):
    acc_to_peptide_dict=defaultdict(set)
    accs_to_remove= set()
    accs_to_keep = set()
    for index, row in df.iterrows():
        for acc in row['Protein']:
            acc_to_peptide_dict[acc].add(row["Peptide"])
    for acc, pep_set in acc_to_peptide_dict.items():
        if len(pep_set)<nb_peptides:
            accs_to_remove.add(acc)
        else:
            accs_to_keep.add(acc)
    df_no_one_hits=df[remove_acc_row(df.Protein, accs_to_remove)]
    df_no_one_hits.Protein = df_no_one_hits.Protein.apply(lambda acc_list: set(acc_list).intersection(accs_to_keep))
    df_no_one_hits.Protein = df_no_one_hits.Protein.apply(lambda acc_set: sorted(list(acc_set)))
    df_no_one_hits = df_no_one_hits[remove_empty_rows(df_no_one_hits.Protein)]
    df_no_one_hits[f"taxID_{level}"] = df_no_one_hits.Protein.apply(lambda acc_list: add_tax_information(acc_list, acc2tax_dict, level))

    return df_no_one_hits

def remove_unrelated_accs_and_taxa(df, taxids, level, acc_2_taxid_dict):
    df.Protein = df.Protein.apply(lambda acc_list: [acc for acc in acc_list if \
                                                    get_taxid_of_acc(acc, level, taxon_graph, acc_2_taxid_dict) in taxids])

    df[f"taxID_{level}"] = df.Protein.apply(lambda acc_list: [get_taxid_of_acc(acc, level, taxon_graph, acc_2_taxid_dict) for acc in acc_list])
    return df

def taxon_exclusive_spectra(df, df_no_one_hits, taxon, level, acc_2_taxid_dict):
    df_identified_spectra_for_taxa = df[get_taxa_rows(df[f"taxID_{level}"], taxon)]
    if df_identified_spectra_for_taxa.empty:
        df_identified_spectra_for_taxa = pd.DataFrame(columns=["Title", "Protein", f"taxID_species"])
    df_exclusive_spectra = df_identified_spectra_for_taxa[get_exclusive_rows(df_identified_spectra_for_taxa[f"taxID_{level}"], taxon)]
    if df_exclusive_spectra.empty:
        df_exclusive_spectra = pd.DataFrame(columns=["Title", "Protein", f"taxID_species"])
    exclusive_spectra = set(df_exclusive_spectra.Title)
    df_identified_spectra_for_taxa = remove_unrelated_accs_and_taxa(df_identified_spectra_for_taxa, taxon, level, acc_2_taxid_dict)
    all_identified_accs = set([item for sublist in df_identified_spectra_for_taxa.Protein for item in sublist])
    all_identified_spectra = set(df_identified_spectra_for_taxa.Title)

    df_identified_spectra_no_one_hits_for_taxa = df_no_one_hits[get_taxa_rows(df_no_one_hits[f"taxID_{level}"], taxon)]
    if df_identified_spectra_no_one_hits_for_taxa.empty:
        df_identified_spectra_no_one_hits_for_taxa = pd.DataFrame(columns=["Title", "Protein", f"taxID_species"])
    df_exclusive_spectra_no_one_hits = df_identified_spectra_no_one_hits_for_taxa[get_exclusive_rows(df_identified_spectra_no_one_hits_for_taxa[f"taxID_{level}"], taxon)]
    if df_exclusive_spectra_no_one_hits.empty:
        df_exclusive_spectra_no_one_hits = pd.DataFrame(columns=["Title", "Protein", f"taxID_species"])
    exclusive_spectra_no_one_hits = set(df_exclusive_spectra_no_one_hits.Title)
    df_identified_spectra_no_one_hits_for_taxa = remove_unrelated_accs_and_taxa(df_identified_spectra_no_one_hits_for_taxa, taxon, level, acc_2_taxid_dict)
    all_identified_spectra_no_one_hits = set(df_identified_spectra_no_one_hits_for_taxa.Title)
    all_identified_accs_no_one_hits = set([item for sublist in df_identified_spectra_no_one_hits_for_taxa.Protein for item in sublist])

    return all_identified_spectra, all_identified_spectra_no_one_hits, all_identified_accs, all_identified_accs_no_one_hits, exclusive_spectra, exclusive_spectra_no_one_hits


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

def spectra_identified_in_df(df_species, df_family, taxon_species, taxon_family):
    species_spectra = set(df_species[get_taxa_rows(df_species.taxID_species, taxon_species)].Title)
    family_spectra = set(df_family[get_taxa_rows(df_family.taxID_family, taxon_family)].Title)
    new_identified_spectra = family_spectra.difference(species_spectra)
    all_species_spectra = set(df_species.Title)
    spectra_in_species_df = new_identified_spectra.intersection(all_species_spectra)
    return new_identified_spectra, spectra_in_species_df

def load_spec_and_fam_df():
    fdr=0.05
    df_in_fdr_uniprot_species = get_df_in_fdr_without_decoy(uniprot_nr_reduced_tsv['species'], fdr, columns=['taxID_species'])
    df_in_fdr_uniprot_family =  get_df_in_fdr_without_decoy(uniprot_nr_reduced_tsv['family'], fdr, columns=['taxID_family'])
    df_in_fdr_uniprot_species.Protein = df_in_fdr_uniprot_species.Protein.apply(lambda acc_set: {acc.split('|')[1] \
                                                                                                 for acc in acc_set})
    df_in_fdr_uniprot_family.Protein = df_in_fdr_uniprot_family.Protein.apply(lambda acc_set: {acc.split('|')[1] \
                                                                                               for acc in acc_set})

    return df_in_fdr_uniprot_species, df_in_fdr_uniprot_family

def write_header(out):
    out.write(f"name\ttaxon_species\ttaxon_family\t")
    out.write(f"all identified_spectra for taxa level species\tall identified_spectra for taxa level family\t")
    out.write(f"all identified_spectra for taxa level species no one hits\tall identified_spectra for taxa level family no one hits\t")
    out.write("all_identified_accs_spe\tall_identified_accs_spe_no_one_hits\t")
    out.write("all_identified_accs_fam\tall_identified_accs_fam_no_one_hits\t")
    out.write(f"exclusive spectra for taxa level species\texclusive spectra for taxa level family\t")
    out.write(f"intersection species and family exclusive spectra\t")
    out.write(f"exclusive spectra for taxa level species no one hits\texclusive spectra for taxa level family no one hits\t")
    out.write(f"intersection species and family exclusive spectra no one hits\t")
    out.write(f"new_identified_spectra in fam\tnew_identified_spectra in fam bereits in species_df\t")
    out.write(f"new_identified_spectra in fam no one hits\tnew_identified_spectra in fam bereits in species_df no one hit\n")

def write_results(out, name, taxon_species, taxon_family,
                  all_identified_spectra_spe, all_identified_spectra_fam,
                  all_identified_spectra_no_one_hits_spe, all_identified_spectra_no_one_hits_fam,
                  all_identified_accs_spe, all_identified_accs_spe_no_one_hits,
                  all_identified_accs_fam, all_identified_accs_fam_no_one_hits,
                  exclusive_spectra_spe, exclusive_spectra_fam,
                  exclusive_spectra_no_one_hits_spe, exclusive_spectra_no_one_hits_fam,
                  new_identified_spectra, spectra_in_species_df,
                  new_identified_spectra_no_one_hits, spectra_in_species_df_no_one_hits):
    out.write(f"{name}\t{taxon_species}\t{taxon_family}\t")
    out.write(f"{len(all_identified_spectra_spe)}\t{len(all_identified_spectra_fam)}\t")
    out.write(f"{len(all_identified_spectra_no_one_hits_spe)}\t{len(all_identified_spectra_no_one_hits_fam)}\t")
    out.write(f"{len(all_identified_accs_spe)}\t{len(all_identified_accs_spe_no_one_hits)}\t")
    out.write(f"{len(all_identified_accs_fam)}\t{len(all_identified_accs_fam_no_one_hits)}\t")
    out.write(f"{len(exclusive_spectra_spe)}\t{len(exclusive_spectra_fam)}\t")
    out.write(f"{len(exclusive_spectra_spe.intersection(exclusive_spectra_fam))}\t")
    out.write(f"{len(exclusive_spectra_no_one_hits_spe)}\t{len(exclusive_spectra_no_one_hits_fam)}\t")
    out.write(f"{len(exclusive_spectra_no_one_hits_spe.intersection(exclusive_spectra_no_one_hits_fam))}\t")
    out.write(f"{len(new_identified_spectra)}\t{len(spectra_in_species_df)}\t")
    out.write(f"{len(new_identified_spectra_no_one_hits)}\t{len(spectra_in_species_df_no_one_hits)}\t")
    out.write("\n")



def main():
    df_in_fdr_uniprot_species, df_in_fdr_uniprot_family = load_spec_and_fam_df()
    acc_to_taxid_dict = get_acc2taxid_dict(df_in_fdr_uniprot_species, df_in_fdr_uniprot_family)
    df_in_fdr_uniprot_species = sort_taxid_and_acc_in_df(df_in_fdr_uniprot_species, 'species', acc_to_taxid_dict)
    df_in_fdr_uniprot_family = sort_taxid_and_acc_in_df(df_in_fdr_uniprot_family, 'family', acc_to_taxid_dict)

    df_in_fdr_uniprot_species_no_one_hits = remove_all_accs_with_less_then_x_peptides(df_in_fdr_uniprot_species.copy(deep=True),
                                                                                      'species', 2, acc_to_taxid_dict)
    df_in_fdr_uniprot_family_no_one_hits = remove_all_accs_with_less_then_x_peptides(df_in_fdr_uniprot_family.copy(deep=True),
                                                                                 'family', 2, acc_to_taxid_dict)

    with open("/home/jules/Documents/Metaproteomics/Tax2Proteome/benchmarking/results_analysis/results_per_taxa_more_then_1_hits.tsv", "w") as out:
        write_header(out)

    with open("/home/jules/Documents/Metaproteomics/Tax2Proteome/benchmarking/results_analysis/results_per_taxa_more_then_1_hits.tsv", "a") as out:
        for name, taxon_species in Kleiner_2:
            taxon_family = [taxon_graph.find_level_up(taxon_species[0], "family")]

            all_identified_spectra_spe, all_identified_spectra_no_one_hits_spe, all_identified_accs_spe, all_identified_accs_spe_no_one_hits, exclusive_spectra_spe, exclusive_spectra_no_one_hits_spe = \
                taxon_exclusive_spectra(df_in_fdr_uniprot_species, df_in_fdr_uniprot_species_no_one_hits, taxon_species, 'species', acc_to_taxid_dict)
            all_identified_spectra_fam, all_identified_spectra_no_one_hits_fam, all_identified_accs_fam, all_identified_accs_fam_no_one_hits, exclusive_spectra_fam, exclusive_spectra_no_one_hits_fam = \
                taxon_exclusive_spectra(df_in_fdr_uniprot_family, df_in_fdr_uniprot_family_no_one_hits, taxon_family, 'family', acc_to_taxid_dict)
            new_identified_spectra, spectra_in_species_df = spectra_identified_in_df(df_in_fdr_uniprot_species,
                                                                                     df_in_fdr_uniprot_family, taxon_species, taxon_family)

            new_identified_spectra_no_one_hits, spectra_in_species_df_no_one_hits = spectra_identified_in_df(df_in_fdr_uniprot_species_no_one_hits,
                                                                                                             df_in_fdr_uniprot_family_no_one_hits, taxon_species, taxon_family)
            write_results(out, name, taxon_species, taxon_family,
                          all_identified_spectra_spe, all_identified_spectra_fam,
                          all_identified_spectra_no_one_hits_spe, all_identified_spectra_no_one_hits_fam,
                          all_identified_accs_spe, all_identified_accs_spe_no_one_hits,
                          all_identified_accs_fam, all_identified_accs_fam_no_one_hits,
                          exclusive_spectra_spe, exclusive_spectra_fam,
                          exclusive_spectra_no_one_hits_spe, exclusive_spectra_no_one_hits_fam,
                          new_identified_spectra, spectra_in_species_df,
                          new_identified_spectra_no_one_hits, spectra_in_species_df_no_one_hits)

        for name, taxon_species, taxon_family in Kleiner_3:
            all_identified_spectra_spe, all_identified_spectra_no_one_hits_spe, all_identified_accs_spe, all_identified_accs_spe_no_one_hits, exclusive_spectra_spe, exclusive_spectra_no_one_hits_spe = \
                taxon_exclusive_spectra(df_in_fdr_uniprot_species, df_in_fdr_uniprot_species_no_one_hits, taxon_species, 'species', acc_to_taxid_dict)
            all_identified_spectra_fam, all_identified_spectra_no_one_hits_fam, all_identified_accs_fam, all_identified_accs_fam_no_one_hits, exclusive_spectra_fam, exclusive_spectra_no_one_hits_fam = \
                taxon_exclusive_spectra(df_in_fdr_uniprot_family, df_in_fdr_uniprot_family_no_one_hits, taxon_family, 'family', acc_to_taxid_dict)
            new_identified_spectra, spectra_in_species_df = spectra_identified_in_df(df_in_fdr_uniprot_species,
                                                                                     df_in_fdr_uniprot_family, taxon_species, taxon_family)

            new_identified_spectra_no_one_hits, spectra_in_species_df_no_one_hits = spectra_identified_in_df(df_in_fdr_uniprot_species_no_one_hits,
                                                                                                             df_in_fdr_uniprot_family_no_one_hits, taxon_species, taxon_family)
            write_results(out, name, taxon_species, taxon_family,
                          all_identified_spectra_spe, all_identified_spectra_fam,
                          all_identified_spectra_no_one_hits_spe, all_identified_spectra_no_one_hits_fam,
                          all_identified_accs_spe, all_identified_accs_spe_no_one_hits,
                          all_identified_accs_fam, all_identified_accs_fam_no_one_hits,
                          exclusive_spectra_spe, exclusive_spectra_fam,
                          exclusive_spectra_no_one_hits_spe, exclusive_spectra_no_one_hits_fam,
                          new_identified_spectra, spectra_in_species_df,
                          new_identified_spectra_no_one_hits, spectra_in_species_df_no_one_hits)

