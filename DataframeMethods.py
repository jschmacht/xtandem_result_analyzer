

def get_hit_rows2(decoy_column):
    return [True if d_set in  [{True, False}, {False}] else False for d_set in decoy_column]


def get_all_exclusive_rows(level_taxa_column):
    l=[]
    for t_set in level_taxa_column:
        if type(t_set) == list:
            t_set=set(t_set)
        if 'CRAP' in t_set:
            t_set.remove('CRAP')
        if 'DECOY' in t_set:
            t_set.remove('DECOY')
        if len(t_set)== 1:
            l.append(True)
        else:
            l.append(False)
    return l


def get_all_non_exclusive_rows(level_taxa_column):
    l=[]
    for t_set in level_taxa_column:
        if type(t_set) == list:
            t_set=set(t_set)
        if 'CRAP' in t_set:
            t_set.remove('CRAP')
        if 'DECOY' in t_set:
            t_set.remove('DECOY')
        if len(t_set) > 1:
            l.append(True)
        else:
            l.append(False)
    return l


def get_taxa_rows(column, taxID):
    if type(taxID)==int:
        return [True if taxID in t_set else False for t_set in column]
    elif type(taxID)==list:
        return [True if len(set(taxID).intersection(set(taxa_set)))>0 else False for taxa_set in column]


def get_taxon_rows(column, taxon):
    return [True if taxon in t_list else False for t_list in column]

def remove_acc_row(column, acc_to_remove_set):
    return [False if len(set(accs).difference(acc_to_remove_set))==0 else True for accs in column]


def remove_empty_rows(column):
    return [False if len(accs)==0 else True for accs in column]

remove_empty_rows, remove_acc_row, get_taxon_rows, get_taxa_rows, get_all_non_exclusive_rows, get_all_exclusive_rows, get_hit_rows2