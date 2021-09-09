import pandas as pd


class PSM_FDR:
    def __init__(self, path_to_file, path_to_crap, decoy_tag):
        """
        :param path_to_file: path to xtandem output .xml converted to tsv
        one line one header with all accessions seperated by \t
        """
        self.path_to_file = path_to_file
        self.crap_acc_set = self.get_all_crap_accs(path_to_crap)
        self.decoy_tag = decoy_tag
        self.decoy_list = ['DECOY', 'CRAP', 'DECOY/CRAP']
        self.sorted_xtandem_df = None
        self.fdr_pos = 0
        self.number_psms = 0
        self.decoys = 0

    def flatten_set(self, s):
        return {item for sublist in s for item in sublist}

    @staticmethod
    def get_all_crap_accs(path_to_crap):
        crap_acc_set = set()
        with open(str(path_to_crap)) as crap:
            for line in crap.readlines():
                if line.startswith('>'):
                    crap_acc_set.add(line[1:].strip())
        return crap_acc_set

    def add_level_specific_taxid_column(self, taxon_graph, level):
        self.sorted_xtandem_df[f'taxID_{level}'] = self.sorted_xtandem_df[f'taxID'].apply(lambda taxid:
                                                taxon_graph.find_level_up(taxid, level)
                                                if taxid not in self.decoy_list else taxid)
        return self.sorted_xtandem_df

    def group_df(self, level):
        reduced_df = self.sorted_xtandem_df.groupby(['Title', 'Peptide', 'Hyperscore'], as_index=False).agg(
            {'Protein': lambda acc: set(acc), 'decoy': lambda x: set(x),
             'taxID': lambda taxid: set(taxid), f'taxID_{level}': lambda x: set(x)})
        return reduced_df

    def get_taxid_of_prot_acc(self, protein_acc, acc2tax_dict, db_type):
        if db_type=='ncbi':
            if protein_acc in self.crap_acc_set:
                return {'CRAP'}
            if self.decoy_tag in protein_acc:
                return {'DECOY'}
            try:
                # get taxa set
                return acc2tax_dict[protein_acc]
            except KeyError:
                print('DECOY/CRAP- protein accs: ', protein_acc)
                return {'DECOY/CRAP'}
        else:
            if protein_acc in self.crap_acc_set:
                return 'CRAP'
            if self.decoy_tag in protein_acc:
                return 'DECOY'
            try:
                if db_type == 'uniprot' or db_type == 'swissprot':
                    return int(acc2tax_dict[protein_acc.split('|')[1]])
                elif db_type == 'custom':
                    return int(acc2tax_dict[protein_acc.strip()])
            except KeyError:
                # custom acc
                if protein_acc.startswith('CRAP'):
                    return 'CRAP'
                else:
                    print('DECOY/CRAP- protein accs: ', protein_acc)
                    return 'DECOY/CRAP'

    def create_PSM_dataframe_for_uniprot_accs(self, acc2tax_dict, taxon_graph, level):
        self.sorted_xtandem_df['Protein'] = self.sorted_xtandem_df['Protein'].apply(lambda protein_acc: protein_acc.split()[0])
        self.sorted_xtandem_df['taxID'] = self.sorted_xtandem_df['Protein'].apply(lambda protein_acc:
                                                                self.get_taxid_of_prot_acc(protein_acc, acc2tax_dict, 'uniprot'))

        self.sorted_xtandem_df = self.add_level_specific_taxid_column(taxon_graph, level)
        reduced_df = self.group_df(level)
        return reduced_df

    def get_first_acc(self, acc):
        if not acc.startswith('generic'):
            return acc
        decoy = True if self.decoy_tag in acc else False
        if decoy:
            return acc.split('|', maxsplit=2)[2].split()[0] + '_' + self.decoy_tag
        else:
            return acc.split('|', maxsplit=2)[2].split()[0]

    def get_taxa_set_of_specified_level(self, taxon_graph, taxID_set, level):
        level_up_set=set()
        for taxID in taxID_set:
            if taxID in self.decoy_list:
                level_up_set.add(taxID)
            elif taxID == '0' or taxID == 0: # for some accs in prot2accs file
                continue
            else:
                level_up_set.add(taxon_graph.find_level_up(int(taxID), level))
        return level_up_set

    def create_PSM_dataframe_for_ncbi_accs(self, acc2tax_dict, taxon_graph, level):
        """
        :param acc2tax_dict: for ncbi multiple accs per accs posiible, so multiple taxa posiible: {acc string: {set of taxa(int)}
        :param taxon_graph: loaded TaxonGraph
        :param level: 'species', 'genus', ...
        :return: reduced_df same spectra ID and Hyperscore entries combined into one, all values into sets
        """
        # Multiaccs in Protein column: 'PNW76085.1 BAB64417.1 BAB64413.1 XP_001693987.1'
        print('get first accs in protein column')
        self.sorted_xtandem_df['Protein'] = \
            self.sorted_xtandem_df['Protein'].apply(lambda acc: self.get_first_acc(acc))
        print('get taxIDs of protein accs')
        self.sorted_xtandem_df['taxID'] = self.sorted_xtandem_df['Protein'].apply(lambda protein_acc:
                                                                                  self.get_taxid_of_prot_acc(protein_acc, acc2tax_dict, 'ncbi'))
        # 'Protein': generic|AC:4260012|KKY45073.1 WP_046834534.1, generic|AC:6172351|WP_086660554.1 ...
        # NCBI: 'taxID': {83334}, {562}, takes lot of time
        print('get taxIDs of specified level')
        self.sorted_xtandem_df[f'taxID_{level}'] = \
            self.sorted_xtandem_df['taxID'].apply(lambda taxID_set: self.get_taxa_set_of_specified_level(taxon_graph, taxID_set, level))
        print('reducing df')
        reduced_df = self.sorted_xtandem_df.groupby(['Title', 'Peptide', 'Hyperscore'], as_index=False).agg(
            {'Protein': lambda acc: set(acc), 'decoy': lambda decoy: set(decoy),
             'taxID': lambda taxid_sets: self.flatten_set(taxid_sets),
             f'taxID_{level}': lambda taxid_sets: self.flatten_set(taxid_sets)})
        return reduced_df

    def create_PSM_dataframe_for_custom_accs(self, acc2tax_dict, taxon_graph, level):
        self.sorted_xtandem_df['taxID'] = self.sorted_xtandem_df['Protein'].apply(
            lambda acc: self.get_taxid_of_prot_acc(acc, acc2tax_dict, 'custom'))
        self.sorted_xtandem_df = self.add_level_specific_taxid_column(taxon_graph, level)
        reduced_df = self.group_df(level)
        return reduced_df

    def create_PSM_dataframe(self, db_type, level, taxon_graph, acc2tax_dict):
        """

        :return: Hyperscore sorted dataframe, Protein with highest Hyperscore first, Title = spectra_file_spectraID
        """
        xtandem_df = pd.read_csv(str(self.path_to_file), delimiter='\t')
        xtandem_df['Protein'] = xtandem_df['Protein'].apply(lambda acc: acc.strip())
        # change spectra Title
        xtandem_df['Title'] = xtandem_df['Title'].apply(lambda row: row.split(' File')[0])
        print('Sorting panda dataframe by hyperscore. Find taxa, and level taxa')
        self.sorted_xtandem_df = xtandem_df.sort_values(by=['Hyperscore', 'Title'], ascending=False).reset_index(drop=True)
        self.sorted_xtandem_df['decoy'] = self.sorted_xtandem_df.apply(lambda row: True if self.decoy_tag in row['Protein'] else False, axis=1)
        print('create dataframe with level taxa')
        if db_type == 'uniprot' or db_type == 'swissprot':
            reduced_df = self.create_PSM_dataframe_for_uniprot_accs(acc2tax_dict, taxon_graph, level)
        elif db_type == 'ncbi':
            reduced_df = self.create_PSM_dataframe_for_ncbi_accs(acc2tax_dict, taxon_graph, level)
        elif db_type == 'custom':
            reduced_df = self.create_PSM_dataframe_for_custom_accs(acc2tax_dict, taxon_graph, level)

        # sort by Hyperscore
        reduced_df = reduced_df.sort_values(by=['Hyperscore', 'Title'], ascending=False).reset_index(drop=True)

        print(f'entries in sorted df: {len(self.sorted_xtandem_df)} '
              f'decoys in sorted df: {len(self.sorted_xtandem_df[self.sorted_xtandem_df.decoy==True])} '
              f'hits in sorted df: {len(self.sorted_xtandem_df[self.sorted_xtandem_df.decoy==False])}')
        print(f'entries in reduced_df: {len(reduced_df)} '
              f'decoys in reduced_df: {len(reduced_df[reduced_df.decoy=={True}])} '
              f'hits in reduced_df: {len(reduced_df[reduced_df.decoy=={False}])} '
              f'mixed in reduced_df: {len(reduced_df[reduced_df.decoy=={False, True}])}')
        return reduced_df

    @staticmethod
    def get_all_PSMs(decoy_column):
        return [True if False in decoy_set else False for decoy_set in decoy_column]

    @staticmethod
    def get_all_decoys(decoy_column):
        return [True if True in decoy_set else False for decoy_set in decoy_column]

    @staticmethod
    def determine_FDR_position(sorted_xtandem_df, fdr, is_decoy_column_set, decoy_column_name='decoy'):
        """
        :param fdr: false discovery rate, for example 0.01
        """
        repeatedly_identified_spectra = set()
        number_multiple_identified_spectra = 0
        hits = decoy = 0
        FDR_position_not_set = True
        title_set = set()
        spectra_header = [(x, y) for x, y in zip(sorted_xtandem_df['Title'], sorted_xtandem_df[decoy_column_name])]
        fdr_pos = len(sorted_xtandem_df)-1
        is_ref_file = True if 'Ref_decoy' in sorted_xtandem_df.columns else False
        if is_ref_file:
            number_psms = len(set(sorted_xtandem_df[PSM_FDR.get_all_PSMs(sorted_xtandem_df.Ref_decoy)].Title))
            decoys = len(set(sorted_xtandem_df[PSM_FDR.get_all_decoys(sorted_xtandem_df.Ref_decoy)].Title))
        else:
            number_psms = len(set(sorted_xtandem_df[PSM_FDR.get_all_PSMs(sorted_xtandem_df.decoy)].Title))
            decoys = len(set(sorted_xtandem_df[PSM_FDR.get_all_decoys(sorted_xtandem_df.decoy)].Title))

        for elem in spectra_header:
            if elem[0] in title_set:
                repeatedly_identified_spectra.add(elem[0])
                number_multiple_identified_spectra += 1
                continue
            else:
                title_set.add(elem[0])
            if not is_decoy_column_set:
                if elem[1]:
                    decoy += 1
                else:
                    hits += 1
            else:
                ## TODO or elem[1] == {TRUE} ? what to do with mixed
                if True in elem[1]:
                    decoy += 1
                else:
                    hits += 1
            if decoy / (hits + decoy) > fdr and FDR_position_not_set:
                fdr_pos = hits + decoy - 1 + number_multiple_identified_spectra
                number_psms = hits
                decoys = decoy - 1
                break
            if decoy / (hits + decoy) <= fdr:
                continue
        if is_ref_file:
            score_last_item = sorted_xtandem_df['Ref_Hyperscore'][fdr_pos]
        else:
            score_last_item = sorted_xtandem_df['Hyperscore'][fdr_pos]

        # repeatedly_identified_spectra of reduced_df: different Proteins, same spectra,
        # print('Number of PSMs: %d' % number_psms)
        # print('Number of decoys: %d' % decoys)
        # print(f"double identified spectra {number_multiple_identified_spectra}")
        # print('Position FDR border/Number of PSMs: %d' % fdr_pos)
        # print('score last item: %d' % score_last_item)
        return fdr_pos, number_psms, decoys, number_multiple_identified_spectra, score_last_item


