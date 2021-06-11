import tarfile
import logging.config
from collections import defaultdict


class TaxonGraph():
    def __init__(self):
        self.parent_child_graph = defaultdict(list)
        self.child_rank_dict = {}
        self.oldtax_newtax_dict = {}
        self.taxon_name_dict = {}
        self.scientific_alternative_name_dict = {}
        self.child_parent_dict = {}
        self.synonym_taxon_dict = {}
        self.order = {'no rank': 0, 'varietas': 1, 'forma': 1, 'subspecies': 2, 'strain':2, 'species': 3, 'species subgroup': 4,
                 'species group': 5, 'series': 6, 'subsection': 7, 'section': 8, 'subgenus': 9, 'genus': 10,
                 'subtribe': 11, 'tribe': 12, 'subfamily': 13, 'family': 14, 'superfamily': 15, 'parvorder': 16,
                 'infraorder': 17, 'suborder': 18, 'order': 19, 'superorder': 20, 'subcohort': 21, 'cohort': 22,
                 'infraclass': 23, 'subclass': 24, 'class': 25, 'superclass': 26, 'subphylum': 27, 'phylum': 28,
                 'superphylum': 29, 'subkingdom': 30, 'kingdom': 31, 'superkingdom': 32}

    # create taxon graph from nodes.dmp (in taxdump.tar.gz archive), a dictionary with taxon-ID-species name 
    # and a dictionary with changed/merged taxon IDs
    def create_graph(self, path_to_taxdump):
        """
        :param path_to_taxdump: path to location of the taxdump.tar.gz file (contains information of ancestry tree, source:NCBI)
          """
        # all lines split at | take 1. as parent, 0. as child and 3. as rank
        try:
            tar = tarfile.open(str(path_to_taxdump), "r:gz")
            nodes_handle = tar.extractfile("nodes.dmp")
            names_handle = tar.extractfile("names.dmp")
            merged_handle = tar.extractfile("merged.dmp")

            for line in nodes_handle:
                fields = [item.strip('\t') for item in line.decode(tarfile.ENCODING).split('|')]
                self.child_rank_dict[int(fields[0])] = fields[2]
                self.child_parent_dict[int(fields[0])] = int(fields[1])
            for line in merged_handle:
                fields = [item.strip('\t') for item in line.decode(tarfile.ENCODING).split('|')]
                self.oldtax_newtax_dict[int(fields[0])] = int(fields[1])
            for line in names_handle:
                fields = [item.strip('\t') for item in line.decode(tarfile.ENCODING).split('|')]
                if fields[3] == 'scientific name':
                    self.taxon_name_dict[int(fields[0])] = fields[1]
                else:
                    self.synonym_taxon_dict[fields[1]] = int(fields[0])

            for key, value in self.child_parent_dict.items():
                self.parent_child_graph[value].append(key)

            # delete circles in graph for top node
            for child, parent in self.child_parent_dict.items():
                if child == parent:
                    self.parent_child_graph[parent].remove(parent)
            # add leaves
            leaves = set(list(self.child_parent_dict.keys()) + list(self.child_parent_dict.values())) - set(self.parent_child_graph.keys())
            self.parent_child_graph.update(dict.fromkeys(leaves, []))

        except tarfile.ReadError as e:
            print("Can not read taxdump.tar.gz File. Without taxon tree, program quits.", exc_info=True)
            exit(1)

    # find for taxID all child tax IDs and returns them
    def find_taxIDs (self, taxID, childs_until_species=False):
        """
         :param taxID: tax ID for searching all children
         :type taxID:  positive int
         :return: set of all child taxon IDs of one taxon ID
         """
        if self.parent_child_graph == {}:
            raise Exception("The parent_child_graph is empty. Check if file is correct and if first function "
                            "\'create_graph\' is called.")
        species_taxIDs = set()
        # check if taxon IDs exist and have valid values
        if taxID not in self.child_rank_dict.keys():
            if self.oldtax_newtax_dict.get(taxID) is None:
                print(
                    'User given taxon ID %d does not exist and is excluded from further analysis.' % taxID)
                return species_taxIDs
            taxID = self.oldtax_newtax_dict.get(taxID)

        if not childs_until_species:
            temp_children = {taxID}
            while temp_children:
                taxID = temp_children.pop()
                species_taxIDs.add(taxID)
                temp_children = temp_children.union(set(self.parent_child_graph.get(taxID)))
        else:
            # if rank of taxon ID is above species level, return same taxID, else go up until level 'species'
            taxID = self.find_level_up(taxID, 'species')
            temp_children = {taxID}
            while temp_children:
                taxID=temp_children.pop()
                if self.child_rank_dict[taxID] == 'species':
                    species_taxIDs.add(taxID)
                else:
                    species_taxIDs.add(taxID)
                    temp_children = temp_children.union(set(self.parent_child_graph.get(taxID)))

        return species_taxIDs

    # find for every taxID the given level up, if level not exist (this is regular) take the next smaller one
    # specified in order. If level of given taxID is higher as searched level, give taxID back unchanged
    def find_level_up(self, taxID, level):
        """
            :param taxID: tax ID for searching all children , positive int
            :param level: user given level for taxon search: string, one out of the strings in order dict
            :return taxon ID of specified level
                 """
        #check if level is really up, else take original taxID
        taxIDnew = taxID
        while self.child_rank_dict[taxIDnew] == 'no rank':
            taxIDnew = self.child_parent_dict.get(taxIDnew)
        # start taxID level higher, returned unchanged taxID
        if self.order[self.child_rank_dict[taxIDnew]] > self.order[level]:
            return taxID
        # start taxID level equals
        elif self.order[self.child_rank_dict[taxIDnew]] == self.order[level]:
            return taxIDnew
        # start taxID level is lower as wished
        else:
            last_specified_level_taxID = None
            while self.order[self.child_rank_dict[taxIDnew]] < self.order[level]:
                if self.child_rank_dict[taxIDnew] != 'no rank':
                    last_specified_level_taxID = (self.child_rank_dict[taxIDnew], taxIDnew)
                taxIDnew = self.child_parent_dict.get(taxIDnew)
            if self.order[self.child_rank_dict[taxIDnew]] == self.order[level]:
                return taxIDnew

            # sometimes not all pylogenetic ranks are filled (101570(=species) -> 196894(=no rank) -> 10699(=family))
            # rank genus and subfamily not present, take next highest with unique classification
            if self.order[self.child_rank_dict[taxIDnew]] > self.order[level]:
                if last_specified_level_taxID:
                    print("No ancestor with level %s of taxID %d exists. TaxID %d of level %s is "
                                   "returned." % (level, taxID, last_specified_level_taxID[1],
                                                  self.child_rank_dict[last_specified_level_taxID[1]]))
                    return last_specified_level_taxID[1]
                else:
                    print("Taxon ID %d has higher level than specified. Unchanged taxon ID %d of level %s is "
                                   "returned." % (taxID, taxID, self.child_rank_dict[taxID]))
                    return taxID

    def get_all_taxids(self, taxonIDs, level):
        taxids_of_level=[]
        all_taxids = []
        for taxon in taxonIDs:
            taxids_of_level.append(self.find_level_up(taxon, level))
        for taxon in taxids_of_level:
            all_taxids.append(self.find_taxIDs(taxon))
        return all_taxids

    # finds the next common ancestror of a set of taxon IDs
    def find_next_common_ancestor(self, taxonIDs):
        """
        :param taxonIDs: set or list of taxon IDs
        :return common ancestor of all taxon IDs
                 """
        # 1. find rank
        ranks=[]
        for taxon in taxonIDs:
            taxIDnew = taxon
            while self.child_rank_dict[taxIDnew] == 'no rank':
                taxIDnew = self.child_parent_dict.get(taxIDnew)
            ranks.append(self.order[self.child_rank_dict[taxIDnew]])
        index = ranks.index(max(ranks))
        highest_taxon = taxonIDs[index]
        del ranks[index]
        del taxonIDs[index]

        #list all parents:
        ancestor_list = [highest_taxon]
        while highest_taxon != 1:
            highest_taxon = self.child_parent_dict[highest_taxon]
            ancestor_list.append(highest_taxon)

        while taxonIDs:
            taxon = taxonIDs.pop()
            while taxon not in ancestor_list:
                taxon = self.child_parent_dict[taxon]
            # delete all list elem before found common ancestor
            index = ancestor_list.index(taxon)
            ancestor_list=ancestor_list[index:]

        return taxon