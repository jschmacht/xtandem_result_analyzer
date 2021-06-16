from pathlib import Path
from ReadAccTaxon import ReadAccTaxon
from SearchAccessions import SearchAccessions

class Multiaccs():
    def __init__(self, path_to_multiacc_file, path_to_multiaccs_taxids_file):
        self.path_to_multiaccs_taxids_file = Path(path_to_multiaccs_taxids_file)
        if not self.path_to_multiaccs_taxids_file.exists():
            self.write_multaccs_to_taxids_file(path_to_multiacc_file)
        self.multiacc_marks = ['WP', 'XP']
        self.multiaccs_to_taxid_dict = self.read_multiaccs_file()

    def read_multiaccs_file(self):
        multiaccs_to_taxid_dict = {}
        with open(self.path_to_multiaccs_taxids_file) as multiacc_tax:
            for line in multiacc_tax:
                fields = line.split()
                multiaccs_to_taxid_dict[fields[0]] = fields[1:]
        return multiaccs_to_taxid_dict

    def is_multiacc(self, acc):
        if any(self.multiacc_marks) in acc or '...' in acc:
            return True
        else:
            return False

    def get_multiacc(self, acc):
        accs = acc.split()
        acc_marks = [acc for acc in accs if acc.split('_')[0] in self.multiacc_marks]
        if len(acc_marks) == 1:
            return acc_marks[0]
        else:
            print(accs)

    def get_taxids_of_multiaccs(self, acc):
        taxid_list = self.multiaccs_to_taxid_dict(acc)
        return taxid_list

    def write_multaccs_to_taxids_file(self, path_to_multiacc_file):
        accessions = []
        all_marks = set()
        with open(path_to_multiacc_file, 'r') as multiacc:
            for line in multiacc:
                accessions = accessions.extend(line.split())
        accessions = set(accessions)
        path_tax2proteome_db = '/home/jules/Documents/databases/databases_tax2proteome/'
        acc2taxon = ReadAccTaxon(path_tax2proteome_db, 'ncbi')
        acc_taxon_dict = acc2taxon.read_acc2tax(accessions)
        with open(path_to_multiacc_file, 'r') as multiacc:
            with open(self.path_to_multiaccs_taxids_file, 'r') as multiacc_tax:
                for line in multiacc:
                    fields = line.split()
                    multiacc = fields[0]
                    all_marks.add(multiacc.split('_')[0])
                    taxa = {acc_taxon_dict[acc] for acc in fields}
                    accessions = accessions.union(set(line.split()))
                    multiacc_tax.write(multiacc + '\t' + ('\t').join(accessions))
        print('all marks', all_marks)