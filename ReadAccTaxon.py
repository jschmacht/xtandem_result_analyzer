from pathlib import Path
import multiprocessing as mp
import gzip


class ReadAccTaxon:
    def __init__(self, path_to_folder, db_type):
        """
        :param path_to_folder: path to folder with accession2taxonID files
        (Names NCBI prot.accession2taxid & pdb.accession2taxid.gz Uniprot: generated file from uniprot database: acc2tax_uniprot
        one line one header with all accessions seperated by \t
        """

        self.path_to_folder = Path(path_to_folder)
        self.db_type = db_type
        self.acc_taxon_dict = {}

    @staticmethod
    def get_acc2taxonID_dict(path_to_db):
        """
        preotein accession to taxon ID from custom db file
        :param path_to_db:
        :return:
        """
        acc_taxon_dict={}
        with open(path_to_db) as acc2tax:
            for line in acc2tax:
                fields = line.split()
                acc_taxon_dict[fields[1]] = int(fields[-1])
        return acc_taxon_dict

    # return generators, chunks of complete lines
    def read_in_chunks(self, file, size=1024 * 1024):
        f = open(str(file), "rb")
        f.readline()
        while True:
            start = f.tell()
            f.seek(size + start)
            s = f.readline()
            end = (f.tell() - start)
            f.seek(start + end)
            yield start, end
            if not s:
                f.close()
                break


    def read_pdb2accession(self):
        with gzip.open(str(self.path_to_folder / 'pdb.accession2taxid.gz'), "rt") as accession_handle:
            accession_handle.readline()
            part_dict = {}
            for line in accession_handle:
                fields = [item.strip('\t') for item in line.split('\t')]
                if fields[1] in accs:
                    part_dict[fields[1]] = fields[2]
        return part_dict


    def read_chunks_into_part_dict_ncbi(self, chunk):
        """
                :param chunk: tuple of numbers, chunk[0]=entry point, chunk[1]= size to be readed
                  """
        try:
            f = open(str(self.path_to_folder / 'prot.accession2taxid'), "rb")
            f.seek(chunk[0])
            part_dict = {}
            for line in f.read(chunk[1]).decode("utf-8").splitlines():
                fields = [item.strip('\t') for item in line.split('\t')]
                if fields[1] in accs:
                    part_dict[fields[1]] = fields[2]
            f.close()
            return part_dict
        except FileNotFoundError:
            print("Path to database does not exist.")
            exit(1)


    def read_chunks_into_part_dict_uniprot(self, chunk):
        """
                :param chunk: tuple of numbers, chunk[0]=entry point, chunk[1]= size to be readed
                  """
        try:
            f = open(str(self.path_to_folder / 'acc2tax_uniprot'), "rb")
            f.seek(chunk[0])
            part_dict = {}
            for line in f.read(chunk[1]).decode("utf-8").splitlines():
                fields = line.split('\t')
                if fields[0] in accs:
                    part_dict[fields[0]] = fields[1]
            f.close()
            return part_dict
        except FileNotFoundError:
            print("Path to database does not exist.")
            exit(1)

    # read uniprot_acc2tax or NCBI accession2taxid file and put matches in acc_dict with accession as key
    def read_acc2tax(self, accessions, threads=None):
        acc_taxon_dict={}
        global accs
        accs = accessions

        if not threads:
            threads = mp.cpu_count()
        pool = mp.Pool(threads)
        # chunks = list of all chunks list tuples: (entry point to db, number of chars to be read)
        chunks = []
        if self.db_type == 'uniprot':
            for chunk in self.read_in_chunks(self.path_to_folder / 'acc2tax_uniprot'):
                chunks.append(chunk)
        elif self.db_type == 'ncbi':
            for chunk in self.read_in_chunks(self.path_to_folder / 'prot.accession2taxid'):
                chunks.append(chunk)
        print('Start reading accession2prot database file with %s threads.' % str(threads))
        # here multiprocessing starts, results saved in self.accessionIDs
        i, j = 0, 0
        ten = int(len(chunks) / 10)
        if self.db_type == 'uniprot':
            for part_dict in pool.imap(self.read_chunks_into_part_dict_uniprot, chunks):
                if i == ten:
                    j += 1
                    print('%d0%% readed.' % j)
                    i = -1
                acc_taxon_dict.update(part_dict)
                i += 1
        elif self.db_type == 'ncbi':
            for part_dict in pool.imap(self.read_chunks_into_part_dict_ncbi, chunks):
                if i == ten:
                    j += 1
                    print('%d0%% readed.' % j)
                    i = -1
                acc_taxon_dict.update(part_dict)
                i += 1
            acc_taxon_dict.update(self.read_pdb2accession())
        pool.close()
        self.acc_taxon_dict = acc_taxon_dict
        return acc_taxon_dict
