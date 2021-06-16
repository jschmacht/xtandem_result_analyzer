import gzip
import re
import logging
from pathlib import Path
import multiprocessing as mp


class SearchAccessions:
    def __init__(self, path_to_db, db):
        """
        :param path_to_xml: path to the protein database
        :param db: database format, 'ncbi' or 'uniprot' )
        """
        self.path_to_file = Path(path_to_db)
        self.db = db
        self.acc_dict = {}
        self.accs = []

    # ncbi peptide database
    def read_ncbi_mp(self, chunk):
        """
        :param chunk: tuple of numbers, chunk[0]=entry point in database file, chunk[1]= size to be readed
        """
        f = open(str(self.path_to_file), "rb")
        f.seek(chunk[0])
        part_dict={}
        for line in f.read(chunk[1]).decode("utf-8").splitlines():
            fields = [item.strip('\t') for item in line.split('\t')]
            if fields[0] in accessionIDs:
                part_dict[fields[0]]= fields
        return part_dict

    def read_nr_mp(self, chunk):
        """
        :param chunk: tuple of numbers, chunk[0]=entry point in database file, chunk[1]= size to be readed
        """
        f = open(str(self.path_to_file), "rb")
        f.seek(chunk[0])
        part_dict={}
        for line in f.read(chunk[1]).decode("utf-8").splitlines():
            if line.startswith('>'):
                if '\x01' in line:
                    line = line[1:]
                    accs = [item.split(' ')[0] for item in line.split('\x01')]
                    part_dict[accs[0]] = accs
        return part_dict

    def read_in_chunks(self, size=1024 * 1024):
        f = open(str(self.path_to_file), "rb")
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

    # divide database file in chunks for parallel processing
    def divide_into_chunks(self):
        chunks = []
        print('Start dividing database file in chunks.')
        for chunk in self.read_in_chunks():
            chunks.append(chunk)
        print('Database divided into %d chunks.' % len(chunks))
        return chunks

    # read database and store position information in list (value) to taxon_ID
    # if database is without taxon IDs, before matching scientifc name to ID
    def read_database(self, accessions, threads=None):
        """
         :param accessions: set of accessionIDs matching searched taxon IDs (for ncbi db)
         """
        if not threads:
            threads = mp.cpu_count()
        global accessionIDs
        accessionIDs = set()
        for i in accessions:
            accessionIDs.update(set(i))
        print('Set accessions length %d' % len(accessionIDs))
        chunks = self.divide_into_chunks()
        pool = mp.Pool(threads)
        # here multiprocessing starts
        i, j = 0, 0
        ten = int(len(chunks) / 10)
        if self.db == 'ncbi':
            for part_dict in pool.imap_unordered(self.read_ncbi_mp, chunks):
                if i == ten:
                    j += 1
                    print('%d0%% readed.' % j)
                    i = -1
                self.acc_dict.update(part_dict)
                i += 1
        if self.db == 'uniprot':
            for part_dict in pool.imap_unordered(self.read_nr_mp, chunks):
                if i == ten:
                    j += 1
                    print('%d0%% readed.' % j)
                    i = 0
                self.acc_dict.update(part_dict)
                i += 1
        pool.close()
        pool.join()
